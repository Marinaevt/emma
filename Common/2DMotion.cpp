#include "StdAfx.h"
#include "2DMotion.h"

IOIMPL (CIO2DMotion, C2DMOTION)

// По времени находим пройденное расстояние (скорости известны)
DBL CIO2DMotion::TimeToLength(size_t i, DBL t)
{
	//TODO: работа с контурами (пока 0)
	//int curveNum = GetContour(0)->GetCurveIndexByIndex(i); //dynamic_cast<C2DContour &>(*(m_contours()[0]))[i];
	
	C2DCurve* curve = GetContour(0)->GetCurve(i);	//dynamic_cast<C2DCurve *>(m_curves()[curveNum]);

	C2DMotionParams& p0 = m_params()[curve->GetStart()];//,
					 //p1 = m_params()[curve->GetEnd()];

	C2DFunction* vel = dynamic_cast<C2DFunction *>(m_vels()[0]); //(C2DFunction&)*(IO::CInterface*)m_vel;
	
	DBL dLength = vel->GetIntegral(p0.GetExitTime(), t);	// L = Integral_t1_t2 (V) dt

	return dLength;
}

bool CIO2DMotion::GetPosOnCurve(size_t i, DBL dTime, C2DPosition& pos)
{
	C2DContour* contour = GetContour(0);
	if (!contour)
		return false;

	//TODO: разобраться с контурами
	C2DCurve* curve = contour->GetCurve(i);
	if (!curve)
		return false;

	Math::C2DPoint tau;	//нормализованный вектор (от начальной до конечной точки)

	bool ret = curve->GetPoint(TimeToLength(i, dTime), pos.m_pos, tau); // .GetPoint(m_nodes(), TimeToLength(i, t), pos.pos, tau);

	if (!ret)
		return false;

	// точка посчитана, нужно рассчитать остальные параметры
	
	C2DMotionParams& p0 = m_params()[curve->GetStart()],
					 p1 = m_params()[curve->GetEnd()];

	C2DFunction* vel = dynamic_cast<C2DFunction *>(m_vels()[0]); //линейная скорость
	C2DFunction* angVel = dynamic_cast<C2DFunction *>(m_vels()[1]); //угловая скорость

	//Нужно ли?
	//DBL alpha = Math::iABto01(p0.GetExitTime(), p1.GetEntrTime(), t); // время, интерполированное в отрезок [0,1]

	pos.m_angvel = angVel->GetValue(dTime); //Math::i01toAB(p0.exitAngvel, p1.entrAngvel, alpha);

	pos.m_angle = m_params()[i].m_angle + angVel->GetIntegral(p0.GetExitTime(), dTime);

	// Length * DBL Time
	DBL dSpeed = vel->GetValue(dTime);
	pos.m_vel = tau * dSpeed; //Math::i01toAB(p0.exitVel, p1.entrVel, alpha);

	return true;
}

bool CIO2DMotion::GetPos(DBL time, C2DPosition& pos)
{
	// эта функция подразумевает, что:
	// - число узлов соответствует числу кривых;
	// - времена в узлах идут в порядке увеличения с номером узла.
	//
	
	//Проверка кол-ва параметров (время, пауза, угол)
	if (m_params().size() < 1)
		return false;

	//Проверка диапазона и указанного времени
	if ((time < m_params()[0].GetEntrTime()) || (time > m_params()[m_params().size() - 1].GetExitTime()))
		return false;

	//Проходим по всем параметрам, ищем нужное время, вызываем метод расчёта позиции 
	for (size_t i = 1; i < m_params().size(); i++)
	{
		if (time <= m_params()[i].GetEntrTime())
		{
			if (GetPosOnCurve(i - 1, time, pos)) 
				return true;
			else
				return false;
		}
		
		if (time <= m_params()[i].GetExitTime())
		{
			pos.m_pos = dynamic_cast<C2DNode *>(m_nodes()[i])->GetPoint();
			pos.m_angvel = 0;
			pos.m_vel = Math::C2DPoint::Zero;
			return true;
		}
	}
	

	return false;
}


bool CIO2DMotion::InitParams(size_t num){

	m_params().resize(num);
	C2DPieceLinearFunction* vel = dynamic_cast<C2DPieceLinearFunction *>(m_vels()[0]); //0 - скорость движения, 1 - скорость вращения (угловая)
	vel->m_args().resize(num);
	vel->m_values().resize(num);

	return true;
}

//! Проверяет параметры траектории
bool CIO2DMotion::CheckParams(){

	if(m_params().size() != m_nodes().size()) return false;

	for (size_t i = 1; i < m_params().size(); i++)
	{
		if (m_params()[i].GetEntrTime() < m_params()[i - 1].GetExitTime()) {
			return false;
		}
		if (m_params()[i].GetExitTime() <= m_params()[i].GetEntrTime()) {
			return false;
		}
	}

	return false;

}

//! Считает время для каждой точки траектории
void CIO2DMotion::CalcTime(){

	C2DPieceLinearFunction* vel = dynamic_cast<C2DPieceLinearFunction *>(m_vels()[0]);	// 0 - линейная скорость
	
	for (size_t i = 1; i < GetNodeCount(); i++){

		//Находим расстояние между соседними точками
		double dLength = (GetNode(i)->GetPoint() - GetNode(i-1)->GetPoint()).Norm();

		//Находим аргумент (время) с индексом i по формуле:
		// t_i = t_i-1 + (2*L)/(V_i + V_i-1); расстояние, делённое на среднюю скорость
		double dTime_i = vel->m_args().at(i-1) + 2*dLength/(vel->m_values()[i-1] + vel->m_values()[i]);
		vel->m_args().at(i) = dTime_i;

		//Устанавливаем вычисленное время в массив параметров C2DMotionParams
		m_params()[i].SetEntrTime(vel->m_args().at(i));
		m_params()[i].SetExitTime(m_params()[i].GetEntrTime());

	}
}

bool CIO2DMotion::GetVel(DBL dLen, DBL &dVelx, DBL &dVely){
	
	//TODO: проблема с контурами до сих пор
	C2DContour* pContour = GetContour(0);
	if (!pContour)
		return false;

	DBL dContourLength = pContour->CalcLength();

	if(dLen < EPS || dLen > dContourLength) {
		return false;
	}

	DBL dTempLen = 0.0;	//суммируем все длины кривых в контуре сюда

	if (pContour->GetCurveCount() == 1){
		return GetVelOnCurve(0, dLen, dVelx, dVely);	//одна кривая - индекс 0
	}

	for(size_t i=0; i<pContour->GetCurveCount(); i++){
		C2DCurve *pcurve = pContour->GetCurve(i);
		dTempLen += pcurve->CalcLength();
		
		if (dLen <= dTempLen){
			return GetVelOnCurve(i, dLen-dTempLen, dVelx, dVely);
		}
	}

	return false;
}

bool CIO2DMotion::GetVelOnCurve(size_t i, DBL dLen, DBL &dVelx, DBL &dVely){
	
	//TODO: проблема с контурами до сих пор
	C2DContour* pContour = GetContour(0);
	if (!pContour)
		return false;

	C2DCurve *pCurve = pContour->GetCurve(i);
	size_t nStart = pCurve->GetStart();
	size_t nEnd = pCurve->GetEnd();	//для определения времени

	//Подразумеваем, что время известно
	DBL t1 = m_params()[nStart].GetExitTime();	//выход, т.к. возможна пауза
	DBL t2 = m_params()[nEnd].GetEntrTime();	//вход, т.к. дальше другая кривая
	if(t1 > t2) {	//если не правильная ориентация кривой
		t1 = m_params()[nEnd].GetExitTime();
		t2 = m_params()[nStart].GetEntrTime();
	}
	DBL dTime = t2-t1;

	DBL dStep = pCurve->CalcLength();
	DBL dCurveLen = dStep;

	while(dStep >= dCurveLen / 10){		//TODO: Magick Number 10
		if(LENGTHUNIT < 1.0){			//TODO: зависит от константы, как зададим
			dStep = dStep*LENGTHUNIT;	
			dTime = dTime*LENGTHUNIT;	// т.к. L = ((V0 + V1)*t)/2, значит dStep & dTime пропорционаяльно изменяется
		}else{
			dStep = dStep/LENGTHUNIT;	
			dTime = dTime/LENGTHUNIT;	// т.к. L = ((V0 + V1)*t)/2, значит пропорционально изменяется
		}
		// V0 -> V1, а V1 -> V0 (сходятся в точке)
		// EPS использовать плохо, т.к. деление DBL\DBL ниже
	}

	Math::C2DPoint pPoint1, pPoint2, pTau;
	//Отступаем слева и справа от точки, чтобы отрезок между оставался равным dStep

	//!!! проверка на выезд за пределы длин +\-
	DBL dMinus = dLen - dStep/2;
	DBL dPlus = dLen + dStep/2;

	if(dMinus < 0.0) {
		dMinus = 0.0;
	}
	if(dPlus > dCurveLen){
		dPlus = dCurveLen;
	}


	if (pCurve->GetPoint(dMinus, pPoint1, pTau) && pCurve->GetPoint(dPlus, pPoint2, pTau)){
		dVelx = fabs(pPoint1.x - pPoint2.x)/dTime;
		dVely = fabs(pPoint1.y - pPoint2.y)/dTime;
		return true;
	}

	return false;
}

bool CIO2DMotion::GetVelOnCurveInPoint(size_t i, DBL dTime, DBL &dVelx, DBL &dVely){
	
	C2DPosition pos1, pos2;
	DBL dt = dTime + EPS;

	C2DCurve *pCurve = GetCurve(i);
	
	bool bFlag1 = GetPosOnCurve(i, dTime, pos1);	//Основная позиция
	C2DMotionParams& pEnd = m_params()[pCurve->GetEnd()];
	if (dt > pEnd.GetEntrTime()){
		dt = dTime - EPS;
	}
	bool bFlag2 = GetPosOnCurve(i, dt, pos2);	//если не вышел за пределы, прогноз

	pos1.m_vel = (pos2.m_pos - pos1.m_pos) / EPS;

	return false; //TODO:

}

// задаёт всем узлам линейную скорость
void CIO2DMotion::SetLinVelToAll(DBL dNewVel){
	
	C2DPieceLinearFunction* vel = dynamic_cast<C2DPieceLinearFunction *>(m_vels()[0]);	// 0 - линейная скорость

	for (size_t i = 1; i < GetNodeCount(); i++){
		vel->m_values().at(i) = dNewVel;
	}
}

/*
void Test4C2DMotion()
{
	using namespace Math;

	//LOGGER.Init(CString(_T("..\\..\\Logs\\C2DMotion.cpp_Test4C2DMotion.txt")));
	//LOGGER.SetTypeStamp(false);
	//DLOG(CString(_T("Testing C2DMotion")), log_info);

	CIO2DMotion mot;
	// создаем узлы
	int n1 = mot.CreateNode(C2DNODE, C2DPoint(0, 0));
	int n2 = mot.CreateNode(C2DNODE, C2DPoint(0, 10));
	int n3 = mot.CreateNode(C2DNODE, C2DPoint(10, 10));

	// и кривые между ними
	int c1 = mot.CreateCurve(C2DLINESEG, n1, n2);
	//int c2 = mot.CreateCurve(C2DLINESEG, n2, n3);
	int c2 = mot.CreateCurve(C2DCIRCARC, n2, n3);
	C2DCircleArc* curve2 = dynamic_cast<C2DCircleArc*>(mot.GetCurve(c2));
	curve2->SetRadius(-10);
	curve2->SetCW(false);


	// создаем контур
	int con1 = mot.CreateContour(C2DCONTOUR);

	// и добавляем в него эти кривые
	mot.GetContour(con1)->UnClose();
	mot.AddCurveToContour(c1, con1);
	mot.AddCurveToContour(c2, con1);

	// создаем (кусочно-)линейную функцию с двумя точками с запасом по времени
	// (0,1) - (100,1)
	C2DPieceLinearFunction fun1;
	fun1.m_args().Add(0);
	fun1.m_args().Add(100);
	fun1.m_values().Add(1);
	fun1.m_values().Add(1);
	// добавляем эту функцию в объект "движение"
	mot.m_vels().Add(dynamic_cast<IO::CInterface *> (&fun1));
	mot.m_vels().Add(dynamic_cast<IO::CInterface *> (&fun1));

	C2DMotionParams par1;
	
	par1.SetEntrTime(0);
	par1.SetPause(0);
	par1.SetAngle(0);
	mot.m_params().Add(par1);

	DBL t1 = 5; //2 * mot.GetCurve(c1)->CalcLength() / (fun1.GetValue(5) - fun1.GetValue(0));
	par1.SetEntrTime(t1);
	par1.SetPause(0);
	par1.SetAngle(0);
	mot.m_params().Add(par1);

	DBL t2 = 10; //2 * mot.GetCurve(c2)->CalcLength() / (fun1.GetValue(10) - fun1.GetValue(5)) + t1;
	par1.SetEntrTime(t2);
	par1.SetPause(0);
	par1.SetAngle(0);
	mot.m_params().Add(par1);

	C2DPosition pos;
	for (double di = 0; di <= t1; di += 0.1)
	{
		if (mot.GetPos(di, pos))
		{
			/*
			DLOG(   CString(_T("T:")) + AllToString(di) +
				CString(_T(" POS:")) + AllToString(pos.m_pos.m_x) + 
				    CString(_T(",")) + AllToString(pos.m_pos.m_y) +
				CString(_T(" VEL:")) + AllToString(pos.m_vel.m_x) +
				    CString(_T(",")) + AllToString(pos.m_vel.m_y) + 
				CString(_T(" TTL:")) + AllToString(mot.TimeToLength(0, di)), log_info);
			///
		}
	}

}
*/
