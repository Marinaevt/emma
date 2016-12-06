#include "StdAfx.h"
#include "2DRigid.h"


IOIMPL (C2DRigid, T2DRIGID)

C2DRigid::C2DRigid()
{
	RegisterMember(&m_shape);
	RegisterMember(&m_motion);
	RegisterMember(&m_stickpoint);
	RegisterMember(&m_FrictionCoeff);
	m_oldpos = Math::C2DPoint::Zero;
	m_IsStuck = false;
	m_IsFirstStuck = true;	//истина, т.к. будет первый контакт, далее по обстоятельствам
	ResetTime();
}

void C2DRigid::ResetTime()
{
	m_time = 0;
}

bool C2DRigid::Init()
{
	//LOGGER.Init(CString(_T("..\\..\\Logs\\C2DRigid.cpp_Init.txt")));
	Math::C2DPoint dpos;
	
	// Находим начальную точку траектории, это будет первоначальная "старая" позиция
	m_oldpos = m_motion.GetNode(0)->GetPoint();
	
	// Смещение Инструмента
	dpos = m_oldpos - m_stickpoint();

	// Меняем положение Инструмента
	for (size_t i = 0; i < m_shape.GetNodeCount(); i++){
		m_shape.GetNode(i)->SetPoint(m_shape.GetNode(i)->GetPoint() += dpos);
	}

	//m_shape.WriteToLog();
	return true;
}

void C2DRigid::Move(double dt)
{
	
	Math::C2DPoint pos, tau, dpos;
	
	DBL l = dynamic_cast<C2DFunction*>(m_motion.m_vels()[0])->GetIntegral(0, m_time + dt);
	m_motion.GetPoint(0, l, pos, tau);

	//Изменение позиции Инструмента
	dpos = pos - m_oldpos;

	for (size_t i = 0; i < m_shape.GetNodeCount(); i++){
		tau = m_shape.GetNode(i)->GetPoint();
		m_shape.GetNode(i)->SetPoint(tau + dpos);
	}

	m_oldpos = pos;
	m_time += dt;
	
}

void C2DRigid::Calc(double dt)
{
}

//! Точка внутри инструмента или нет (точность len)
bool C2DRigid::IsInside(const Math::C2DPoint& point, Math::C2DPoint& closep, DBL len)
{
	//точность пока не учитывается, ищется ближайший
	
	short inside = m_shape.IsInside(point);
	if (inside == 0 || inside == 1) {
		
		DBL dist;
		int nNode = m_shape.GetClosestNode(point,dist);	//получаем ближайший узел
		if(nNode == -1) return false;

		closep = m_shape.GetNode(nNode)->GetPoint();
		Math::C2DPoint minim;
		
		//находим все кривые с этим узлом
		for(size_t i=0; i<m_shape.GetCurveCount(); i++){

			C2DCurve *pCur = m_shape.GetCurve(i);
			if (pCur->GetStart() == nNode || pCur->GetEnd() == nNode){
				
				//находим ближайшую точку на кривой и сравниваем с предыдущей
				int p = pCur->GetClosestPoint(point,minim);
				if (p == -1) return false;
				if (dist > point.Len(minim)) {
					dist = point.Len(minim);
					closep = minim;
				}

				/*
				// если точность достаточна
				if(dist <= len){

				}*/
			}
		}
		
		return true;
	}
	return false;
}


bool C2DRigid::GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc)
{
	/*// не используем этот метод (раньше - да)
	if (!pMesh) return false;

	// Если с прошлого раза мы были присоединены, то это уже не первая стыковка
	// (значит трение или отцепились)
	if (m_IsStuck){
		m_IsFirstStuck = false;
	}

	m_IsStuck = false;	//заново проходим по граничным узлам (вдруг отцепились)

	//Устанавливаем ГУ только точкам внутри\на границе инуструмента
	if (m_shape.IsInside(pMesh->GetBorderNode(nBoundaryNode)) > 0) {
			
		m_IsStuck = true; //Пристыковались

		C2DPosition tmp;
		m_motion.GetPos(m_time, tmp);

		// <-- для осесимметричной задачи (х == 0, y == 0 - заделка)
		if (m_IsStuck && m_IsFirstStuck 
			|| fabs(m_FrictionCoeff - 1.0) < EPS
			|| (fabs(pMesh->GetBorderNode(nBoundaryNode).x) < EPS 
			//	&&	fabs(pMesh->GetBorderNode(nBoundaryNode).y) < EPS
			)){ 
			//Заделка
			bc.setKinematic(tmp.m_vel.m_x, tmp.m_vel.m_y);

		}else if (m_IsStuck && !m_IsFirstStuck){
	
			DBL dSigmaNormal_x = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::sigma_x);
			DBL dSigmaNormal_y = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::sigma_y);
			DBL dSigmaInt = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::int_s);
			DBL dSquare = pMesh->GetCircleSquare(nBoundaryNode);

			double dRes = Friction(dSigmaNormal_x,dSigmaNormal_y,dSigmaInt,dSquare);
			//dRes /= (1.000000001 - m_FrictionCoeff);

			//Случай для одноосного сжатия
			bc.setSymX(tmp.m_vel.m_y, dRes);
		}
	}

	//Если контакта не было (после перебора), то сбрасываем значения и отбой
	// отцепились
	if (!m_IsStuck){
		m_IsFirstStuck = true;	//следующая стыковка будет "первой"
		return false;
	}

	return true;
	//*/
	return false;	//тут остановки быть не должно
}


//используем этот метод
bool C2DRigid::GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom> *bc){
	
	if (!pMesh || pMesh == nullptr) {
		CDlgShowError cError(ID_ERROR_2DRIGID_MESHINTERFACE_NULL); //_T("Указатель на C2DMeshInterface null"));
		return false;
	}

	// Если с прошлого раза мы были присоединены, то это уже не первая стыковка
	// (значит трение или отцепились)
	if (m_IsStuck){
		m_IsFirstStuck = false;
	}

	m_IsStuck = false;	//заново проходим по граничным узлам (вдруг отцепились)

	size_t nSize = pMesh->m_bordernodes.GetSize();
	for(size_t nBoundaryNode=0; nBoundaryNode < nSize; nBoundaryNode++){
		//Устанавливаем ГУ только точкам внутри\на границе инуструмента
		if (m_shape.IsInside(pMesh->GetBorderNode(nBoundaryNode)) > 0) {
			
			m_IsStuck = true; //Пристыковались

			C2DPosition tmp;
			m_motion.GetPos(m_time, tmp);

			// <-- для осесимметричной задачи (х==0 - симметрия)
			if (m_IsStuck && m_IsFirstStuck 
				|| fabs(m_FrictionCoeff - 1.0) < EPS
				|| (fabs(pMesh->GetBorderNode(nBoundaryNode).x) < EPS 
				//	&& fabs(pMesh->GetBorderNode(nBoundaryNode).y) < EPS
				)){ 
				//Заделка
				bc->at(nBoundaryNode).setKinematic(tmp.m_vel.m_x, tmp.m_vel.m_y);

			}else if (m_IsStuck && !m_IsFirstStuck){
	
				DBL dSigmaNormal_x = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::sigma_x);
				DBL dSigmaNormal_y = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::sigma_y);
				DBL dSigmaInt = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::int_s);
				DBL dSquare = pMesh->GetCircleSquare(nBoundaryNode);

				double dRes = Friction(dSigmaNormal_x,dSigmaNormal_y,dSigmaInt,dSquare);
				//dRes /= (1.000000001 - m_FrictionCoeff);

				//Случай для одноосного сжатия
				bc->at(nBoundaryNode).setSymX(tmp.m_vel.m_y, dRes);
			}
		}
	}

	//Если контакта не было (после перебора), то сбрасываем значения и отбой
	// отцепились
	if (!m_IsStuck){
		m_IsFirstStuck = true;	//следующая стыковка будет "первой"
		return false;
	}

	return true;
}

Math::C2DRect C2DRigid::GetBoundingBox()
{
	Math::C2DRect rect;

	rect.AddRect(m_shape.GetBoundingBox());
	rect.AddRect(m_motion.GetBoundingBox());
	
	return rect;
}

bool C2DRigid::GetBoundingBox(CRect2D &rect){
	rect.AddRect(m_shape.GetBoundingBox());
	rect.AddRect(m_motion.GetBoundingBox());
	return true;
}

void C2DRigid::DrawGL(GLParam &parameter){ 
	
	glColor3ub(138, 51, 36);
	m_shape.DrawGL(parameter);
	//m_motion.DrawGL(parameter);
}


void C2DRigid::DrawGL3D(GLParam &parameter) {

	//m_shape.DrawGL(parameter);
	DBL j, z, z1;
	DBL  step;

	DBL angle = parameter.GetAngle(); //угол поворота
	int n = parameter.GetStep3D();   // число разбиений
	step = angle / DBL(n);            // "шаг" отрисвки
									  //DBL angle_1;
									  //angle_1= angle*0;

									  //glColor3ub(138, 51, 36);
									  //m_shape.DrawGL(parameter);//Отрисовываем первую заглушку
									  //m_motion.DrawGL(parameter);
	if (parameter.IsDrawInstrument()) {

		if (m_shape.GetContourCount())
		{
			//проверка на наличие нулевого контура , если есть , то заполняем кэш
			C2DContour *pContour = m_shape.GetContour(0);
			if (!pContour->IsCache())
			{

				pContour->FillCache();


			}
			//дополнительно сохраняем кэш
			size_t f = pContour->GetCache().size();

			glColor3ub(138, 51, 36);
			//m_shape.DrawGL(parameter);
			glBegin(GL_POLYGON);
			for (size_t i = 0; i < f; i++)
				glVertex2dv(pContour->GetCache()[i]);
			glEnd();

			//отрисовка второй заглушки
			glMatrixMode(GL_MODELVIEW);
			//glLoadIdentity();
			glPushMatrix();
			GLfloat gl_angle = float(-angle);
			GLfloat	gl_x = 0.0;
			GLfloat gl_y = 1.0;
			GLfloat gl_z = 0.0;
			glRotatef(gl_angle, gl_x, gl_y, gl_z);
			glBegin(GL_POLYGON);
			for (size_t i = 0; i < f; i++)
				glVertex2dv(pContour->GetCache()[i]);
			glEnd();

			glPopMatrix();

			//glLoadIdentity();

			//циклы для отрисовки инструмента - аналогично C2DPlaneFEM

			glClearDepth(1.0f);
			glEnable(GL_DEPTH_TEST);
			glDepthFunc(GL_LEQUAL);
			glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

			for (size_t i = 0; i < f; i++)
			{
				Math::C2DPoint node = pContour->GetCache()[i];
				Math::C2DPoint node1 = pContour->GetCache()[(i + 1) % f];

				for (j = 0;j<n;j++) {

					z = j*step * M_PI_180;
					z1 = ((j + 1)*step) * M_PI_180;

					glBegin(GL_QUADS);              //отрисовываем полигоны  квадратами
					glColor3d(0.7, 0.7, 0.7);
					glVertex3d((node.x)*cos(z), (node.y), (node.x)*sin(z));
					glVertex3d((node1.x)*cos(z), (node1.y), (node1.x)*sin(z));
					glVertex3d((node1.x)*cos(z1), (node1.y), (node1.x)*sin(z1));
					glVertex3d((node.x)*cos(z1), (node.y), (node.x)*sin(z1));

					glEnd();

					//glLineWidth(0.7f);
					//glColor3d(0,0,0);
					//glBegin(GL_LINES);  //отрисовываем черным линии уровня //стыки между полигонами
					//  glVertex3d( (node.x)*cos(z),  (node.y),(node.x)*sin(z));
					//  glVertex3d( (node.x)*cos(z1),  (node.y),(node.x)*sin(z1));
					//   glEnd();
				}
			}

			for (size_t i = 0; i < m_shape.GetContour(0)->GetCurveCount(); i++)
			{
				
				Math::C2DPoint node = m_shape.GetContour(0)->GetCurve(i)->GetStartNode()->GetPoint();

				for (j = 0; j<n; j++) {



					z = j*step * M_PI_180;
					z1 = ((j + 1)*step) * M_PI_180;



					glLineWidth(0.7f);
					glColor3d(0, 0, 0);
					glBegin(GL_LINES);  //отрисовываем черным линии уровня //стыки между полигонами
					glVertex3d((node.x)*cos(z), (node.y), (node.x)*sin(z));
					glVertex3d((node.x)*cos(z1), (node.y), (node.x)*sin(z1));
					glEnd();

				}
			}

		}

		

		//INT_PTR f = m_shape.GetContour(0)->GetCount();

		//циклы для отрисовки рёбер инструмента
		
	}


}
//------------------------------------------------------	
//bool C2DRigid::IsDrawInstrument() const
//{
//	return m_DrawInstrument;
//}
//
//void C2DRigid::SetDrawInstrument(bool draw_inst)
//{
//	m_DrawInstrument = draw_inst;
//}
//-------------------------------------------------------


void C2DRigid::Preparations(const ITask *task){

}

C2DRigid::~C2DRigid()
{
}

/**************************
	Трение / Friction
*************************/

/*	fric_type - тип закона трения
	1 - Кулона-Амонтона
	2 - Зибеля
	3 - Леванова-Колмогорова
//*/
DBL C2DRigid::ChooseFrictionMethod(short fric_type, DBL dFrictCoeff, DBL dSigma_s, DBL dNormalPressure){
	
	switch(fric_type){
		case 1: return FrictLawKulon(dFrictCoeff, dNormalPressure); break;
		case 2: return FrictLawZibel(dFrictCoeff, dSigma_s); break;
		case 3: return FrictLawLevanov(dFrictCoeff, dSigma_s, dNormalPressure); break;
		default: return 0.0;
	}
}

DBL C2DRigid::FrictLawKulon(DBL dFrictCoeff, DBL dNormalPressure){
	
	// F = mu*P
	// P - нормальное давление
	DBL Ftr = dFrictCoeff*dNormalPressure;
	return Ftr;
}

DBL C2DRigid::FrictLawZibel(DBL dFrictCoeff, DBL dSigma_s){
	
	// F = mu*sigma_s 
	// sigma_s - интенсивность напряжений
	DBL Ftr = dFrictCoeff*dSigma_s;
	return Ftr;
}

//всё по y
DBL C2DRigid::FrictLawLevanov(DBL dFrictCoeff, DBL dSigma_s, DBL dNormalPressure){
	
	// F = K*(sigma_s/sqrt(3))*(1 - EXP(-1.25*P/sigma_s))
	// K - фактор трения

	if (fabs(dSigma_s) < EPS) return 0.0;	//чтобы не \ на 0
	DBL Ftr = dFrictCoeff*(dSigma_s/sqrt(3))*(1.0 - exp(-1.25*dNormalPressure/dSigma_s));
	return Ftr;
}


DBL C2DRigid::Friction(DBL dSn_x, DBL dSn_y, DBL dS_int, DBL dSquare){
			
		int nType = 3; //тип трения	

		//Трение
		DBL dSigmaNormal_x = dSn_x,	//нормальные напряжения
			dSigmaNormal_y = dSn_y,
			dSigmaInt = dS_int;		//интенсивность напряжения

		//Пропускаем первый шаг
		if (fabs(m_time) < EPS){
			dSigmaNormal_x = 0.0;
			dSigmaNormal_y = 0.0;
			dSigmaInt = 0.0;
		}

		// Пока считаем
		DBL dFrictForce_x = ChooseFrictionMethod(nType,m_FrictionCoeff,dSigmaInt,dSigmaNormal_x),
			dFrictForce_y = ChooseFrictionMethod(nType,m_FrictionCoeff,dSigmaInt,dSigmaNormal_y);
				
		//DBL dRes = dFrictForce_y*dSquare;	//Если давим сверху, то сила трения по OY
		DBL dRes = dFrictForce_y;
				
		//Указываем направление
		/*
		if (dSigmaNormal_y < 0.0){ 
			dRes = abs(dRes);
		}else{
			dRes = (-1)*abs(dRes);
		}//*/
	return dRes;
}