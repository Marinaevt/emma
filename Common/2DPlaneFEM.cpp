#include "stdafx.h"
#include "2DPlaneFEM.h"


IOIMPL(C2DPlaneFEM, T2DPLANEFEM)

//! Временная реализация быстрой смены температуры и 
//! скорости деформации через init.txt
class CMuInfo {
public:
	CMuInfo() {
		m_dTemperature = 1253.15; //Кельвины, или 980 по Цельсию
		m_dee = 0.1; //с-1
		
		CString path = _T("init.txt");
		CFileException e;
		CStorage file;

		/*
		FILE *fStream;
		int e = _wfopen_s(&fStream, path, _T("rt,ccs=UNICODE"));
		if (e != 0) {
			CDlgShowError cError(ID_ERROR_2DPLANEFEM_INIT_NOT_LOADED);
		}

		CStorage file(fStream);  // open the file from this stream
		*/
		bool bRead1, bRead2;
		if (file.Open(path, CStorage::modeRead | CStorage::typeBinary, &e)) {
			bRead1 = ReadValue(file, "<Temperature>", &m_dTemperature);
			bRead2 = ReadValue(file, "<ee>", &m_dee);
			file.Close();
		}
		else {
			CDlgShowError cError(ID_ERROR_2DPLANEFEM_INIT_NOT_LOADED); //_T("init.txt not loaded (C2DPlaneFEM)"));
		}
	}
	~CMuInfo() {}

	DBL m_dTemperature;	//Температура
	DBL m_dee;			//Скорость деформации
};
CMuInfo g_MuInf;


///////////////////////////	C2DPlaneFEM

C2DPlaneFEM::C2DPlaneFEM()
{
	RegisterMember(&m_mesh_adapt);
	RegisterMember(&m_mat);
	RegisterMember(&m_bc);
	m_FirstVolSum = 0.0;
	m_CurrentVolSum = 0.0;
	m_FirstVolSumStart = false;
}

//! Пересохранение КЭ сетки при сборке в C2DBuilder
bool C2DPlaneFEM::InitMesh(const CMesh *pInitMesh) {

	C2DMesh3* mesh3 = new C2DMesh3(*pInitMesh);	//создание объекта C2DMesh3 на основе CMesh
	m_mesh_adapt.SetObject(mesh3);				//Сохранение созданного объекта сетки в задачу

	mesh3->InitNodeFields(0xffffffff);	//Инициализация полей в узлах
	mesh3->InitElemFields(0xffffffff);	//Инициализация полей в эл-ах

	mesh3->FillBorderNodes(); //Заполнение массива граничных эл-ов

	//m->GenerateMarks();	//Установление горизонтальных рисок

	return true;
}

//! Пересохранение Материала при сборке в C2DBuilder
bool C2DPlaneFEM::InitMaterial(const CString& path) {

	//Загрузка материала  //MLEG_corrected.spf, VT - 4.spf, VT6_old.spf, Amg6.spf
	CMaterial &mater = m_mat.Material();
	//int nLoad = mater.Load(CString("..\\..\\EMMAPrePost\\Data\\Properties\\Amg6.spf"));
	bool bLoad = m_mat.ReLoadMaterial(path);

	if (!mater.IsLoaded()) {
		CDlgShowError cError(ID_ERROR_2DPLANEFEM_MATERIAL_NOT_LOADED); //_T("Материал не загружен."));
		return false;
	}

	//Получаем уже загруженную сетку
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());
	if (mesh3 == nullptr) {
		CDlgShowError cError(ID_ERROR_2DPLANEFEM_FE_NOT_LOADED); //_T("КЭ сетка не загружена."));
		return false;
	}

	//Инициализация полей для расчётов
	//DBL mu = mater.Mu(0, g_MuInf.m_dee, g_MuInf.m_dTemperature);
	DBL mu = mater.Mu(0, g_MuInf.m_dee, g_MuInf.m_dTemperature, 1);// Albakri only (grain size=1)
	for (size_t i = 0; i < mesh3->m_elements().size(); i++) {
		mesh3->GetEField(i, eFields::mju) = mu;
		mesh3->GetEField(i, eFields::int_d) = 0.0;
		mesh3->GetEField(i, eFields::avg_d) = 0.0;
		mesh3->GetEField(i, eFields::smo_d) = 0.0;

		mesh3->GetEField(i, eFields::int_ss) = 1.0; // Albakri only (grain size)
	}
	return true;
}

void C2DPlaneFEM::InitBC() {

	//LOGGER.Init(CString(_T("..\\..\\Logs\\C2DPlaneFEM.cpp_InitBC.txt")));
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	/*for (int i = 0; i < m_bc.GetSize(); i++)
	{

	if (m->m_nodes()[m->m_bordernodes[i]].y <= 0.0001)	//низ
	{
	m_bc[i].setSymX();
	if (m->m_nodes()[m->m_bordernodes[i]].x == 0.0)	//левый нижний узел
	{
	m_bc[i].setKinematic(0,0);
	}
	}

	if (m->m_nodes()[m->m_bordernodes[i]].x <= 0.00001 && m->m_nodes()[m->m_bordernodes[i]].y != 0.0)
	{
	m_bc[i].setSymY();
	}


	}*/
}

//! Инициализация в Srv
bool C2DPlaneFEM::Init()
{
	//LOGGER.Init(CString(_T("..\\..\\Logs\\C2DPlaneFEM.cpp_Init.txt")));

	//Материал //MLEG_corrected.spf, VT6_old.spf, Amg6.spf
	//Теперь (12.02.2016) материал грузится из файла задачи (.emma_task)
#ifndef __LOCAL
	//m_mat.Material().Load(CString("..\\..\\EMMAPrePost\\Data\\Properties\\Amg6.spf"));
#else
	//m_mat.Material().Load(CString("./Properties.spf"));
#endif

	if (!m_mat.Material().IsLoaded()) {
		
		//TODO: загрузка временного файла, его загрузка в материал, удаление
		
		CDlgShowError cError(ID_ERROR_2DPLANEFEM_MATERIAL_NOT_LOADED);
		
		//TODO: загрузка временного файла, его загрузка в материал, удаление

		return false;
	}

	//Сетка и СЛАУ
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3 *>(m_mesh_adapt());
	m_mesh_adapt.SetObject(mesh3);

	m_slae.Init(mesh3->m_nodes().size() * 2, mesh3->m_nodes().size() * 2);

	//Инициализация ГУ
	m_bc.SetSize(mesh3->m_bordernodes().size()); //Задание массива граничных условий
	InitBC();

	//В лог
	//m->WriteToLog();
	//WriteBCToLog();

	return true;
}

//! Расчёт
double g_dError = -1;
int g_nK = -1;

class CErrInfo {
public:
	CErrInfo() {
		//количество перерасчётов
		m_N = 100;
		//точность вычислений
		m_dErr = 0.1;
		CStorage file;

		if (file.Open(_T("init.txt"), CStorage::modeRead | CStorage::typeBinary)) {
			double dVal;
			ReadValue(file, "<N_mju>", &dVal);
			m_N = static_cast<int>(dVal);
			ReadValue(file, "<MaxErr>", &m_dErr);
			file.Close();
		}
	}
	~CErrInfo() {}

	int m_N;//Количество перерасчётов;
	double m_dErr;//максимальное отклонение;
};
CErrInfo g_ErrInf;

void C2DPlaneFEM::Calc(DBL dt)
{

	m_CurTime() += dt;
	// проверка окончания расчета <- внешний процесс (+)
	// InitMju(); // инициализация mju - тупо заполнить значение mju в элементах числом 0.01 (+)
	// расчет линейной задачи CalcFEM в цикле -> нужна МЖ, правая часть и решение
	// расчет полей в контактных узлах
	// расчет полей в элементах
	// сглаживание средних деформаций
	// приращение времени <- внешний процесс (+)
	InitStep();
	CalcFEM(dt);

	CalcAvgDefSpeed();
	MakeSmoothDefSpeed();

	CalcAvgDef(dt);
	MakeSmoothDef();

	FieldsToNodes(dt);
	//CalcForce(dt);

	WriteEFieldsToLog();




	
	C2DMesh3* m = dynamic_cast<C2DMesh3 *>(m_mesh_adapt());

	double minY=0, maxY=0;
	ptrdiff_t iMin=0, iMax=0;
	bool bFirst = true;

	for (ptrdiff_t i = 0; i < m->GetBorderSize(); i++) {
		Math::C2DPoint p = m->GetBorderNode(i);
		if (p.x < 0.01) {
			if (bFirst || p.y < minY) {
				minY = p.y;
				iMin = i;
			}
			if (bFirst || p.y > maxY) {
				maxY = p.y;
				iMax = i;
			}
			bFirst = false;
		}
	}
	//double e_min = m->GetNField(m->m_bordernodes[iMin], int_d);
	double e_cp = m->GetNField(m->m_bordernodes[iMax], avg_d);
	double e_max = m->GetNField(m->m_bordernodes[iMax], int_d);
	//double ee_min = m->GetNField(m->m_bordernodes[iMin], int_ds);
	double dMju = m->GetNField(m->m_bordernodes[iMax], mju);
	double ee_max = m->GetNField(m->m_bordernodes[iMax], int_ds);
	//double s_min = m->GetNField(m->m_bordernodes[iMin], int_s);
	double d_max = m->GetNField(m->m_bordernodes[iMax], int_ss);
	double matMju = m_mat.Material().Mu(e_max, ee_max, 420, d_max);
	double s_max = m_mat.Material().Si(e_max, ee_max, 420, d_max);

	ALOGI("SA", AllToString(GetParent()->GetElapsedTime()) + _T(", ") +
	AllToString(dt) +   _T(", ") +
	AllToString(maxY) +   _T(", ") +
	AllToString(maxY-minY) + _T(", ")+
	//AllToString(e_min) + _T(", ") +
	AllToString(e_cp) + _T(", ") +
	AllToString(e_max) + _T(", ") +
	//AllToString(ee_min) + _T(", ") +
	AllToString(dMju) + _T(", ") +
	AllToString(ee_max) + _T(", ") +
	AllToString(matMju) + _T(", ") +
	//AllToString(s_min) + _T(", ") +
	AllToString(s_max) + _T(", ") +
	AllToString(g_dError) + _T(", ") +
	AllToString(g_nK) + _T(", ") +
	AllToString(d_max) + _T(", ")
	);
	
}

void C2DPlaneFEM::InitStep()
{
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());
	m_slae.ZeroAll();

	//LOGGER.Init(CString("..\\..\\Logs\\C2DPlaneFEM.cpp_InitStep.txt"));
	//m->WriteToLog();
	//WriteBCToLog();

	//площадь и средние точки (sqr, cm_x, cm_y)
	mesh3->CalcFormCoeffs();

	m_CurrentVolSum = 0.0;

	size_t nElemSize = mesh3->m_elements().size();
	for (size_t i = 0; i < nElemSize; i++) {
		//m->GetEField(i, eFields::mju) = m_mat.Material().Mu(m->GetEField(i, eFields::int_d), m->GetEField(i, eFields::int_ds), TEMPERATURE);// Расходится((
		//mesh3->GetEField(i, eFields::mju) = m_mat.Material().Mu(mesh3->GetEField(i, eFields::int_d), mesh3->GetEField(i, eFields::int_ds), g_MuInf.m_dTemperature);
		mesh3->GetEField(i, eFields::mju) = m_mat.Material().Mu(mesh3->GetEField(i, eFields::int_d), mesh3->GetEField(i, eFields::int_ds), g_MuInf.m_dTemperature, mesh3->GetEField(i, eFields::int_ss)); // Albakri only
		//m->GetEField(i, eFields::smo_d) = 0;
		if (!m_FirstVolSumStart) {
			m_FirstVolSum += mesh3->GetEField(i, eFields::sqr) * 2 * M_PI*mesh3->GetEField(i, eFields::cm_x);
		}
		m_CurrentVolSum += mesh3->GetEField(i, eFields::sqr) * 2 * M_PI*mesh3->GetEField(i, eFields::cm_x);
	}

	m_FirstVolSumStart = true;	//установили, чтобы не обновлять

								/*for (int i = 0; i < m->m_nodes().GetCount(); i++){

								m->GetNField(i, eFields::smo_d) = 0.0;	//Среднее давление в узле
								m->GetNField(i, eFields::sqr) = 0.0;	//Площадь в узле
								} //*/

}

void C2DPlaneFEM::CalcFEM(DBL dt)
{
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	//ALOGGER.AddLog(CString("..\\..\\Logs\\Log_mSLAE_matr.txt"), CString(_T("YD1")), false);
	//ALOGGER.AddLog(CString("..\\..\\Logs\\Log_mSLAE_rp.txt"), CString(_T("YD2")), false);
	LOGGER.Init(CString("..\\..\\Logs\\C2DPlaneFEM.cpp_CalcFEM.txt"));
	LOGGER.SetTypeStamp(false);
	LOGGER.SetTimeStampFilter(0);
	//m->WriteToLog();
	//WriteBCToLog();//*/

	std::vector<DBL> prevSol;

	//m->CalcFormCoeffs();	//уже есть в Init()
	//WriteEFieldsToLog();

	for (int k = 0; k < g_ErrInf.m_N; k++)		//Итерации по нелинейности
	{

		clock_t start, end;
		

		m_slae.ZeroAll();
		//m_slae.WriteToLog();		

		//LOGGER.Init(CString("..\\..\\Logs\\C2DPlaneFEM.cpp_CalcFEM_FILLSTM.txt"));
		mesh3->FillSTM(m_slae.m_matr, m_slae.m_rp, m_mat.Material(), dt, 0, mesh3->m_elements().size());
		//m_slae.WriteToLog();

		// внесение ГУ
		size_t bordernodes_size = mesh3->m_bordernodes().size();
		for (size_t i = 0; i < bordernodes_size; i++)
		{
			m_slae.SetBC(mesh3->m_bordernodes[i] * 2, m_bc[i]); // не уверен насчет умножения на 2
		}
		//задали правую часть, записываем в лог
		//m_slae.WriteToLogRightPart();
		//m_slae.WriteToLogMatrix();
		
		//m_slae.m_matr.Symmetrize();

		// решение
		start = clock();
		
		//LOGGER.Init(CString("..\\..\\Logs\\C2DPlaneFEM.cpp_CalcFEM_Gauss.txt"));
		m_slae.Gauss();
		
		//m_slae.WriteToLog();
		//m_slae.WriteToLogSolution();
		//WriteToLogNewLines(2);
		
		end = clock();
		// вычисление IntDef
		CalcIntDef(dt);

		//вычисление Mju
		size_t elem_size = mesh3->m_elements().size();
		for (size_t i = 0; i < elem_size; i++) {
			//Было m->GetEField(i, eFields::int_ds) -> g_MuInf.dee
			mesh3->GetEField(i, eFields::mju) = m_mat.Material().Mu(mesh3->GetEField(i, eFields::int_d), mesh3->GetEField(i, eFields::int_ds), g_MuInf.m_dTemperature, mesh3->GetEField(i, eFields::int_ss)); // Albakri only
																																																			  //mesh3->GetEField(i, eFields::mju) = m_mat.Material().Mu(mesh3->GetEField(i, eFields::int_d), mesh3->GetEField(i, eFields::int_ds), g_MuInf.m_dTemperature);
		}

		g_nK = k;
		if (k > 3) {
			g_dError = MaxErrorP(prevSol, m_slae.m_sol);
			if (g_dError < g_ErrInf.m_dErr) break;
		}
		prevSol = m_slae.m_sol; // копируем предыдущее решение

		
		//OutputDebugString(_T("The above code block was executed in %.4f second(s)\n") + AllToString(((double)end - start) / ((double)CLOCKS_PER_SEC)));
	}

	//if(g_nK<7)	m_parent->RecommendDt(dt*1.1);
	//if(g_nK>15)m_parent->RecommendDt(dt*0.9);

	// заполнение МЖ (+)
	// расчет трения (везде устанавливаем режим slide)
	// выталкивание на штамп
	// обнуление столбцов у известных компонент (внесение ГУ)
	// решение линейной задачи методом Гаусса (решение получаем в правой части)
	// восстановление решения
	// перенос решения из правой части в вектор решения
	// вычисление полей интенсивности скорости деформации
	// DBL CalcMju(); // вычисление mju
}

void C2DPlaneFEM::WriteToLogNewLines(int numNewLines) {
	for (int i = 0; i < numNewLines; i++) {
		DLOG(CString(_T("")), log_info);
	}
}


//! для CalcFEM
void C2DPlaneFEM::CalcIntDef(DBL dt)
{
	/*
	// ijk - индексы узлов текущего элемента
	// b[], c[] - коэффициенты функции формы
	// m_C[] - решение
	// sm - площадь элемента

	double exx = (b[0]*m_C[2*ijk[0]]+b[1]*m_C[2*ijk[1]]+b[2]*m_C[2*ijk[2]])/2/sm;
	double eyy = (c[0]*m_C[2*ijk[0]+1]+c[1]*m_C[2*ijk[1]+1]+c[2]*m_C[2*ijk[2]+1])/2/sm;
	double exy = (c[0]*m_C[2*ijk[0]]+c[1]*m_C[2*ijk[1]]+c[2]*m_C[2*ijk[2]]+b[0]*m_C[2*ijk[0]+1]+b[1]*m_C[2*ijk[1]+1]+b[2]*m_C[2*ijk[2]+1])/4/sm;
	*/
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	DBL B[3], C[3];
	DBL exx, eyy, exy;

	size_t elem_size = mesh3->m_elements().size();
	for (size_t ei = 0; ei < elem_size; ei++)
	{
		C2DElement3& e = mesh3->m_elements()[ei];

		Math::C2DPoint& n0 = mesh3->m_nodes()[e._0]; // 0 - i
		Math::C2DPoint& n1 = mesh3->m_nodes()[e._1]; // 1 - j
		Math::C2DPoint& n2 = mesh3->m_nodes()[e._2]; // 2 - m k

		B[0] = n1.m_y - n2.m_y; //Y[j]-Y[k];
		B[1] = n2.m_y - n0.m_y; //Y[k]-Y[i];
		B[2] = n0.m_y - n1.m_y; //Y[i]-Y[j];
		C[0] = n2.m_x - n1.m_x; //X[k]-X[j];
		C[1] = n0.m_x - n2.m_x; //X[i]-X[k];
		C[2] = n1.m_x - n0.m_x; //X[j]-X[i];

		DBL sm = abs(n0.m_x * B[0] + n2.m_x * B[2] + n1.m_x * B[1]) * 0.5; //(X[i]*(B[0]) + X[k]*(B[2]) + X[j]*(B[1]))/2;

		exx = (B[0] * m_slae.m_sol[2 * e.n[0]] + B[1] * m_slae.m_sol[2 * e.n[1]] + B[2] * m_slae.m_sol[2 * e.n[2]]) / 2 / sm;
		eyy = (C[0] * m_slae.m_sol[2 * e.n[0] + 1] + C[1] * m_slae.m_sol[2 * e.n[1] + 1] + C[2] * m_slae.m_sol[2 * e.n[2] + 1]) / 2 / sm;

		exy = (C[0] * m_slae.m_sol[2 * e.n[0]] + C[1] * m_slae.m_sol[2 * e.n[1]] + C[2] * m_slae.m_sol[2 * e.n[2]] +
			B[0] * m_slae.m_sol[2 * e.n[0] + 1] + B[1] * m_slae.m_sol[2 * e.n[1] + 1] + B[2] * m_slae.m_sol[2 * e.n[2] + 1]) / 4 / sm;


		mesh3->GetEField(ei, eFields::int_ds) = sqrt((exx - eyy) * (exx - eyy) + eyy * eyy + exx * exx + 6 * (exy * exy)) * sqrt(2.0) / 3.0;

		//Для осесимметричных задач{{
		DBL e_fi = (m_slae.m_sol[2 * e.n[0]] + m_slae.m_sol[2 * e.n[1]] + m_slae.m_sol[2 * e.n[2]]) / (n0.m_x + n1.m_x + n2.m_x);
		mesh3->GetEField(ei, eFields::int_ds) = sqrt((exx - eyy) * (exx - eyy) + (eyy - e_fi) * (eyy - e_fi) + (exx - e_fi) * (exx - e_fi) + 6 * (exy * exy)) * sqrt(2.0) / 3.0;
		//сохраняем компоненты "х" и "у"
		mesh3->GetEField(ei, eFields::sf_a) = exx;
		mesh3->GetEField(ei, eFields::sf_b) = eyy;
		mesh3->GetEField(ei, eFields::sf_c) = e_fi;
		//Для осесимметричных задач}}*/

		//НЕ прибавлять
		//m->GetEField(ei, eFields::avg_d) += (exx + eyy) * dt * 0.5;

	}

}

void C2DPlaneFEM::CalcAvgDefSpeed() {

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());
	DBL B[3], C[3];
	DBL exx, eyy, e_fi;

	size_t elem_size = mesh3->m_elements().size();
	for (size_t ei = 0; ei < elem_size; ei++)
	{
		C2DElement3& e = mesh3->m_elements()[ei];

		Math::C2DPoint& n0 = mesh3->m_nodes()[e._0]; // 0 - i
		Math::C2DPoint& n1 = mesh3->m_nodes()[e._1]; // 1 - j
		Math::C2DPoint& n2 = mesh3->m_nodes()[e._2]; // 2 - m k

		B[0] = n1.m_y - n2.m_y; //Y[j]-Y[k];
		B[1] = n2.m_y - n0.m_y; //Y[k]-Y[i];
		B[2] = n0.m_y - n1.m_y; //Y[i]-Y[j];
		C[0] = n2.m_x - n1.m_x; //X[k]-X[j];
		C[1] = n0.m_x - n2.m_x; //X[i]-X[k];
		C[2] = n1.m_x - n0.m_x; //X[j]-X[i];

								//Площадь эл-та
		DBL sm = abs(n0.m_x * B[0] + n2.m_x * B[2] + n1.m_x * B[1]) * 0.5; //(X[i]*(B[0]) + X[k]*(B[2]) + X[j]*(B[1]))/2;
		if (fabs(sm) < EPS) {
			CDlgShowError cError(ID_ERROR_2DPLANEFEM_SQ_TRIANGLE_LESS_EPS); //_T("Площадь треугольника меньше EPS \n"));
			return;
		}

		exx = (B[0] * m_slae.m_sol[2 * e.n[0]] + B[1] * m_slae.m_sol[2 * e.n[1]] + B[2] * m_slae.m_sol[2 * e.n[2]]) / 2 / sm;
		eyy = (C[0] * m_slae.m_sol[2 * e.n[0] + 1] + C[1] * m_slae.m_sol[2 * e.n[1] + 1] + C[2] * m_slae.m_sol[2 * e.n[2] + 1]) / 2 / sm;

		/*
		// Отдельные компоненты для отладки
		DBL xx_0 = m_slae.m_sol[2 * e.n[0]],
		xx_1 = m_slae.m_sol[2 * e.n[1]],
		xx_2 = m_slae.m_sol[2 * e.n[2]];
		DBL yy_0 = m_slae.m_sol[2 * e.n[0] + 1],
		yy_1 = m_slae.m_sol[2 * e.n[1] + 1],
		yy_2 = m_slae.m_sol[2 * e.n[2] + 1];
		//*/

		/*//Для плоских задач{{
		m->GetEField(ei, eFields::avg_ds) += (exx + eyy) / 2.0;
		//Для плоских задач}}
		//
		/*///Для осесимметричных задач{{
		e_fi = (m_slae.m_sol[2 * e.n[0]] + m_slae.m_sol[2 * e.n[1]] + m_slae.m_sol[2 * e.n[2]]) / (n0.m_x + n1.m_x + n2.m_x);

		mesh3->GetEField(ei, eFields::avg_ds) = (exx + eyy + e_fi) / 3.0;
		//Для осесимметричных задач}}*/
	}

}

void C2DPlaneFEM::MakeSmoothDefSpeed() {

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	size_t nodes_size = mesh3->m_nodes().size();
	for (size_t i = 0; i < nodes_size; i++) {
		mesh3->GetNField(i, eFields::smo_ds) = 0.0;
		mesh3->GetNField(i, eFields::sqr) = 0.0;
	}

	/*//умножаем на х
	for (int i = 0; i < m->m_elements().GetCount(); i++){
	m->GetEField(i, eFields::avg_ds) *= m->GetEField(i, eFields::cm_x);
	}//*/

	size_t elem_size = mesh3->m_elements().size();
	for (size_t ei = 0; ei < elem_size; ei++)
	{
		C2DElement3& e = mesh3->m_elements()[ei];

		//Для плоских задач{{
		/*DBL x =1;
		//Для плоских задач}}
		/*///Для осесимметричных задач{{ Не уверен тут..
		DBL x = (e.n[0] + e.n[1] + e.n[2]) / 3.0;
		//Для осесимметричных задач}}*/
		for (int j = 0; j < 3; j++) {
			mesh3->GetNField(e.n[j], eFields::smo_ds) += mesh3->GetEField(ei, eFields::avg_ds) * mesh3->GetEField(ei, eFields::sqr) * mesh3->GetEField(ei, eFields::cm_x);
			//m->GetNField(e.n[j], eFields::smo_ds) += m->GetEField(ei, eFields::avg_ds) * m->GetEField(ei, eFields::sqr);
			//NField Volume
			mesh3->GetNField(e.n[j], eFields::sqr) += mesh3->GetEField(ei, eFields::sqr) * mesh3->GetEField(ei, eFields::cm_x);
			//m->GetNField(e.n[j], eFields::sqr) += m->GetEField(ei, eFields::sqr);
		}
	}

	// делим полученное поле на (площади * x) == объёмы
	for (size_t i = 0; i < nodes_size; i++) {
		mesh3->GetNField(i, eFields::smo_ds) /= mesh3->GetNField(i, eFields::sqr);
	}

	// усредняем итоговое поле из узлов на элементы
	for (size_t i = 0; i < elem_size; i++)
	{
		mesh3->GetEField(i, eFields::smo_ds) = 0;
		C2DElement3& es = mesh3->m_elements()[i];
		for (int j = 0; j < 3; j++) {
			mesh3->GetEField(i, eFields::smo_ds) += mesh3->GetNField(es.n[j], eFields::smo_ds);
		}
		mesh3->GetEField(i, eFields::smo_ds) = mesh3->GetEField(i, eFields::smo_ds) / 3.0;
	}

	/*//умножаем на х
	for (int i = 0; i < m->m_elements().GetCount(); i++){
	m->GetEField(i, eFields::avg_ds) /= m->GetEField(i, eFields::cm_x);
	//m->GetEField(i, eFields::smo_ds) /= m->GetEField(i, eFields::cm_x);
	}//*/
}

void C2DPlaneFEM::CalcAvgDef(DBL dt)
{
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	size_t elem_size = mesh3->m_elements().size();
	for (size_t ei = 0; ei < elem_size; ei++)
	{
		mesh3->GetEField(ei, eFields::avg_d) += mesh3->GetEField(ei, eFields::smo_ds) * dt;
		//m->GetEField(ei, eFields::avg_d) = (m_CurrentVolSum-m_FirstVolSum)/m_FirstVolSum;
	}
}

void C2DPlaneFEM::MakeSmoothDef()
{
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	//LOGGER.Init(CString("..\\..\\Logs\\C2DPlaneFEM.cpp_MakeSmooth.txt"));
	//m->WriteToLog();
	//WriteBCToLog();
	//WriteEFieldsToLog();

	size_t nodes_size = mesh3->m_nodes().size();
	for (size_t i = 0; i < nodes_size; i++) {
		mesh3->GetNField(i, eFields::smo_d) = 0.0;
		mesh3->GetNField(i, eFields::sqr) = 0.0;
	}

	/*//умножаем на х
	for (int i = 0; i < m->m_elements().GetCount(); i++){
	m->GetEField(i, eFields::avg_d) *= m->GetEField(i, eFields::cm_x);
	}//*/

	// распихиваем по узлам среднюю деформацию в элементе * площадь элемента, а также суммируем площади прилегающих к узлу элементов
	size_t elem_size = mesh3->m_elements().size();
	for (size_t i = 0; i < elem_size; i++)
	{
		C2DElement3& e = mesh3->m_elements()[i];

		//Для плоских задач{{
		/*DBL x =1;
		//Для плоских задач}}
		/*/
		//Для осесимметричных задач{{ Не уверен тут..
		DBL x = (e.n[0] + e.n[1] + e.n[2]) / 3.0;
		//Для осесимметричных задач}}*/
		for (int j = 0; j < 3; j++)
		{
			//was *x (at the end)
			mesh3->GetNField(e.n[j], eFields::smo_d) += mesh3->GetEField(i, eFields::avg_d) * mesh3->GetEField(i, eFields::sqr) * mesh3->GetEField(i, eFields::cm_x);
			//m->GetNField(e.n[j], eFields::smo_d) += m->GetEField(i, eFields::avg_d) * m->GetEField(i, eFields::sqr);
			//nField Volume !
			mesh3->GetNField(e.n[j], eFields::sqr) += mesh3->GetEField(i, eFields::sqr) * mesh3->GetEField(i, eFields::cm_x);
			//m->GetNField(e.n[j], eFields::sqr) += m->GetEField(i, eFields::sqr);
		}
	}

	// делим полученное поле на суммарные площади 
	// и на 3 (т.к. суммировали)
	for (size_t i = 0; i < nodes_size; i++) {
		mesh3->GetNField(i, eFields::smo_d) /= (mesh3->GetNField(i, eFields::sqr));		//включает уже *х
	}

	// усредняем итоговое поле из узлов на элементы
	for (size_t i = 0; i < elem_size; i++)
	{
		mesh3->GetEField(i, eFields::smo_d) = 0;
		C2DElement3& es = mesh3->m_elements()[i];
		for (int j = 0; j < 3; j++) {
			mesh3->GetEField(i, eFields::smo_d) += mesh3->GetNField(es.n[j], eFields::smo_d);
		}
		mesh3->GetEField(i, eFields::smo_d) = mesh3->GetEField(i, eFields::smo_d) / 3.0;
	}

	/*//умножаем на х
	for (int i = 0; i < m->m_elements().GetCount(); i++){
	m->GetEField(i, eFields::avg_d) /= m->GetEField(i, eFields::cm_x);
	//m->GetEField(i, eFields::smo_d) /= m->GetEField(i, eFields::cm_x);
	}//*/

}

void C2DPlaneFEM::FieldsToNodes(DBL dt) {

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());
	DBL d_int_d, d_int_ds, d_avg_d;	//поля
	DBL d_sigma_x, d_sigma_y, d_sigma_xy;	//поля напряжений по х, у, ху.
	DBL d_V_el, d_V_node;	//объём одного эл-та и всех эл-ов вокруг узла
	DBL d_mju, d_avg_ds;
	size_t nEl_count = 0;	//кол-во эл-ов, содержащих нашу точку

	//! Проходим по всем узлам, находим содержащие их элементы
	size_t nodes_size = mesh3->m_nodes().size();
	for (size_t i = 0; i < nodes_size; i++)
	{
		//! Сброс всех сохранённых ранее значений
		d_V_el = 0.0;
		d_V_node = 0.0;
		d_int_d = d_int_ds = 0.0;
		d_avg_d = d_avg_ds = 0.0;
		d_sigma_x = d_sigma_y = d_sigma_xy = 0.0;
		nEl_count = 0;

		DBL d_ss = 0; // Albakri only (SA)

		DBL dIntS = 0.0;
		DBL dMju = 0.0;

		//По всем элементам
		size_t elem_size = mesh3->m_elements().size();
		for (size_t j = 0; j<elem_size; j++) {
			if (mesh3->m_elements()[j].IsHaveNode(i)) {

				// V = S*(Pi*r*r)
				d_V_el = mesh3->GetEField(j, eFields::sqr) * (mesh3->GetEField(j, eFields::cm_x));
				d_int_d += mesh3->GetEField(j, eFields::int_d) * d_V_el;
				d_int_ds += mesh3->GetEField(j, eFields::int_ds) * d_V_el;

				d_avg_d += mesh3->GetEField(j, eFields::avg_d) * d_V_el;
				d_avg_ds += mesh3->GetEField(j, eFields::avg_ds) * d_V_el;
				
				d_ss += mesh3->GetEField(j, eFields::int_ss) * d_V_el; // Albakri only (SA)

				d_mju = mesh3->GetEField(j, eFields::mju);

				//sigma_x = 2*mu*e_dot_x + (K*dt - 2/3*mu)*tetta_dot + sigma(*);
				// CalcK(dt, mu) = (K*dt - 2/3*mu);
				d_sigma_x += (2.0*d_mju*mesh3->GetEField(j, eFields::sf_a) + m_mat.Material().CalcK(dt, d_mju)*(3.0*mesh3->GetEField(j, eFields::avg_ds)) + mesh3->GetEField(j, eFields::avg_d)) * d_V_el;
				d_sigma_y += (2.0*d_mju*mesh3->GetEField(j, eFields::sf_b) + m_mat.Material().CalcK(dt, d_mju)*(3.0*mesh3->GetEField(j, eFields::avg_ds)) + mesh3->GetEField(j, eFields::avg_d)) * d_V_el;
				//d_int_sxy += (2.0*d_mju*m->GetEField(j, eFields::int_ds)  * d_V_el;

				//sigma_int = 3*mu*e_dot_int;
				dIntS += 3.0*d_mju* mesh3->GetEField(j, eFields::int_ds) * d_V_el;
				dMju += d_mju*d_V_el;

				nEl_count++;
				d_V_node += d_V_el;
			}
		}

		//d_V_node не может быть == 0, иначе есть висячие узлы - расчёт невозможен
		mesh3->GetNField(i, eFields::int_d) = d_int_d / d_V_node;
		mesh3->GetNField(i, eFields::int_ds) = d_int_ds / d_V_node;

		mesh3->GetNField(i, eFields::avg_d) = d_avg_d / d_V_node;
		mesh3->GetNField(i, eFields::avg_ds) = d_avg_ds / d_V_node;

		//у напряжения пока интересует только "y" компонента
		mesh3->GetNField(i, eFields::sigma_x) = d_sigma_x / d_V_node;
		mesh3->GetNField(i, eFields::sigma_y) = d_sigma_y / d_V_node;

		mesh3->GetNField(i, eFields::int_s) = dIntS / d_V_node;
		mesh3->GetNField(i, eFields::mju) = dMju / d_V_node;

		mesh3->GetNField(i, eFields::int_ss)  = d_ss/ d_V_node; // Albakri only (SA)
	}

	/*
	CString sTmp = CString(_T(""));
	LOG(CString(_T("Writing EFields:")), log_info);
	LOG(AllToString(dt), log_info);
	LOG(CString(_T("No | mju | eyy_s | int_d | int_ds | avg_d | smo_d | exx | eyy | efi | CalcK(mju,dt) | s_y")), log_info);

	//LOG
	sTmp = AllToString(i)+CString(_T("-"))+AllToString(j)+CString(_T(": "));

	sTmp += AllToString(m->GetEField(j,mju))+CString(_T(" | "));
	sTmp += AllToString(m->GetEField(j,nju))+CString(_T(" | "));

	//sTmp += AllToString(m_mat.Material().GetK())+CString(_T(" | "));
	sTmp += AllToString(m->GetEField(j,int_d))+CString(_T(" | "));
	sTmp += AllToString(m->GetEField(j,int_ds))+CString(_T(" | "));
	sTmp += AllToString(m->GetEField(j,avg_d))+CString(_T(" | "));
	sTmp += AllToString(m->GetEField(j,smo_d))+CString(_T(" | "));

	sTmp += AllToString(m->GetEField(j,sf_a))+CString(_T(" | "));
	sTmp += AllToString(m->GetEField(j,sf_b))+CString(_T(" | "));
	sTmp += AllToString(m->GetEField(j,sf_c))+CString(_T(" | "));
	sTmp += AllToString(m_mat.Material().CalcK(dt,d_mju))+CString(_T(" | "));


	sTmp += AllToString(d_int_sy)+CString(_T(" | "));

	LOG(sTmp,log_info); //*/
}

void C2DPlaneFEM::CalcForce(DBL dt) {

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	//Считаем силу, сигмы посчитаны в FieldsToNodes()

	DBL d_FirstSigma = 0.0, //для одноосного сжатия sigma_y
		d_SecondSigma = 0.0,
		d_Sigma_Sqr = 0.0, d_V_Sigma = 0.0,
		d_Force = 0.0,
		d_Sigma = 0.0;
	size_t nCount = 0; //сколько нижних эл-ов

	DBL dVolumeTest = 0.0,
		dSqrTest = 0.0;

	size_t elem_size = mesh3->m_elements().size();
	for (size_t j = 0; j<elem_size; j++) {

		Math::C2DPoint p0 = mesh3->m_nodes()[mesh3->m_elements()[j].n0],
			p1 = mesh3->m_nodes()[mesh3->m_elements()[j].n1],
			p2 = mesh3->m_nodes()[mesh3->m_elements()[j].n2];

		d_Sigma_Sqr = 0.0;	// зануляем sqr, чтобы не суммировать ненужные эл-ты
							// находим эл-ты с узлами с 0 по оси "OY"

							// находим те из них, где есть 2 таких узла
							// вычисляем dF, потом суммируем
		if (fabs(p0.y) < EPS) {
			if (fabs(p1.y) < EPS) {
				//считаем площадь
				d_FirstSigma = mesh3->GetNField(mesh3->m_elements()[j].n0, eFields::sigma_y);	//было int_s
				d_SecondSigma = mesh3->GetNField(mesh3->m_elements()[j].n1, eFields::sigma_y);
				d_Sigma_Sqr = 0.5*(d_FirstSigma + d_SecondSigma)*abs(p0.x - p1.x);
				nCount++;

			}
			else if (fabs(p2.y) < EPS) {
				d_FirstSigma = mesh3->GetNField(mesh3->m_elements()[j].n0, eFields::sigma_y);	//было int_s
				d_SecondSigma = mesh3->GetNField(mesh3->m_elements()[j].n2, eFields::sigma_y);
				d_Sigma_Sqr = 0.5*(d_FirstSigma + d_SecondSigma)*abs(p0.x - p2.x);
				nCount++;
			}
		}
		else if (fabs(p1.y) < EPS) {
			if (fabs(p0.y) < EPS) {
				d_FirstSigma = mesh3->GetNField(mesh3->m_elements()[j].n0, eFields::sigma_y);	//было int_s
				d_SecondSigma = mesh3->GetNField(mesh3->m_elements()[j].n1, eFields::sigma_y);
				d_Sigma_Sqr = 0.5*(d_FirstSigma + d_SecondSigma)*abs(p0.x - p1.x);
				nCount++;

			}
			else if (fabs(p2.y) < EPS) {
				d_FirstSigma = mesh3->GetNField(mesh3->m_elements()[j].n1, eFields::sigma_y);	//было int_s
				d_SecondSigma = mesh3->GetNField(mesh3->m_elements()[j].n2, eFields::sigma_y);
				d_Sigma_Sqr = 0.5*(d_FirstSigma + d_SecondSigma)*abs(p1.x - p2.x);
				nCount++;

			}
		}
		//p2 проверять не надо, т.к. иначе p0 || p1 != 0 по "y"

		d_V_Sigma = d_Sigma_Sqr*(2 * M_PI*mesh3->GetEField(j, eFields::cm_x));	// V = Sigma*(2*pi*r)
		d_Force += d_V_Sigma;
		d_Sigma += 0.5*(d_FirstSigma + d_SecondSigma);

		//log
		dVolumeTest += 2 * M_PI*mesh3->GetEField(j, eFields::cm_x) * mesh3->GetEField(j, eFields::sqr);
		dSqrTest += mesh3->GetEField(j, eFields::sqr);
		//dAVGD += m->GetEField(j, eFields::avg_d);
	}

	DBL d_avgSigma = d_Sigma / nCount;	// Средние интенсивности напряжений в узлах

										//Расчёт давления (F/S), вычисляем площадь круга основания заготовки
	DBL dMax_x = 0.0;
	DBL dMax_y = 0.0;
	DBL dCircleSqr = 0.0;

	DBL dMu = 0.0,
		dCalcK = 0.0;

	//Граничные точки заготовки
	size_t bordernodes_size = mesh3->m_bordernodes().size();
	for (size_t j = 0; j < bordernodes_size; j++) {
		Math::C2DPoint pBorderPoint = mesh3->GetBorderNode(j);
		if (pBorderPoint.m_y < EPS) {	// т.е. m_y == 0

			if (pBorderPoint.m_x > dMax_x) {
				dMax_x = pBorderPoint.m_x;	//Запоминаем крайнюю правую точку у основания заготовки
			}

			//if (m->GetBorderNode(j).m_y == 0){
			dMu = mesh3->GetNField(j, eFields::mju);
			dCalcK = m_mat.Material().CalcK(dt, dMu);
			//}	//*/
		}

		if (pBorderPoint.m_y > dMax_y) {
			dMax_y = pBorderPoint.m_y;	//Запоминаем крайнюю верхнюю точку заготовки
		}
	}

	dCircleSqr = M_PI*dMax_x*dMax_x;	// S = pi*r*r
										//m->GetNField(0, eFields::sf_b) = dCircleSqr; 
										//m->GetNField(0, eFields::sf_c) = d_Force / dCircleSqr; //TODO: Сохраняем давление в sf_c временно
	DBL dPressure = d_Force / dCircleSqr;

	CString tmp = AllToString(m_CurTime()) + CString(_T(" | "));
	tmp += AllToString(dt) + CString(_T(" | "));
	tmp += AllToString(d_Force) + CString(_T(" | "));
	tmp += AllToString(d_avgSigma) + CString(_T(" | "));

	tmp += AllToString(dCircleSqr) + CString(_T(" | "));
	tmp += AllToString(dPressure) + CString(_T(" | "));

	tmp += AllToString(dMax_y) + CString(_T(" | "));
	tmp += AllToString(dMax_x) + CString(_T(" | "));

	tmp += AllToString(g_dError) + CString(_T(" | "));
	tmp += AllToString(g_nK) + CString(_T(" | "));

	//tmp += AllToString(dMu) + CString(_T(" | "));
	//tmp += AllToString(dCalcK) + CString(_T(" | "));

	tmp += AllToString(dVolumeTest) + CString(_T(" | "));
	tmp += AllToString(dSqrTest) + CString(_T(" | "));
	tmp += AllToString(m_CurrentVolSum) + CString(_T(" | "));
	tmp += AllToString(m_FirstVolSum) + CString(_T(" | "));

	//tmp += AllToString(m_mat.Material().GetK(dt, dMu)) + CString(_T(" | "));
	ALOGI("YD", tmp);
}


void C2DPlaneFEM::Move(DBL dt)
{
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	//LOGGER.Init(CString("..\\..\\Logs\\C2DPlaneFEM.cpp_Move.txt"));
	//m->WriteToLog();
	//WriteBCToLog();

	for (size_t i = 0; i < mesh3->m_nodes().size(); i++)
	{
		mesh3->m_nodes()[i].m_x += m_slae.m_sol[2 * i] * dt;
		mesh3->m_nodes()[i].m_y += m_slae.m_sol[2 * i + 1] * dt;
	}

	for (size_t ei = 0; ei < mesh3->m_elements().size(); ei++)
	{
		//Вариант m->GetEField(ei, eFields::int_ds) -> g_MuInf.dee
		mesh3->GetEField(ei, eFields::int_d) += mesh3->GetEField(ei, eFields::int_ds) * dt; // накопленная интенсивность деформаций
																							//m->GetEField(ei, eFields::int_d) = fabs(m_CurrentSqrSum-m_FirstSqrSum)/m_FirstSqrSum;
		
		mesh3->GetEField(ei, eFields::int_ss) = m_mat.Material().GetGrainSize(mesh3->GetEField(ei, eFields::int_ss), mesh3->GetEField(ei, eFields::int_ds), dt); //  Albakri only
		//ALOGI("SA1", AllToString(ei) + _T(" ") +AllToString(mesh3->GetEField(ei, eFields::int_ss)));
	}

	/*
	for(int i=0; i < m->m_bordernodes().GetCount(); i++){
	if(m_bc[i].getType() == C2DBCAtom::symX || m_bc[i].getType() == C2DBCAtom::symY){
	int n = m->m_bordernodes()[i];
	}
	} //*/

}


Math::C2DRect C2DPlaneFEM::GetBoundingBox()
{
	Math::C2DRect rect;

	/*rect.left = 0;
	rect.bottom = -5;
	rect.top = 45;
	rect.right = 60;
	return rect; //костыль
	*/
	if (m_mesh_adapt()->GetType() == CIO2DMESH3) {
		rect = dynamic_cast<C2DMesh3*>(m_mesh_adapt())->GetBoundingBox();
	}

	return rect;
}

bool C2DPlaneFEM::GetBoundingBox(CRect2D &rect)
{
	if (m_mesh_adapt()->GetType() == CIO2DMESH3) {
		rect = dynamic_cast<C2DMesh3*>(m_mesh_adapt())->GetBoundingBox();
		return true;
	}
	return false;
}

void C2DPlaneFEM::Preparations(const ITask *pTask) {

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());
	if (!mesh3) {
		//Нет КЭ сетки, ничего не можем больше
		CDlgShowError cError(ID_ERROR_2DPLANEFEM_NO_FE_MESH); //_T("Нет КЭ сетки \n"));
		return;
	}

	//Освобождаем все граничные узлы
	size_t bc_size = mesh3->m_bordernodes.GetSize();
	for (size_t i = 0; i < bc_size; i++) {
		m_bc[i].freePoint();
	}

	//Проходим по границе и смотрим, как на неё влияют все объекты
	size_t nObjNum = pTask->GetObjNum();
	for (size_t j = 0; j<nObjNum; j++) {
		I2DInteractiveObject *obj2D = dynamic_cast<I2DInteractiveObject *>(pTask->GetObj(j));
		if (obj2D && obj2D != this) {

			std::vector<C2DBCAtom> temp_bc(bc_size);
			//Суммируем граничные условия со всех инструментов\объектов
			if (obj2D->GetBC(mesh3, &temp_bc)) {
				for (size_t k = 0; k < bc_size; k++) {
					m_bc[k] += temp_bc[k];
				}
			}
		}
	}
}

void C2DPlaneFEM::DrawGL(GLParam &parameter)
{
	//Чтобы записать в лог значения полей (в отладке)
	//LOGGER.Init(CString(_T("..\\..\\Logs\\C2DPlaneFEM.cpp_DrawGL.txt")));
	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	//! Отрисовка шкалы
	if (!mesh3) {
		return;
	}
	mesh3->DrawGL(parameter);

	//! Отрисовка сетки
	GLParam p;
	p.SetField(-1);
	mesh3->DrawGL(p); // Граница

	//! Отображение обычных узлов сетки
	/*glPointSize(3);
	for (int i = 0; i < mesh3->m_nodes().GetCount(); i++) {
		glColor3d(0.5, 0.5, 0.5);
		glBegin(GL_POINTS);
		glVertex2d(mesh3->m_nodes()[i].x, mesh3->m_nodes()[i].y);
		glEnd();
	}*/



	 //Отрисовка граничных условий
	 double pixelsize = 1;// валюнтаризм, по идее это должно передаваться в DrawGL()
	 double max_force=0;
	 size_t bordernodes_count = mesh3->m_bordernodes().size();
	 size_t mbc_count = m_bc.GetSize();

	 for(size_t i=0; i < min(bordernodes_count, mbc_count); i++){
		if(m_bc[i].getType() == C2DBCAtom::load){
			double force = sqrt(m_bc[i].getQx()*m_bc[i].getQx() + m_bc[i].getQy()*m_bc[i].getQy());
			if(force>max_force)max_force = force;
		}
	 }

	 double kff = 10 * pixelsize;
	 glPointSize(5);
	 for (size_t i = 0; i < min(bordernodes_count, mbc_count); i++) {
		 Math::C2DPoint& node = mesh3->m_nodes()[mesh3->m_bordernodes()[i]];
		 switch (m_bc[i].getType()) {
		 case C2DBCAtom::free:
			 glColor3d(0.7, 0.7, 0.7);
			 break;
		 case C2DBCAtom::symY:
			 glColor3d(0, 1, 0);	//green
			 break;
		 case C2DBCAtom::symX:
			 glColor3d(0, 0, 1);	//blue
			 break;
		 case C2DBCAtom::kinematic:
			 glColor3d(1, 0, 0);	//red
			 break;
		 case C2DBCAtom::load:
			 glColor3d(0, 0, 0);
			 glBegin(GL_LINES);
			 glVertex2d(node.x, node.y);
			 glVertex2d(node.x - kff*m_bc[i].getQx(), node.y - kff*m_bc[i].getQy());
			 glEnd();
			 break;
		 default:
			 glColor3d(1, 1, 1);
			 break;
		 }
		 glBegin(GL_POINTS);
		 glVertex2d(node.x, node.y);
		 glEnd();
	 }//*/

	 /* уже готово, но придётся перестраивать сборки - не залил в svn
	 //Отрисовка горизонтальных рисок
	 for (int i=0; i < m->m_bordernodes_marks().GetCount(); i++){
	 Math::C2DPoint& node = m->m_nodes()[m->m_bordernodes_marks()[i]];
	 glColor3d(0.1,0.2,0.5);
	 glPointSize(3);
	 glLineWidth(3);
	 glBegin(GL_LINES);
	 glVertex2d(node.x, node.y);
	 glVertex2d(node.x + 3.0, node.y);
	 glEnd();
	 } //*/

	 //WriteEFieldsToLog();

	glPointSize(1);
}


void C2DPlaneFEM::DrawGL3D(GLParam &parameter)
{

	DBL j, z, z1;
	DBL  step;
	//float m_Arg=270;             //значени угла без учета шага
	//float arg= m_Arg - step;    //значение угла c учетом шага для полигонов

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());
	if (!mesh3) { return; }

	DBL angle = parameter.GetAngle();
	int n = parameter.GetStep3D();
	step = angle / DBL(n);



	//! Отрисовка шкалы
	mesh3->DrawGL3D(parameter);

	//! Отрисовка сетки
	GLParam p;
	//m->DrawGL(p);//Сетка
	p.SetField(-1);
	mesh3->DrawGL3D(p); // Граница боковая 1

	glMatrixMode(GL_MODELVIEW);
	//glPushMatrix();
	//glLoadIdentity();
	GLfloat gl_angle = float(angle);
	GLfloat gl_x = 0.0;
	GLfloat gl_y = 1.0;
	GLfloat gl_z = 0.0;
	glRotatef(gl_angle, gl_x, gl_y, gl_z);
	mesh3->DrawGL3D(parameter);
	mesh3->DrawGL3D(p);  // Граница боковая 2
	//glPopMatrix();


	auto& bnodes = mesh3->m_bordernodes(); //граничные узлы



	glLineWidth(1);
	size_t bordernodes_count = mesh3->m_bordernodes().size();

	std::vector<C2DVertex> bverts;
	bverts.clear();
	parameter.m_PolygonBuilder.UpdateBorderNodes(*mesh3, parameter.GetField(), *(parameter.GetScale()), bverts);

	//parameter.m_PolygonBuilder.SetScaleBounds(*mesh3, parameter.GetField(), *(parameter.GetScale()));


	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	for (size_t i = 0; i <bverts.size() - 1; i++) {
		Math::C2DPoint& node = bverts[i]; //m->m_nodes()[m->m_bordernodes()[i]];
		Math::C2DPoint& node1 = bverts[(i + 1)]; /*% bverts.size()*/ //m->m_nodes()[m->m_bordernodes()[i+1]];

		//DBL p1=mesh3->GetNField(m->m_bordernodes()[i], parameter.GetField());
		//устанавливаем цвет 
		C2DColor color = parameter.GetScale()->GetColor(bverts[i].m_value);

		for (j = 0;j<n;j++) {
			glColor3d(color.GetR(), color.GetG(), color.GetB());

			z = j*step * M_PI_180;
			z1 = ((j + 1)*step) * M_PI_180;
			//glColor3d(0.7,0.7,0.7);

			//подключаем z буфер	
			if (parameter.IsFillColor()) {
				//2DColor color = parameter.GetScale()->GetColor(m->GetNField(m->m_bordernodes()[i], parameter.GetField()));
				//	glColor3d(color.GetR(), color.GetG(), color.GetB());
				glBegin(GL_QUADS);              //отрисовываем полигоны  квадратами
				glVertex3d((node.x)*cos(z), (node.y), (node.x)*sin(z));
				glVertex3d((node1.x)*cos(z), (node1.y), (node1.x)*sin(z));
				glVertex3d((node1.x)*cos(z1), (node1.y), (node1.x)*sin(z1));
				glVertex3d((node.x)*cos(z1), (node.y), (node.x)*sin(z1));

				glEnd();

				if (!bverts[i].level)
					continue;

				glLineWidth(0.7f);
				glColor3d(0, 0, 0);
				glBegin(GL_LINES);  //отрисовываем черным линии уровня //стыки между полигонами
				glVertex3d((node.x)*cos(z), (node.y), (node.x)*sin(z));
				glVertex3d((node.x)*cos(z1), (node.y), (node.x)*sin(z1));
				glEnd();


			}
			else {

				glColor3d(0.7, 0.7, 0.7);
				glBegin(GL_QUADS);              //отрисовываем полигоны  квадратами
				glVertex3d((node.x)*cos(z), (node.y), (node.x)*sin(z));
				glVertex3d((node1.x)*cos(z), (node1.y), (node1.x)*sin(z));
				glVertex3d((node1.x)*cos(z1), (node1.y), (node1.x)*sin(z1));
				glVertex3d((node.x)*cos(z1), (node.y), (node.x)*sin(z1));

				glEnd();

				/*if (!bverts[i].level)
				continue;*/

				glLineWidth(0.7f);
				glColor3d(0, 0, 0);
				glBegin(GL_LINE_STRIP);  //отрисовываем черным линии уровня //стыки между полигонами
				glVertex3d((node.x)*cos(z), (node.y), (node.x)*sin(z));
				glVertex3d((node1.x)*cos(z), (node1.y), (node1.x)*sin(z));
				glVertex3d((node1.x)*cos(z1), (node1.y), (node1.x)*sin(z1));
				glVertex3d((node.x)*cos(z1), (node.y), (node.x)*sin(z1));
				glEnd();
			}

		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	//}
	/* glCullFace (GL_FRONT);
	glPolygonMode (GL_BACK, GL_FILL);
	glDepthFunc (GL_LEQUAL);*/

	//}


	//	
	//for (int i=0; i<m->m_elements().GetCount(); i++){
	//	
	//	Math::C2DPoint p0 = m->m_nodes()[m->m_elements()[i].n0],
	//				p1 = m->m_nodes()[m->m_elements()[i].n1],
	//				p2 = m->m_nodes()[m->m_elements()[i].n2]; 
	//
	//	for(j=arg;60<=j;j-=10){
	//		if(j==arg || j ==60){
	//			z= j * M_PI_180;
	//			glColor3d(0.7,0.7,0.7);
	//		     glBegin(GL_TRIANGLES);
	//		      glVertex3f((p0.x)*cos(z), (p0.y),(p0.x)*sin(z));
	//		        glVertex3f((p1.x)*cos(z), (p1.y),(p1.x)*sin(z));
	//				  glVertex3f((p2.x)*cos(z), (p2.y),(p2.x)*sin(z));
	//			  
	//		     glEnd();
	//		}
	//	       
	//		}
	//
	//	
	//	
	//	}
	//	glLineWidth(3);
	//	INT_PTR bordernodes_count = m->m_bordernodes().GetCount();
	//	
	//	for(INT_PTR i=1; i <bordernodes_count; i++){
	//		Math::C2DPoint& node = m->m_nodes()[m->m_bordernodes()[i]];
	//			for(j=arg;60<=j;j--){
	//			//glPolygonMode(GL_FRONT,GL_LINE);
	//		    //glPolygonMode(GL_BACK,GL_FILL);
	//
	//			z= j * M_PI_180;
	//			k=(j+1)* M_PI_180;
	//			glColor3d(0.7,0.7,0.7);
	//			glBegin(GL_QUAD_STRIP);
	//		    glVertex3d( (m->m_nodes()[m->m_bordernodes()[i]].x),  (m->m_nodes()[m->m_bordernodes()[i]].y),0);
	//			glVertex3d((m->m_nodes()[m->m_bordernodes()[i-1]].x), (m->m_nodes()[m->m_bordernodes()[i-1]].y),0);
	//
	//		    glVertex3d((m->m_nodes()[m->m_bordernodes()[i]].x)*cos(z),(m->m_nodes()[m->m_bordernodes()[i]].y),(m->m_nodes()[m->m_bordernodes()[i]].x)*sin(z));
	//			glVertex3d((m->m_nodes()[m->m_bordernodes()[i-1]].x)*cos(z), (m->m_nodes()[m->m_bordernodes()[i-1]].y),(m->m_nodes()[m->m_bordernodes()[i-1]].x)*sin(z));
	//			glEnd();
	//
	//			//glBegin(GL_POLYGON);
	//		   // glVertex3d( (m->m_nodes()[m->m_bordernodes()[i]].x),  (m->m_nodes()[m->m_bordernodes()[i]].y),0);
	//			
	//			//glVertex3d((m->m_nodes()[m->m_bordernodes()[i-1]].x)*cos(z), (m->m_nodes()[m->m_bordernodes()[i-1]].y),(m->m_nodes()[m->m_bordernodes()[i-1]].x)*sin(z));
	//		   // glVertex3d((m->m_nodes()[m->m_bordernodes()[i]].x)*cos(z),(m->m_nodes()[m->m_bordernodes()[i]].y),(m->m_nodes()[m->m_bordernodes()[i]].x)*sin(z));
	//			
	//			//glEnd();
	//		
	//	}  
	//	}
	//	
	//	
	//	/*
	//	for(int i=0; i < m->m_nodes().GetCount(); i++){
	//		for(j=arg;358<=j;j--){
	//				z= j * M_PI_180;
	//		    glColor3d(1.0,0.0,0.0);
	//		      glBegin(GL_POLYGON);
	//		        glVertex3d((m->m_nodes()[i].x)*cos(z), (m->m_nodes()[i].y),(m->m_nodes()[i].x)*sin(z));
	//		      glEnd();
	//			}
	//	}*/
	//
	//
	//		//! Отображение обычных узлов сетки
	//	/*glPointSize(30);
	//	for(int i=0; i < m->m_nodes().GetCount(); i++){
	//		glColor3d(0.0,1.0,0.0);
	//		glBegin(GL_POINTS);
	//		glVertex2d(m->m_nodes()[i].x, m->m_nodes()[i].y);
	//		glEnd();
	//		
	//	}//*/
	//
	//	/*for(int i=0; i < m->m_nodes().GetCount(); i++){
	//			for(j=arg;359<=j;j--){
	//				z= j * M_PI_180;
	//			     glBegin(GL_POLYGON);
	//				 glColor3d(0.7,0.7,0.7);
	//			      glVertex3f((m->m_nodes()[i].x)*cos(z), (m->m_nodes()[i].y),(m->m_nodes()[i].x)*sin(z));
	//				
	//			     // glVertex3f((m->m_nodes()[i].x)*cos(z), (m_nodes[m_elements[i].n1].m_y),(m_nodes[m_elements[i].n1].m_x)*sin(z));
	//			     // glVertex3f((m_nodes[m_elements[i].n2].m_x)*cos(z), (m_nodes[m_elements[i].n2].m_y),(m_nodes[m_elements[i].n2].m_x)*sin(z));
	//				  
	//			     glEnd();
	//
	//		       
	//			}
	//	}
	//	for(int i=0; i < m->m_nodes().GetCount(); i++){
	//			for(j=arg;359<=j;j--){
	//				z= j * M_PI_180;
	//			     glBegin(GL_POLYGON);
	//				 glColor3d(0.7,0.7,0.7);
	//			      glVertex3f((m->m_nodes()[i].x)*cos(z), (m->m_nodes()[i].y),(m->m_nodes()[i].x)*sin(z));
	//				
	//			     // glVertex3f((m->m_nodes()[i].x)*cos(z), (m_nodes[m_elements[i].n1].m_y),(m_nodes[m_elements[i].n1].m_x)*sin(z));
	//			     // glVertex3f((m_nodes[m_elements[i].n2].m_x)*cos(z), (m_nodes[m_elements[i].n2].m_y),(m_nodes[m_elements[i].n2].m_x)*sin(z));
	//				  
	//			     glEnd();
	//
	//		       
	//			}
	//	}
	//	
	//	//Отрисовка граничных условий
	//	double pixelsize = 1;// валюнтаризм, по идее это должно передаваться в DrawGL()
	//	double max_force=0;
	////	INT_PTR bordernodes_count = m->m_bordernodes().GetCount();
	//
	//	
	//
	//	for(INT_PTR i=0; i < min(bordernodes_count, m_bc.GetSize()); i++){
	//		if(m_bc[i].getType() == C2DBCAtom::load){
	//			double force = sqrt(m_bc[i].getQx()*m_bc[i].getQx() + m_bc[i].getQy()*m_bc[i].getQy());
	//			if(force>max_force)max_force = force;
	//		}
	//	}
	//
	//	double kff = 10*pixelsize;
	//	glPointSize(5);
	//	for(INT_PTR i=0; i < min(bordernodes_count, m_bc.GetSize()); i++){
	//		Math::C2DPoint& node = m->m_nodes()[m->m_bordernodes()[i]];
	//		switch(m_bc[i].getType()){
	//		case C2DBCAtom::free:
	//			glColor3d(0.7,0.7,0.7);
	//			break;
	//		case C2DBCAtom::symY:
	//			glColor3d(0,1,0);
	//			glPointSize(1);
	//			break;
	//		case C2DBCAtom::symX:
	//			glColor3d(0,0,1);
	//			break;
	//		case C2DBCAtom::kinematic:
	//			glColor3d(1,0,0);
	//			break; 
	//		case C2DBCAtom::load:
	//			glColor3d(0,0,0);
	//			glBegin(GL_LINE_STRIP);
	//
	//			glVertex3d(node.x, node.y,node.x);
	//			glVertex3d(node.x - kff*m_bc[i].getQx(), node.y - kff*m_bc[i].getQy(),node.x - kff*m_bc[i].getQx());
	//			
	//			glEnd();
	//			break; 
	//		default:
	//			glColor3d(1,1,1);
	//			break;
	//
	//			/*case C2DBCAtom::load:
	//			glColor3d(0.7,0.7,0.7);
	//			glBegin(GL_POLYGON);
	//
	//			
	//			glVertex3f((m->m_nodes()[i].x)*cos(z), (m->m_nodes()[i].y),(m->m_nodes()[i].x)*sin(z));
	//			glEnd();
	//			break; */
	//		}
	//		// Отрисовываем покрытие каркаса из профилей заготовки
	//
	//
	//
	//		for(j=arg;90<=j;j--){
	//			//glPolygonMode(GL_FRONT,GL_LINE);
	//		   // glPolygonMode(GL_BACK,GL_FILL);
	//			z= j * M_PI_180;
	//			k=(j)* M_PI_180;
	//			
	//			glBegin(GL_POLYGON);
	//		    glVertex3d((node.x), node.y,0);
	//			glVertex3d((node.x+1 )*cos(z), node.y,(node.x+1)*sin(z));
	//		//	glVertex3d((node.x+1), node.y,0);
	//		    glVertex3d((node.x)*cos(z), node.y,(node.x)*sin(z));
	//			
	//			glEnd();
	//
	//		
	//			/*for(int i=0; i < m->m_nodes().GetCount(); i++){
	//		    // glColor3d(1.0,0.0,0.0);
	//		      glBegin(GL_POLYGON);
	//		        glVertex3d((m->m_nodes()[i].x)*cos(z), (m->m_nodes()[i].y),(m->m_nodes()[i].x)*sin(z));
	//		      glEnd();
	//			}
	//	}     
	//	} 
	//	
	//	// уже готово, но придётся перестраивать сборки - не залил в svn
	//	//Отрисовка горизонтальных рисок
	//	for (int i=0; i < m->m_bordernodes_marks().GetCount(); i++){
	//		Math::C2DPoint& node = m->m_nodes()[m->m_bordernodes_marks()[i]];
	//		glColor3d(0.1,0.2,0.5);
	//		glPointSize(3);
	//		glLineWidth(3);
	//		glBegin(GL_LINES);
	//		glVertex2d(node.x, node.y);
	//		glVertex2d(node.x + 3.0, node.y);
	//		glEnd();
	//	} 
	//	
	//	WriteEFieldsToLog();
	//	//glPointSize(1);
	//}
	//	
}

//! Записываем ГУ в лог
void C2DPlaneFEM::WriteBCToLog()
{
	DLOG(CString(_T("C2DPlaneFEM<Type:")) + AllToString(GetType()) + CString(_T(">")), log_info);

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());

	CString sTmp = CString(_T(""));

	DLOG(CString(_T("Writing border nodes:")), log_info);
	for (size_t i = 0; i < mesh3->m_bordernodes().size(); i++) {

		int k = mesh3->m_bordernodes()[i];
		sTmp = CString(_T("Node[")) + AllToString(k) + CString(_T("]: ")) + AllToString(mesh3->m_nodes()[k].x) + CString(_T(", ")) + AllToString(mesh3->m_nodes()[k].y);
		DLOG(sTmp, log_info);
	}

	/*
	DLOG(CString(_T("Writing all nodes:")), log_info);
	for (int i=0; i<m->m_nodes().GetCount(); i++){
	Math::C2DPoint& n = m->m_nodes()[i];

	DLOG(AllToString(n._0) + CString(_T(", ")) + AllToString(n._1), log_info);
	} //*/

	for (size_t i = 0; i < m_bc.GetSize(); i++)
	{
		sTmp = CString(_T("BC[")) + AllToString(i) + CString(_T("]: type:")) + AllToString(m_bc[i].getType()) +
			CString(_T(", Qx:")) + AllToString(m_bc[i].getQx()) + CString(_T(", Qy:")) + AllToString(m_bc[i].getQy()) +
			CString(_T(", a:")) + AllToString(m_bc[i].getAngle());
		DLOG(sTmp, log_info);
	}
	DLOG(CString(_T("End of border nodes")), log_info);
}

//! Записываем поля элементов в лог для каждого шага
void C2DPlaneFEM::WriteEFieldsToLog() {

	//DLOG(CString(_T("C2DPlaneFEM<Type:")) + AllToString(GetType()) + CString(_T(">")), log_info);

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());
	//CString strEfield = CString(_T(""));
	//DLOG(CString(_T("Writing EFields:")), log_info);
	//DLOG(CString(_T("cm_x | cm_y | avg_d | smo_d | int_d | int_ds | exx | eyy | efi |")), log_info);

	//for(int i=0; i<m->m_elements().GetCount(); i++){

	//Случайный эл-нт
	size_t i = mesh3->m_elements().size() / 2;
	CString strEfield = AllToString(m_CurTime()) + CString(_T(" | "));

	strEfield += AllToString(mesh3->GetEField(i, cm_x)) + CString(_T(" | "));
	strEfield += AllToString(mesh3->GetEField(i, cm_y)) + CString(_T(" | "));
	//strEfield += AllToString(m->GetEField(i,def))+CString(_T(" | ")); //не используем нигде пока (e+66)
	//strEfield += AllToString(m->GetEField(i,def_s))+CString(_T(" | ")); //не используем нигде пока (e+66)
	//strEfield += AllToString(m->GetEField(i,str))+CString(_T(" | ")); //*/
	strEfield += AllToString(mesh3->GetEField(i, avg_d)) + CString(_T(" | "));
	strEfield += AllToString(mesh3->GetEField(i, avg_ds)) + CString(_T(" | "));
	//sTmp += AllToString(m->GetEField(i,eps))+CString(_T(" | ")); 
	//strEfield += AllToString(m->GetEField(i,nju))+CString(_T(" | ")); //не используем нигде пока
	strEfield += AllToString(mesh3->GetEField(i, mju)) + CString(_T(" | "));

	strEfield += AllToString(mesh3->GetEField(i, smo_d)) + CString(_T(" | "));
	strEfield += AllToString(mesh3->GetEField(i, smo_ds)) + CString(_T(" | "));
	strEfield += AllToString(mesh3->GetEField(i, int_d)) + CString(_T(" | "));
	strEfield += AllToString(mesh3->GetEField(i, int_ds)) + CString(_T(" | "));

	strEfield += AllToString(mesh3->GetEField(i, sf_a)) + CString(_T(" | "));
	strEfield += AllToString(mesh3->GetEField(i, sf_b)) + CString(_T(" | "));
	strEfield += AllToString(mesh3->GetEField(i, sf_c)) + CString(_T(" | "));

	//strEfield += AllToString(m->GetEField(i,int_s))+CString(_T(" | ")); //не используем пока (e+66)
	//sTmp += AllToString(m->GetEField(i,int_ss))+CString(_T(" | "));
	//sTmp += AllToString(m->GetEField(i,hyd_p))+CString(_T(" | ")); //*/
	//sTmp += AllToString(m->GetEField(i,sqr))+CString(_T(" | ")); 
	ALOGI("Elements", strEfield);

	//}
	//}

	CString strNfield;
	//DLOG(CString(_T(" x | y | avg_d | int_d | int_ds | int_s | F |")), log_info);
	//for(int i=0; i<m->m_nodes().GetCount(); i++){

	//Узел выбранного эл-та
	int j0 = mesh3->m_elements()[i].n0;
	strNfield = AllToString(m_CurTime()) + CString(_T(" | "));

	strNfield += AllToString(mesh3->m_nodes()[j0].x) + CString(_T(" | "));
	strNfield += AllToString(mesh3->m_nodes()[j0].y) + CString(_T(" | "));

	strNfield += AllToString(mesh3->GetNField(j0, avg_d)) + CString(_T(" | "));
	strNfield += AllToString(mesh3->GetNField(j0, avg_ds)) + CString(_T(" | "));

	strNfield += AllToString(mesh3->GetNField(j0, int_d)) + CString(_T(" | "));
	strNfield += AllToString(mesh3->GetNField(j0, int_ds)) + CString(_T(" | "));

	strNfield += AllToString(mesh3->GetNField(j0, sigma_x)) + CString(_T(" | "));
	strNfield += AllToString(mesh3->GetNField(j0, sigma_y)) + CString(_T(" | "));

	strNfield += AllToString(mesh3->GetNField(j0, int_s)) + CString(_T(" | "));

	//strNfield += AllToString(m->GetNField(i,sf_a))+CString(_T(" | "));
	ALOGI("Nodes", strNfield);
	//}

	//DLOG(CString(_T("End of EFields")), log_info);
	//LOGGER.Init(CString("..\\..\\Logs\\trash.txt"));
}

//! Записываем поля элементов в лог 
void C2DPlaneFEM::WriteEFieldsToLog(DBL dt) {

	//LOGGER.Init(CString(_T("..\\..\\Logs\\C2DPlaneFEM.cpp_Efileds_"))+AllToString(dt)+CString(_T(".txt")));
	DLOG(CString(_T("C2DPlaneFEM<Type:")) + AllToString(GetType()) + CString(_T(">")), log_info);

	C2DMesh3* mesh3 = dynamic_cast<C2DMesh3*>(m_mesh_adapt());
	CString sTmp = CString(_T(""));
	DLOG(CString(_T("Writing EFields:")), log_info);
	//DLOG(CString(_T("cm_x | cm_y | def | def_s | str | avg_d | mju | sf_a | sf_b | sf_c | int_d | int_ds | int_s | int_ss | hyp_d | smo_d")), log_info);
	DLOG(CString(_T("cm_x | cm_y | def_s | avg_d | mju | int_d | int_ds | smo_d")), log_info);

	for (size_t i = 0; i<mesh3->m_elements().size(); i++) {
		sTmp = AllToString(i) + CString(_T(": "));

		sTmp += AllToString(mesh3->GetEField(i, cm_x)) + CString(_T(" | "));
		sTmp += AllToString(mesh3->GetEField(i, cm_y)) + CString(_T(" | "));//*/
																			//sTmp += AllToString(m->GetEField(i,def))+CString(_T(" | ")); 
		sTmp += AllToString(mesh3->GetEField(i, def_s)) + CString(_T(" | "));
		//sTmp += AllToString(m->GetEField(i,str))+CString(_T(" | ")); //*/
		sTmp += AllToString(mesh3->GetEField(i, avg_d)) + CString(_T(" | "));
		//sTmp += AllToString(m->GetEField(i,eps))+CString(_T(" | ")); 

		//sTmp += AllToString(m->GetEField(i,nju))+CString(_T(" | "));
		sTmp += AllToString(mesh3->GetEField(i, mju)) + CString(_T(" | "));
		//sTmp += AllToString(m->GetEField(i,sf_a))+CString(_T(" | ")); 
		//sTmp += AllToString(m->GetEField(i,sf_b))+CString(_T(" | "));
		//sTmp += AllToString(m->GetEField(i,sf_c))+CString(_T(" | ")); //*/
		//sTmp += AllToString(m->GetEField(i,sqr))+CString(_T(" | ")); 
		sTmp += AllToString(mesh3->GetEField(i, int_d)) + CString(_T(" | "));
		sTmp += AllToString(mesh3->GetEField(i, int_ds)) + CString(_T(" | "));

		/*sTmp += AllToString(m->GetEField(i,int_s))+CString(_T(" | "));
		sTmp += AllToString(m->GetEField(i,int_ss))+CString(_T(" | ")); //*/
		//sTmp += AllToString(m->GetEField(i,hyd_p))+CString(_T(" | ")); 
		sTmp += AllToString(mesh3->GetEField(i, smo_d)) + CString(_T(" | "));
		DLOG(sTmp, log_info);
	}
	DLOG(CString(_T("End of EFields")), log_info);
}

C2DPlaneFEM::~C2DPlaneFEM()
{
}

