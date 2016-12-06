#include "StdAfx.h"
#include "ITask.h"

/*
IInteractiveObject* СInteractiveObjectFactory::createObject(int type, int id, ITask* parent)
{
	if (checkData(type)) {
		IInteractiveObject* ptr = getData(type)(id, parent);
		return ptr;
	}
	return nullptr;
}

СInteractiveObjectFactory ITask::IntObjFactory;
*/

class CStepInfo{
public:
	CStepInfo(){
		t0 = 0.01;
		N = 100;
		max_dt = 1;
		MaxT = 1000;
		CStorage file;

		if(file.Open(_T("init.txt"), CStorage::modeRead | CStorage::typeBinary)){
			ReadValue(file, "<P>", &P);
			ReadValue(file, "<t0>", &t0);
			ReadValue(file, "<MaxDt>", &max_dt);
			double dVal;
			ReadValue(file, "<N>", &dVal);
			N = static_cast<int>(dVal);
			ReadValue(file, "<MinDt>", &MinDt);
			ReadValue(file, "<MaxT>", &MaxT);

			file.Close();
		}else{
			CDlgShowError cError(ID_ERROR_ITASK_INIT_NOT_FOUND); //_T("init.txt not found (ITask)"));
		}
	}
	~CStepInfo(){}

	int N;//Количество шагов для набора давления;
	double t0;//продолжительность набора давления;
	double max_dt;//максимальный шаг;
	double P;//давление;
	double MinDt;
	double MaxT;

};

CStepInfo g_Step;

IInteractiveObject* ITask::GetObj(size_t nIndex) const {

	if (nIndex >= GetObjNum()) {
		CDlgShowError cError(ID_ERROR_ITASK_NINDEX_WRONG); //CDlgShowError cError(_T("nIndex is wrong"));
		return nullptr;
	}

	return dynamic_cast<IInteractiveObject*>(m_objects()[nIndex]);
}

ITask::ITask()
{
	RegisterMember(&m_name);
	RegisterMember(&m_filename);
	RegisterMember(&m_filepath);
	RegisterMember(&m_dt);
	RegisterMember(&m_tStart);
	RegisterMember(&m_tEnd);
	RegisterMember(&m_tCurrent);
	RegisterMember(&m_objects);

	m_DtRecommended = 0.0;
}

ITask::~ITask()
{
	for (size_t i = 0; i < m_objects().size(); i++){
		if (m_objects()[i]) delete m_objects()[i];
	}
}

CString ITask::GetName() const
{
	return m_name();
}

// инициализация всех объектов
bool ITask::Init() {
	for (size_t i = 0; i<GetObjNum(); i++) {
		if (!GetObj(i)->Init()) {
			return false;
		}
	}
	return true;
}

DBL ITask::CalcDt()
{
	
	//if(m_tCurrent > g_Step.t0/g_Step.N*2) m_dt = m_DtRecommended;
	//if(m_tCurrent > g_Step.t0) m_dt = m_DtRecommended;
	//m_DtRecommended = m_dt;
	if(m_tCurrent > g_Step.t0){
		m_dt = min(m_tCurrent/g_Step.N, g_Step.max_dt);
		if(m_tCurrent > g_Step.MaxT){
			m_dt = g_Step.MinDt + (m_dt - g_Step.MinDt)*0.9;
		}
	}else{
		m_dt = g_Step.t0/g_Step.N;
	}//*/
	
	return m_dt;
}

void ITask::SetDt(DBL dt)
{
	m_dt() = dt;
	m_tCurrent() += dt;
}

bool ITask::CheckStop()   // проверка условия окончания расчета
{
	return (m_tCurrent() >= m_tEnd());
}

void ITask::AddObject(IInteractiveObject* obj)
{
	if (obj)
		m_objects().push_back(obj);
}

//! запуск вычислений
void ITask::RunCalc()
{
	int step = 0;
	int nLimiter = 5;
	CString fname;
	CStorage file;
	
	m_tStart = m_tCurrent = 0;
	m_dt = g_Step.t0/g_Step.N;
	m_tEnd = g_Step.MaxT;	//время остановки из файла

	/* //TODO: вытащить время всех траекторий
	for (int i = 0; i < m_objects().GetCount(); i++) {
		GetObj(i)->GetStopTime();
	} //*/

#ifndef __LOCAL
	ALOGGER.AddLog(m_filepath() + CString("Log_SA.txt"), CString(_T("SA")), false);
	ALOGGER.AddLog(m_filepath() + CString("Log_YD.txt"), CString(_T("YD")), false);

	ALOGGER.AddLog(m_filepath() + CString("Log_Elements.txt"), CString(_T("Elements")), false);
	ALOGGER.AddLog(m_filepath() + CString("Log_Nodes.txt"), CString(_T("Nodes")), false);
#else
	ALOGGER.AddLog(CString("Log_SA.txt"), CString(_T("SA")), false);
	ALOGGER.AddLog(CString("Log_YD.txt"), CString(_T("YD")), false);

	ALOGGER.AddLog(CString("Log_Elements.txt"), CString(_T("Elements")), false);
	ALOGGER.AddLog(CString("Log_Nodes.txt"), CString(_T("Nodes")), false);
#endif

	//ALOGI("SA", CString(_T("time, dt, H, s, e_min, e_max, ee_min, ee_max, s_min, s_max,   P=")) + AllToString(g_Step.P));
	ALOGI("YD", CString(_T("time, dt, force, avg_sigma, sqr, press, max_y, max_x, err, nK, dV, Sqr")));
	
	ALOGI("Elements", CString(_T("dt, cm_x, cm_y, avg_d, avg_ds, mju, smo_d, smo_ds, int_d, int_ds, sf_a, sf_b, sf_c")));
	ALOGI("Nodes", CString(_T("dt, x, y, avg_d, avg_ds, int_d, int_ds, sigma_x, sigma_y, int_s")));

	while (!CheckStop())
	{
		fname = m_filepath() + GenFileName(step, 5) + CString(_T("-")) + AllToString(m_tCurrent()) + CString(_T(".emma_task"));

		//if (step % nLimiter == 0) 
		file.Open(fname, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeBinary | CStorage::shareDenyWrite);

		// Двойной слеш уже есть в m_filepath()
		//CString newLogPath = m_filepath() + m_filename().Left(m_filename().GetLength() - 10) + CString("Logs\\");
		//CreateDirectory(newLogPath,nullptr);
		//LOGGER.Init(newLogPath+CString("RunCalc")+GenFileName(step,5)+ CString(_T("-")) + AllToString(m_tCurrent())+CString(_T(".txt")));

		Calculation(file);	//производим расчеты

		//if (step % nLimiter == 0) 
		LOGGER.Init(CString("..\\..\\Logs\\ITask.cpp_RunCalc_SaveFile.txt"));
		Save(file);	//сохраняем расчёты в файл
		file.Close();
		step++;
	}

}

//! расчёт задачи из файла
bool ITask::Calculation(CStorage& file)
{
	DBL dt = CalcDt();	//отрезок времени


	// подготовка к расчётам
	for (size_t i = 0; i < m_objects().size(); i++){
		GetObj(i)->Preparations(this);
	}

	// вычисления
	for (size_t i = 0; i < m_objects().size(); i++){
		GetObj(i)->Calc(dt);
	}
	
	SetDt(dt);
	// перемещения
	for (size_t i = 0; i < m_objects().size(); i++){
		GetObj(i)->Move(dt);
	}

	return true;
}

bool ITask::LoadBody(CStorage& file){
	if(!IO::CClass::LoadBody(file)) return false;

	for(size_t i=0; i<m_objects.GetSize(); i++){
		GetObj(i)->SetParent(this);
	}

	return true;
}

void ITask::DrawGL(GLParam &parameter) {
	for (size_t i = 0; i<m_objects().size(); i++) {
		(dynamic_cast<IInteractiveObject*>(m_objects()[i]))->DrawGL(parameter);
	}
}

void ITask::DrawGL3D(GLParam &parameter) {
	for (size_t i = 0; i<m_objects().size(); i++) {
		(dynamic_cast<IInteractiveObject*>(m_objects()[i]))->DrawGL3D(parameter);
	}
}

/*
bool ITask::Save(CStorage& file)
{
	bool ret = true;

	ret = ret && SaveFileIdAndVer(file);
	ret = ret && SaveTaskData(file);

	if (!ret)
		return false;

	for (size_t i = 0; i < m_objects.GetCount(); i++)
	{
		ret = ret && m_objects[i]->Save(file);
	}

	return ret;
}

bool ITask::Load(CStorage& file)
{
	bool ret = true;

	ret = ret && LoadName(file);
	
	if (!ret)
		return false;

	while (file.GetPosition() != CStorage::SeekPosition::end)
	{
		size_t type, id;

		file.Read((void *)&type, sizeof(type));
		file.Read((void *)&id, sizeof(id));

		IInteractiveObject* ptr = IntObjFactory.createObject(type, id, this);

		ret = ret && (ptr != nullptr);

		if (ptr)
		{
			ret = ret && ptr->Load(file);
			m_objects.Add(ptr);
		}
	}

	return ret;
}
*/