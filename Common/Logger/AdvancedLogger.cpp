#include "StdAfx.h"
#include "AdvancedLogger.h"
#include "atltime.h"

//const CString CLog::defaultLogFile = CString(_T("..\\..\\Logs\\auto_log.txt"));
const CString CLog::defaultLogFile = CString(_T("auto_log.txt"));
const CString CAdvancedLogger::defaultID = CString(_T("auto"));

/*
**	Clog
*/

CLog::CLog(const CString& storageName, bool append, int filter, int timeStamp, bool typeStamp)
{
	m_bOpened = 0;
	Init(storageName, append, filter, timeStamp, typeStamp);
}

bool CLog::Init(const CString& storageName, bool append, int filter, int timeStamp, bool typeStamp)
{
	m_iFilter = filter;
	m_iTimeStamp = timeStamp;
	m_bTypeStamp = typeStamp;
	m_iShift = 0;
	m_bAppend = append;

	return OpenStorage(storageName);
}

bool CLog::OpenStorage(const CString& storageName)
{
	if (m_bOpened == 1)
		m_Storage.Close();

	if (m_bAppend){
		m_bOpened = m_Storage.Open(storageName, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeText | CStorage::shareDenyWrite | CStorage::modeNoTruncate);
		m_Storage.SeekToEnd();
	}
	else {
		m_bOpened = m_Storage.Open(storageName, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeText | CStorage::shareDenyWrite);
	}

	return (m_bOpened != 0);
}

void CLog::CloseStorage()
{
	if (m_bOpened == 1)
	{
		m_Storage.Flush();
		m_Storage.Close();
	}
	m_bOpened = 0;
}

int CLog::GetMessageFilter() const
{
	return m_iFilter;
}

void CLog::SetMessageFilter(int filter)
{
	m_iFilter = filter;
}

int CLog::GetTimeStampFilter() const
{
	return m_iTimeStamp;
}

void CLog::SetTimeStampFilter(int filter)
{
	m_iTimeStamp = filter;
}

bool CLog::GetTypeStamp() const
{
	return m_bTypeStamp;
}

void CLog::SetTypeStamp(bool value)
{
	m_bTypeStamp = value;
}

bool CLog::WriteMessage(const CString& message, int type)
{
	CString writeString;

	int temp;

	SYSTEMTIME sysTime;
	GetLocalTime(&sysTime);

	CTime atlTime(sysTime);

	if (m_bOpened == 0)
		return false;

	if (!(m_iFilter & type))
		return false;

	if (m_bTypeStamp)
	{
		if (type & log_info)
			writeString += _T("[I]");

		if (type & log_error)
			writeString += _T("[E]");

		if (type & log_critical)
			writeString += _T("[C]");
	}

	if (m_iTimeStamp)
		writeString += _T("[");

	if (m_iTimeStamp & log_year)
	{
		temp = atlTime.GetYear();

		if (temp < 10)
			writeString += _T("0");

		writeString += AllToString(temp);
	}

	if ((m_iTimeStamp & log_year) && (m_iTimeStamp & log_month))
		writeString += _T("-");

	if (m_iTimeStamp & log_month)
	{
		temp = atlTime.GetMonth();

		if (temp < 10)
			writeString += _T("0");

		writeString += AllToString(temp);	
	}

	if ((m_iTimeStamp & log_month) && (m_iTimeStamp & log_day))
		writeString += _T("-");

	if (m_iTimeStamp & log_day)
	{
		temp = atlTime.GetDay();

		if (temp < 10)
			writeString += _T("0");

		writeString += AllToString(temp);	
	}

	if ((m_iTimeStamp & log_time) && (m_iTimeStamp & log_date))
		writeString += _T(" ");

	if (m_iTimeStamp & log_hours)
	{
		temp = atlTime.GetHour();

		if (temp < 10)
			writeString += _T("0");

		writeString += AllToString(temp);	
	}

	if ((m_iTimeStamp & log_hours) && (m_iTimeStamp & log_minutes))
		writeString += _T(":");

	if (m_iTimeStamp & log_minutes)
	{
		temp = atlTime.GetMinute();

		if (temp < 10)
			writeString += _T("0");

		writeString += AllToString(temp);	
	}

	if ((m_iTimeStamp & log_minutes) && (m_iTimeStamp & log_seconds))
		writeString += _T(":");

	if (m_iTimeStamp & log_seconds)
	{
		temp = atlTime.GetSecond();

		if (temp < 10)
			writeString += _T("0");

		writeString += AllToString(temp);	
	}

	if (m_iTimeStamp)
		writeString += _T("] ");

	// добавление сдвигов
	for (int i = 0; i < m_iShift; i++)
		writeString += _T("  ");

	writeString += message + _T("\r\n");

	try { 
		m_Storage.WriteString(writeString); 
		m_Storage.Flush();
	}
	catch(CException* pEx) {
		CDlgShowError CError(ID_ERROR_ADVANCEDLOGGER_WRITE_WRONG); //_T("AdvancedLogger WriteMessage Error"));
		pEx->Delete();
		return false; 
	}

		return true;
}

CLog::~CLog(void)
{
	CloseStorage();
}

/**/
///////////////////	Advanced Logger
/**/

CAdvancedLogger& CAdvancedLogger::operator = (const CAdvancedLogger& cadvlog) {
	//defaultID = cadvlog.defaultID;	//TODO: нет оператора присваивания
	return *this;
}

CLog* CAdvancedLogger::AddLog(const CString& storageName, const CString& ID, bool append, int filter, int timeStamp, bool typeStamp)
{
	if (checkData(ID))
		return GetLog(ID);

	CLog* pLog = new CLog(storageName, append, filter, timeStamp, typeStamp);

	if (registerData(ID, pLog))
		return pLog;

	return nullptr;
}

CLog* CAdvancedLogger::GetLog(const CString& ID)
{
	return getData(ID);
}

bool CAdvancedLogger::DelLog(const CString& ID)
{
	return unregisterData(ID);
}

CAdvancedLogger::IdVector CAdvancedLogger::GetAllIDs() const
{
	return getRegisteredIds();
}

bool CAdvancedLogger::WriteMessage(const CString& ID, const CString& message, int type)
{
	if (!checkData(ID))
		return false;

	CLog* pLog = getData(ID);
	if (!pLog) {
		return false;
	}
	return pLog->WriteMessage(message, type);
}

bool CAdvancedLogger::WriteOnDebug(const CString& ID, const CString& message, int type)
{
#ifdef _DEBUG
	return WriteMessage(ID, message, type);
#else
	return true;
#endif
}

CAdvancedLogger::~CAdvancedLogger(void)
{
	DataMap::const_iterator i = map_.begin();

	while (i != map_.end())
	{
		delete i->second;
		++i;
	}
}