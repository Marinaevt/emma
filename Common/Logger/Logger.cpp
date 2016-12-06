#include "StdAfx.h"
#include "Logger.h"
#include "atltime.h"

bool CLogger::Init(const CString& storageName, int filter, int timeStamp, bool typeStamp)
{
	m_iFilter = filter;
	m_iTimeStamp = timeStamp;
	m_bTypeStamp = typeStamp;
	m_iShift = 0;

	return OpenStorage(storageName);
}

bool CLogger::OpenStorage(const CString& storageName)
{
	if (m_bOpened == 1)
		m_Storage.Close();

	m_bOpened = m_Storage.Open(storageName, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeText | CStorage::shareDenyWrite /*| CStorage::modeNoTruncate*/);
	//m_Storage.SeekToEnd();

	return (m_bOpened != 0);
}

int CLogger::GetMessageFilter() const
{
	return m_iFilter;
}

void CLogger::SetMessageFilter(int filter)
{
	m_iFilter = filter;
}

int CLogger::GetTimeStampFilter() const
{
	return m_iTimeStamp;
}

void CLogger::SetTimeStampFilter(int filter)
{
	m_iTimeStamp = filter;
}

bool CLogger::GetTypeStamp() const
{
	return m_bTypeStamp;
}

void CLogger::SetTypeStamp(bool value)
{
	m_bTypeStamp = value;
}


bool CLogger::WriteMessage(const CString& message, int type)
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
		CDlgShowError CError(ID_ERROR_LOGGER_WRITE_ERROR); //_T("CLogger WriteMessage Error"));
		pEx->Delete();
		return false; 
	}

	return true;
}

bool CLogger::WriteOnDebug(const CString& message, int type)
{
#ifdef _DEBUG
	return WriteMessage(message, type);
	//return true;
#else
	return true;
#endif
}

CLogger::CLogger(void)
{
	m_iFilter = log_info | log_error | log_critical;
	m_iTimeStamp = log_date | log_time;
	m_bOpened = 0;
	m_bTypeStamp = true;
	m_iShift = 0;
}

CLogger& CLogger::operator = (const CLogger& clog) {
	
	m_bOpened = clog.m_bOpened;
	m_bTypeStamp = clog.m_bTypeStamp;
	m_iFilter = clog.m_iFilter;
	m_iShift = clog.m_iShift;
	m_iTimeStamp = clog.m_iTimeStamp;
	//m_Storage = clog.m_Storage;	//TODO: оператор присваивания файлов
	return *this;
}

CLogger::~CLogger(void)
{
	if (m_bOpened == 1)
	{
		m_Storage.Flush();
		m_Storage.Close();
	}
}
