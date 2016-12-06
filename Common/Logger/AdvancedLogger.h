#pragma once

#include "../AllTypes.h"
#include "../Singleton.h"
#include "../SimpleManager.h"
#include "Logger.h"
#include "../ResourceCommon.h"
#include "../Strings.h"

#define ALOGGER CSAdvancedLogger::Instance()

// LOG - передает сообщение логгеру всегда, DLOG - только в режиме отладки
#define ALOG(id, m, t) ALOGGER.WriteMessage(CString(_T(id)), m, t)
#define DALOG(id, m, t) ALOGGER.WriteOnDebug(CString(_T(id)), m, t)
#define ALOGI(id, m) ALOGGER.WriteMessage(CString(_T(id)), m, log_info)
#define DALOGI(id, m) ALOGGER.WriteOnDebug(CString(_T(id)), m, log_info)

#define ALOGauto(m, t) ALOGGER.WriteMessage("auto", m, t)
#define DALOGauto(m, t) ALOGGER.WriteOnDebug("auto", m, t)

class CLog
{
	CStorage m_Storage; // файл, куда пишутся диагностические сообщения
	int m_bOpened; // файл открыт?

	int m_iFilter; // ниже какого уровня сообщения выводить

	int m_iTimeStamp; // в каком виде выводить время/дату; битовое множество
	bool m_bTypeStamp; // выводить ли тип сообщения

	int m_iShift; // текущий сдвиг в строке

	bool m_bAppend; // дописывать сообщения в конец файла (файл будет расти, пока не удалишь его вручную)

protected:
	bool OpenStorage(const CString& storageName);
	void CloseStorage();

public:
	static const CString defaultLogFile;

	CLog(const CString& storageName, bool append = true, int filter = log_alltypes, int timeStamp = log_datetime, bool typeStamp = true);

	bool Init(const CString& storageName, bool append = true, int filter = log_alltypes, int timeStamp = log_datetime, bool typeStamp = true);

	int GetMessageFilter() const;
	void SetMessageFilter(int filter);

	int GetTimeStampFilter() const;
	void SetTimeStampFilter(int filter);

	bool GetTypeStamp() const;
	void SetTypeStamp(bool value);

	void IncShift() { m_iShift += 1; }
	void DecShift() { m_iShift -= 1; }
	void ZeroShift() { m_iShift = 0; }

	int GetShift() { return m_iShift; }
	void SetShift(int s) { m_iShift = s; }

	bool GetAppend() { return m_bAppend; }

	bool WriteMessage(const CString& message, int type); // основная функция

	~CLog();
};

#define ADVLOGGER CSAdvancedLogger::Instance()

class CAdvancedLogger;

typedef CSingleton<CAdvancedLogger> CSAdvancedLogger;

class CAdvancedLogger : protected SimpleManager<CString, CLog*>
{
public:
	static const CString defaultID;

	CLog* AddLog(const CString& storageName, const CString& ID, bool append = true, int filter = log_alltypes, int timeStamp = log_datetime, bool typeStamp = true);
	CLog* GetLog(const CString& ID);
	bool DelLog(const CString& ID);
	IdVector GetAllIDs() const;

	bool WriteMessage(const CString& ID, const CString& message, int type); // основная функция
	bool WriteOnDebug(const CString& ID, const CString& message, int type);

	friend CSAdvancedLogger;

private:
	CAdvancedLogger(void)
	{
		AddLog(CLog::defaultLogFile, defaultID);
	}
	CAdvancedLogger(const CAdvancedLogger&);
	CAdvancedLogger& operator = (const CAdvancedLogger&);	//TODO: Доделать
	~CAdvancedLogger(void);
};