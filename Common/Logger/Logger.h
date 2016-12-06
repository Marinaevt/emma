#pragma once

#include "../AllTypes.h"
#include "../Singleton.h"
#include "../IOFiles.h"
#include "../Strings.h"
#include "../ResourceCommon.h"

#define LOGGER CSLogger::Instance()

// LOG - передает сообщение логгеру всегда, DLOG - только в режиме отладки
#define LOG(m, t) LOGGER.WriteMessage(m, t)
#define DLOG(m, t) LOGGER.WriteOnDebug(m, t)

// SLOG = LOG, но только для строк в кавычках, которые сразу же преобразуются к CString
#define SLOG(m, t) LOGGER.WriteMessage(CString(_T(m)), t)
#define DSLOG(m, t) LOGGER.WriteOnDebug(CString(_T(m)), t)

enum eLogTypes// типы сообщений
{
	log_info     = static_cast<int>(0x01),
	log_error    = static_cast<int>(0x02),
	log_critical = static_cast<int>(0x04),
	log_alltypes = log_info | log_error | log_critical // (int)0x07
};

enum eLogDateTime// дата и время 
{
	log_year     = static_cast<int>(0x01),
	log_month    = static_cast<int>(0x02),
	log_day      = static_cast<int>(0x04),
	log_date     = log_year | log_month | log_day, // (int)0x07, // == 0x01 | 0x02 | 0x04
	log_hours    = static_cast<int>(0x08),
	log_minutes  = static_cast<int>(0x10),
	log_seconds  = static_cast<int>(0x20),
	log_time     = log_hours | log_minutes | log_seconds, // (int)0x38, // == 0x08 | 0x10 | 0x20
	log_datetime = log_date | log_time // (int)0x3f  // == 0x07 | 0x38
};

class CLogger;

typedef CSingleton<CLogger> CSLogger;

class CLogger
{
	CStorage m_Storage;	 // файл, куда пишутся диагностические сообщения
	int m_bOpened;		// файл открыт?

	int m_iFilter;		// ниже какого уровня сообщения выводить

	int m_iTimeStamp;	// в каком виде выводить время/дату; битовое множество
	bool m_bTypeStamp;	// выводить ли тип сообщения

	int m_iShift;		// текущий сдвиг в строке

public:

	bool Init(const CString& storageName, int filter = log_alltypes, int timeStamp = log_datetime, bool typeStamp = true);

	bool OpenStorage(const CString& storageName);

	int GetMessageFilter() const;
	void SetMessageFilter(int filter);

	int GetTimeStampFilter() const;
	void SetTimeStampFilter(int filter);

	bool GetTypeStamp() const;
	void SetTypeStamp(bool value);

	//сдвиг блока текста вправо (табуляция)
	void IncShift() { m_iShift += 1; }
	//сдвиг блока текста влево (табуляция)
	void DecShift() { m_iShift -= 1; }

	//int GetShift() { return m_iShift; }
	//void SetShift(int s) { m_iShift = s; }

	bool WriteMessage(const CString& message, int type); // основная функция
	bool WriteOnDebug(const CString& message, int type);

	friend CSLogger;

private:
	CLogger(void);
	CLogger(const CLogger&);
	CLogger& operator = (const CLogger&);	//TODO: нет оператора присваивания файлов
	~CLogger(void);
};
