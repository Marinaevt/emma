#pragma once
#include "emmadoc.h"
#include "..\..\Common\ITask.h"
#include "..\..\Common\2DPoint.h"
#include "Results.h"
#include "..\..\Common\Startup.h"

//Документ для вывода результатов расчётов
class CCalcStage :	public CEMMADoc
{
	GETTYPE (CCALCSTAGE)
public:
	IO::CString m_path;	// путь к загружаемому файлу
	IO::CBool m_IsLoaded;	// загружен ли файл

public:
	CCalcStage(void);
	CCalcStage(CRuntimeClass *pViewClass, size_t nID);
	virtual ~CCalcStage(void);

	virtual void DrawGL(GLParam &parameter);
	virtual void DrawGL3D(GLParam &parameter);

	virtual bool LoadPathToTask(const CString &strFileName);
	//! После выбора результата (click) загружаем файл в m_pTask
	virtual bool LoadTask();
	CString MakeName(const CString &strFileName);
	bool Set3DFlag(bool bFlag);
	bool Get3DFlag();

protected:
	ITask *m_pTask;

protected:
	virtual void Deactivate(){DeleteTask();}
	void DeleteTask();
	virtual CViewDescriptor *GetViewDescriptor();
	void RegisterMembers();	//Регистрирует все нужные объекты в данном классе

};

IODEF (CCalcStage)