#pragma once

#include "afxtempl.h"
#include "afx.h"
//#include "cstringt.h"
#include "FolderWork.h"
#include "./Logger/AdvancedLogger.h"
#include "InteractiveObject.h"
#include "SimpleManager.h"
#include "IOInterface.h"
#include "GLParam.h"
#include "DlgShowError.h"

//!  Класс Фабрика Взаимодействующих Объектов (СInteractiveObjectFactory).
/*!
  Класс Фабрика Взаимодействующих Объектов (СInteractiveObjectFactory) предназначен для
  создания объектов из иерархии Взаимодействующих Объектов с помощью вызова зарегистрированных
  в ней функций-создателей этих объектов.
*/
/*
class СInteractiveObjectFactory : public SimpleManager<int, IInteractiveObject* (*)(int, ITask*) >
{
	public:

	typedef IInteractiveObject* (*Creator)(int, ITask*);

	IInteractiveObject* createObject(int type, int id, ITask* parent);

	СInteractiveObjectFactory() 
	{}

	~СInteractiveObjectFactory()
	{}
};
*/


//!  Класс Задача (ITask).
/*!
  Класс Задача (ITask) предоставляет интерфейс для иерархии всех разновидностей задач
  деформирования - 2х-, 2.5- и 3х-мерных и др.
*/
class ITask : public IO::CClass
{
protected:
	IO::CString m_name;     // имя задачи
	IO::CString m_filename; // [имя_файла], 
	IO::CString m_filepath;	// [путь_к_файлу], при сохранении задачи автоматически будет добавлено [00000x.emma_task]

	IO::CDouble m_dt,       //< шаг по времени
			  m_tStart,   //< начало отсчета времени
			  m_tEnd,     //< конец отсчета времени
			  m_tCurrent; //< текущее время

	IO::CPtrArray<CIOINTPTRARRAY> m_objects;

public:
	//Aksenov{{
	size_t GetObjNum() const { return m_objects().size();}
	IInteractiveObject *GetObj(size_t nIndex) const;

	double GetCurrentT()const{return m_tCurrent();}
	double GetStart()const{return m_tStart();}	
	double GetElapsedTime()const{return m_tCurrent() - m_tStart();}

	double m_DtRecommended;
	void RecommendDt(double dt){m_DtRecommended = dt;}
	//}}Aksenov
		
public:
//	static СInteractiveObjectFactory IntObjFactory;

	ITask();
	
	CString GetName() const;

	void SetName(const CString& name){
		m_name() = name;
	}

	void SetFileName(const CString& fname){
		m_filename() = fname;
	}

	void SetFilePath(const CString& fpath){
		m_filepath() = fpath;
	}

	virtual TYPEID GetType() const = 0;

	// инициализация всех объектов
	virtual bool Init();
	
	virtual void RunCalc();		// запуск вычислений
	virtual void StopCalc() = 0;	// остановка вычислений
	
	virtual DBL CalcDt();       // вычисление дельта t для следующего шага
	virtual void SetDt(DBL dt); // учёт дельта t после выполнения шага
	virtual bool CheckStop();   // проверка условия окончания расчета

	virtual bool Calculation(CStorage& file);	// расчёт задачи из файла
	virtual bool LoadBody(CStorage& file);

	void AddObject(IInteractiveObject* obj);

	virtual void DrawGL(GLParam &parameter);
	virtual void DrawGL3D(GLParam &parameter);

	virtual ~ITask();
};

//typedef IO::CPtrArray<ITask, CIOTASKPTRARRAY> CIOTaskPtrArray;
