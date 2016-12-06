#pragma once

//#include <vector>
//#include <fstream>
//#include <string>
#include "..\..\Common\Factory.h"
#include "..\..\Common\ITask.h"

//!  Класс Менеджер Задач (CTaskManager).
/*!
  Класс Менеджер Задач (CTaskManager) хранит все задачи, позволяет загружать,
  сохранять и взаимодействовать напрямую с задачами, а также включает в себя 
  фабрику задач.
*/

class CTaskManager// : public CIOBase
{
protected:
	//CIOTaskPtrArray m_tasks;
	//CArray<ITask*> m_tasks;
	//typedef Factory<int, ITask> CTaskFactory;

	IO::CPtrArray<CIOTASKPTRARRAY> m_tasks;

public:
	//typedef CTaskFactory::Creator CCreator;

protected:
	//CTaskFactory m_factory;

public:
	CTaskManager();

	bool LoadTask(CStorage& file);
	bool SaveTask(size_t id, CStorage& file);

	int CreateTask(int type);
	
	ITask* GetTask(size_t id) const;

	size_t GetTaskCount() const;

	//bool RegisterTaskCreator(int type, CCreator creator);
	//bool UnRegisterTaskCreator(int type);

	~CTaskManager();
};
