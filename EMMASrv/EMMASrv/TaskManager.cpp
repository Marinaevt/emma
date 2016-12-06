#include "StdAfx.h"
#include "TaskManager.h"

CTaskManager::CTaskManager()
{
}

int CTaskManager::CreateTask(int type)
{
	/*
	CIOManager<ITaskHeader> &manager = CSingleton<CIOManager<ITaskHeader> >::Instance();

	ITaskHeader header;
	header.m_type = type;
	
	ITask* ptr = (ITask*)manager.getData(type)((CIOHeader*)&header); // .createObject(type);

	if (!ptr)
		return -1;
	*/

	IO::CInterface *ptr = IO::CreateObject(type); // надо как-то отделять типы задач от остальных типов
	
	if (!ptr)
		return -1;

	m_tasks().push_back(ptr);

	return (m_tasks().size() - 1);
}

bool CTaskManager::LoadTask(CStorage& file)
{
	// тут нужно открыть файл, прочитать тип задачи, попросить у фабрики создать задачу такого типа,
	// добавить полученный указатель в m_tasks, вызвать у созданной задачи функцию загрузки,
	// вернуть результат загрузки
	/*
	bool ret;

	CIOManager<ITaskHeader> &manager = CSingleton<CIOManager<ITaskHeader> >::Instance();

	ITask* ptr = (ITask*)manager.Load(file);

	if (!ptr)
		return false;
	*/
	IO::CInterface *ptr = IO::LoadObject(file);

	m_tasks().push_back(ptr);

	return ptr->Load(file);
}

bool CTaskManager::SaveTask(size_t id, CStorage& file)
{

	if (id < m_tasks().size()) {
		return m_tasks()[id]->Save(file);
	}
	return false;
}

ITask* CTaskManager::GetTask(size_t id) const
{
	if (id < m_tasks().size()) {
		return dynamic_cast<ITask*>(m_tasks()[id]);
	}
	return nullptr;
}

size_t CTaskManager::GetTaskCount() const
{
	return m_tasks().size();
}
/*
bool CTaskManager::RegisterTaskCreator(int type, CCreator creator)
{
	return m_factory.registerData(type, creator);
}

bool CTaskManager::UnRegisterTaskCreator(int type)
{
	return m_factory.unregisterData(type);
}
*/
CTaskManager::~CTaskManager()
{
	for (size_t i = 0; i < m_tasks().size(); i++)
		if (m_tasks()[i]) {
			delete m_tasks()[i];
		}
}
