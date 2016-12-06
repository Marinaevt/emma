#include "StdAfx.h"
#include "2DTask.h"

IOIMPL (C2DTask, T2DTASK)

bool C2DTask::Save(CStorage& file)
{
	return ITask::Save(file);
}

bool C2DTask::Load(CStorage& file)
{
	return ITask::Load(file);
}

//! запуск вычислений
void C2DTask::RunCalc()
{
	return ITask::RunCalc();
}

//! остановка вычислений
void C2DTask::StopCalc()
{
}

Math::C2DRect C2DTask::GetBoundingBox()
{
	Math::C2DRect rect;

	if (m_objects().size())
	{
		rect = dynamic_cast<I2DInteractiveObject*>(m_objects()[0])->GetBoundingBox();

		for (size_t i = 1; i < m_objects().size(); i++)
			rect.AddRect(dynamic_cast<I2DInteractiveObject*>(m_objects()[i])->GetBoundingBox());
	}

	return rect;
}

bool C2DTask::GetBoundingBox(Math::C2DRect &rect){
	if (m_objects().size())
	{
		for (size_t i = 0; i < m_objects().size(); i++){
			rect.AddRect(dynamic_cast<I2DInteractiveObject*>(m_objects()[i])->GetBoundingBox());
		}
		return true;
	}

	return false;
}

C2DTask::~C2DTask()
{
}