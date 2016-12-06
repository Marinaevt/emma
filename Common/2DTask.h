#pragma once

#include "ITask.h"
#include "I2DInteractiveObject.h"
#include "2DPoint.h"

class C2DTask :	public ITask
{
	GETTYPE (T2DTASK)
public:

	C2DTask()
	{
	}

	virtual void RunCalc();		// запуск вычислений
	virtual void StopCalc();	// остановка вычислений

	virtual bool Save(CStorage& file);
	virtual bool Load(CStorage& file);

	virtual Math::C2DRect GetBoundingBox();
	virtual bool GetBoundingBox(Math::C2DRect &rect);

	virtual ~C2DTask();
};

IODEF (C2DTask)