#pragma once

#include "afx.h" //<fstream>
#include "IOInterface.h"
#include "GLParam.h"

class ITask;

//!  Класс Взаимодействующий Объект (IInteractiveObject).
/*!
  Класс Взаимодействующий объект (IInteractiveObject) предоставляет интерфейс для иерархии
  всех взаимодействующих в процессе деформирования объектов - заготовок, штампов и др.
*/
class IInteractiveObject : public IO::CClass
{
protected:
	ITask *m_parent;

public:
	IInteractiveObject()
	{
		m_parent = nullptr; //parent;
	}
	
	void SetParent(ITask* parent)
	{
		m_parent = parent;
	}

	ITask* GetParent() const
	{
		return m_parent;
	}

	virtual TYPEID GetType() const = 0;
	
	virtual bool Init() = 0;
	virtual void Calc(DBL dt) = 0;
	virtual void Move(DBL dt) = 0;

	virtual void Preparations(const ITask *task) = 0;	// подготовка объектов для расчёта

	virtual void DrawGL(GLParam &parameter) = 0;
	virtual void DrawGL3D(GLParam &parameter) = 0;

	virtual ~IInteractiveObject()
	{}
};
