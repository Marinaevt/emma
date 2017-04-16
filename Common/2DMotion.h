#pragma once

#include "2DOutline.h"
#include "Function.h"

/*!
  Узел траектории C2DMotionNode предназначен для использования классом C2DMotion.
*/
class C2DMotionParams
{
public:
	DBL m_angle, //< угол поворота объекта
		m_time,  //< время прибытия объекта в этот узел
		m_pause; //< время простоя объекта в этом узле

	DBL GetEntrTime() const { return m_time; }
	DBL GetExitTime() const { return m_time + m_pause; }
	DBL GetPause()    const { return m_pause; }

	DBL GetAngle()	  const { return m_angle; }

	//Время "входа", "выхода" и паузы в точке траектории
	void SetEntrTime(DBL t) { m_time = t; }
	void SetExitTime(DBL t) { if (t > m_time) m_pause = t - m_time; }
	void SetPause(DBL t)    { m_pause = t; }

	void SetAngle(DBL a)	{ m_angle = a; }
};


class C2DPosition
{
public:
	Math::C2DPoint m_pos, //< положение объекта
				   m_vel; //< вектор скорости объекта
	
	DBL m_angvel; //< угловая скорость объекта
	DBL m_angle;  //< угол поворота объекта
};


//! Двумерная траектория движения C2DMotion состоит из узлов C2DMotionNode, соединенных кривыми C2DMotionCurve.
class CIO2DMotion : public C2DOutline
{
	GETTYPE (C2DMOTION)
public:
	IO::CSimpleArray<C2DMotionParams, CIO2DMOTIONPARAMSARRAY> m_params;	// угол, время прибытия, пауза
	IO::CPtrArray<CIO2DMOTIONVELARRAY> m_vels; // функции скоростей: 0 - скорость движения, 1 - скорость вращения (угловая)

	CIO2DMotion()
	{
		RegisterMember(&m_params);
		RegisterMember(&m_vels);
	}

	// находит позицию на кривой #i в момент времени time
	bool GetPosOnCurve(size_t i, DBL time, C2DPosition& pos);

	// конвертирует время, которое прошел объект от узла i вдоль кривой i в расстояние вдоль этой кривой
	DBL TimeToLength(size_t i, DBL t); 

	// основная функция; возвращает положение объекта
	bool GetPos(DBL time, C2DPosition& pos); 

	// инициализирует параметры
	bool InitParams(size_t num); 

	// проверяет параметры траектории
	bool CheckParams(); 

	// считает время для каждой точки траектории
	void CalcTime(); 
	
	// получает мгновенную скорость в зависимости от расстояния и времени
	bool GetVel(DBL dLen, DBL &dVelx, DBL &dVely);

	// получает мгновенную скорость на кривой i в зависимости от расстояния и времени
	bool GetVelOnCurve(size_t i, DBL dLen, DBL &dVelx, DBL &dVely);
	bool GetVelOnCurveInPoint(size_t i, DBL dTime, DBL &dVelx, DBL &dVely);

	// задаёт всем узлам линейную скорость
	void SetLinVelToAll(DBL vel);

	virtual ~CIO2DMotion() {}
};

IODEF (CIO2DMotion)