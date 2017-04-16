#pragma once

#include "AllTypes.h"
//#include "math.h"
#include "DlgShowError.h"
#include "ResourceCommon.h"

class C2DBCAtom
{
public:
	void setSymY(DBL dVelX = 0, DBL dForceY = 0, DBL angle = 0);

	void setSymX(DBL dVelY = 0, DBL dForceX = 0, DBL angle = 0);

	void setKinematic(DBL dVx, DBL dVy, DBL angle = 0);

	void setLoad(DBL dPx, DBL dPy, DBL angle = 0);

	void freePoint();
	
	C2DBCAtom(TYPEID t = free) : a(0.0), Qx(0.0), Qy(0.0), type(t) {}

	C2DBCAtom& operator += (const C2DBCAtom& bc){return Add(bc);}
	C2DBCAtom operator + (const C2DBCAtom& bc) const;

	C2DBCAtom& Add(const C2DBCAtom& bc);

	~C2DBCAtom();

	TYPEID getType() const {return type;}
	DBL getQx() const {return Qx;}
	DBL getQy() const {return Qy;}
	DBL getAngle() const {return a;}

protected:
	DBL a; //<- угол ЛСК относительно ГСК
	
	DBL Qx, //<- скорость или давление по оси X в ЛСК
		Qy; //<- скорость или давление по оси Y в ЛСК

	TYPEID type; //<- тип граничного условия (1..3)
public:
	enum{
		//НЕ МЕНЯТЬ ПОСЛЕДОВАТЕЛЬНОСТЬ ИДЕНТИФИКАТОРОВ! (используется в Add)
		free = 0,      // [0] свободная точка
		symY = 1,      // [1] Qx = Px, Qy = Vy, остальное не используется. Симметрия Y (скорость по одной оси, 0 давление по другой)
		symX = 2,      // [2] Qx = Vx, Qy = Py, остальное не используется. Симметрия X (скорость по одной оси, 0 давление по другой)
		load = 3,      // [3] Qx = Px, Qy = Py, Силы в узлах.
		kinematic = 100// [100] Qx = Vx, Qy = Vy, остальное не используется. Прилипание (т.е. движется вместе с границей)
	};	
};
