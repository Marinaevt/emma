#pragma once
#include "AllTypes.h"
#include "SimpleMath.h"
#include "DlgShowError.h"
#include "ResourceCommon.h"

class C3DBCAtom
{
public:
	DBL alpha; //<- 1 угол ЛСК относительно ГСК
	DBL beta; //<- 2 угол ЛСК относительно ГСК 
	
	DBL Qx, //<- скорость или давление по оси X в ЛСК
		Qy, //<- скорость или давление по оси Y в ЛСК
		Qz; //<- скорость или давление по оси Z в ЛСК

	TYPEID type; //<- тип граничного условия (1..3)
	// [1] Qx = Px, Qy = Vy, Qz = Pz
	// [2] Qx = Px, Qy = Py, Qz = Pz
	// [3] Qx = Vx, Qy = Vy, Qz = Vz - прилипание (т.е. движется вместе с границей)

	C3DBCAtom(TYPEID t = 0);

	C3DBCAtom& operator += (const C3DBCAtom& bc);
	C3DBCAtom operator + (const C3DBCAtom& bc) const;

	~C3DBCAtom();
};
