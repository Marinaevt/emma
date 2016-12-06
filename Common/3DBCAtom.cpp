#include "StdAfx.h"
#include "3DBCAtom.h"


// [1] Qx = Px, Qy = Vy, Qz = Pz
// [2] Qx = Px, Qy = Py, Qz = Pz
// [3] Qx = Vx, Qy = Vy, Qz = Vz

C3DBCAtom::C3DBCAtom(TYPEID t) : type(t) {
	alpha = 0.0;
	beta = 0.0;
	Qx = Qy = Qz = 0.0;
}

C3DBCAtom& C3DBCAtom::operator+= (const C3DBCAtom& bc)
{
	
	if (type == 0 || bc.type == 0) {
		CDlgShowError cError(ID_ERROR_3DBCATOM_BCTYPE_WRONG); //_T("Wrong BC type"));
		return *this;
	}

	if (type == 1 && bc.type == 1)
	{
		if (fabs(alpha - bc.alpha) < EPS && fabs(beta - bc.beta) < EPS) //TODO: здесь должно быть корректное сравнение углов
		{
			//type = 1;
			Qx += bc.Qx;
			Qy = max(Qy, bc.Qy);
			Qz += bc.Qz;
			return *this;
		}
		else
		{
			type = 3;
			Qx += bc.Qx * cos(alpha - bc.alpha);
			Qy += bc.Qy * sin(alpha - bc.alpha);
			//не совсем понял, что происходит
			return *this;
		}
	}
	if (type == 1 && bc.type == 2)
	{
		type = 1;
		Qx += bc.Qx * cos(alpha) + bc.Qy * sin(alpha);
		//ай не уверен
		Qz += bc.Qz * cos(beta) + bc.Qy * sin(beta);
		return *this;
	}
	if (type == 2 && bc.type == 1)
	{
		type = 1;
		Qx = bc.Qx + Qx * cos(bc.alpha) + Qy * sin(bc.alpha);
		Qy = bc.Qy;
		//ой не уверен
		Qz = bc.Qz + Qz * cos(bc.beta) + Qy * sin (bc.beta);
		return *this;
	}
	if (type == 3 && bc.type != 3)
	{
		return *this;
	}
	return *this;
}

C3DBCAtom C3DBCAtom::operator+ (const C3DBCAtom& bc) const
{
	return C3DBCAtom();
}

C3DBCAtom::~C3DBCAtom()
{
}