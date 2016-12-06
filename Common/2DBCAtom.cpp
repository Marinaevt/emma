#include "StdAfx.h"
#include "2DBCAtom.h"


// [1] Qx = Px, Qy = Vy
// [2] Qx = Px, Qy = Py
// [3] Qx = Vx, Qy = Vy

void C2DBCAtom::setSymY(DBL dVelX, DBL dForceY) {
	a = 0;
	Qx = dVelX;
	Qy = dForceY;
	type = symY;
}

void C2DBCAtom::setSymX(DBL dVelY, DBL dForceX) {
	a = 0;
	Qx = dForceX;
	Qy = dVelY;
	type = symX;
}

void C2DBCAtom::setKinematic(DBL dVx, DBL dVy) {
	a = 0;
	Qx = dVx;
	Qy = dVy;
	type = kinematic;
}

void C2DBCAtom::setLoad(DBL dPx, DBL dPy) {
	a = 0;
	Qx = dPx;
	Qy = dPy;
	type = load;
}

void C2DBCAtom::freePoint() {
	a = 0;
	Qx = 0;
	Qy = 0;
	type = free;
}

C2DBCAtom C2DBCAtom::operator+ (const C2DBCAtom& bc) const {
	C2DBCAtom ret(bc);
	ret.Add(bc);
	return ret;
}

C2DBCAtom& C2DBCAtom::Add(const C2DBCAtom& bc)
{
	//Свободный узел не меняет объект
	if(bc.type == free) return *this;

	//Если сам объекст является свободным узлом. он принимает значение bc
	if(type == free){
		*this = bc;
		return *this;
	}

	// Обеспечиваем type > bc.getType(), чтобы избавиться от части if-ов
	if(type > bc.getType()){
		C2DBCAtom bc1(bc);
		*this = bc1.Add(*this);
		return *this;
	}

	//bc.type = symY // symY > free
	if (bc.type == symY){ // Оба - Симметрия относительно Y, угол игнорируется
		Qy += bc.Qy;
		Qx = 0.5*(Qx + bc.Qx);
		return *this;
	}
	
	//bc.type = symX // symX > symY > free
	if (bc.type == symX){
		if(type == symX){ // Оба - Симметрия относительно X, угол игнорируется
			Qx += bc.Qx;
			Qy = 0.5*(Qy + bc.Qy);
		}else{// Смешанные условия по разным осям = прилипание
			Qx = bc.Qx;
			type = kinematic;
		}
		return *this;
	}

	//bc.type = load // load > symX > symY > free
	if (bc.type == load){
		if(type == load){ // Оба - силовые. Просто складываем силы.
			Qx += bc.Qx;
			Qy += bc.Qy;
		}else if(type == symX){ // Сила + Симметрия относительно X, угол игнорируется
			Qx += bc.Qx*2;  //по оси X складываем силы, при этом сила bc должна умножаться на 2, чтобы учесть действие давления на площадку, половина которой находится за осью симметрии
			//Qy = bc.Qy;   //оставляем скорость по Y
		}else{					// Сила + Симметрия относительно Y, угол игнорируется
			Qy += bc.Qy*2;  //по оси Y складываем силы, при этом сила bc должна умножаться на 2, чтобы учесть действие давления на площадку, половина которой находится за осью симметрии
			//Qx = bc.Qx;   //оставляем скорость по X
		}
		return *this;
	}

	//bc.type = kinematic // kinematic больше всех остальных значений
	if (bc.type == kinematic){
		*this = bc;
		return *this;
	}
	
	CDlgShowError cError(ID_ERROR_2DBCATOM_WRONG_BC); //_T("Неверный тип ГУ \n"));
	return *this;
}

C2DBCAtom::~C2DBCAtom()
{
}
