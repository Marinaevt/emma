#include "StdAfx.h"
#include "SimpleMath.h"

namespace Math {
	
	//Аппроксимация Zn по 4-м опорным значениям
	double APPR(double z1, double w1, double z2, double w2, double zn) {
		double dval;
		if (fabs(z1 - z2) < EPS) {
			dval = w2;
		}
		else {
			//не делим на 0
			dval = (w1*(z2 - zn) - w2*(z1 - zn)) / (z2 - z1);
		}
		return dval;
	}

}