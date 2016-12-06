#pragma once
#include "Math.h"

class C3DElement4
{
public:
	union // индексы соседних (прилегающих) элементов
	{
		int n[4];
		struct { int _0, _1, _2, _3; };
		struct { int n0, n1, n2, n3; };
		struct { int first, second, third, fourth; };
	};

	double GetVolume(); 
};