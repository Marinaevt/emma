#pragma once
#include "../EmmaDoc.h"
#include "../Resource.h"
#include "../2DSketch.h"

// Прокатка, валки (нижний и верхний)
class CGrooves : public CEMMADoc
{
	GETTYPE (CGROOVES)
public:
	CGrooves(void);
	virtual ~CGrooves(void);
};

IODEF (CGrooves)