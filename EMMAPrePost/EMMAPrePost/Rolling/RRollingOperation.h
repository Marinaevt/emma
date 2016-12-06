#pragma once
#include "../emmadoc.h"
#include "Grooves.h"
#include "../Resource.h"
#include "../2DResults.h"
#include "../2DBuilder.h"


class CRRollingOperation : public CEMMADoc
{
	GETTYPE (CR_ROLLING_OPERATION)
public:
	CRRollingOperation(void);
	virtual ~CRRollingOperation(void);
};

IODEF (CRRollingOperation)