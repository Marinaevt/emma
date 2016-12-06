#pragma once
#include "../Gen2Doperation.h"

class C2DFormingOperation :	public CGen2DOperation{
	GETTYPE (C2D_FORMING_OPERATION)
public:
	C2DFormingOperation(CString strName = 0);
	~C2DFormingOperation(void);

	virtual CEMMADoc *Add(){return nullptr;}
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc){return false;}
	virtual bool CanAdd() {return false;}
};

IODEF (C2DFormingOperation)