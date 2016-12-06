#pragma once
#include "../Gen2Doperation.h"
class C2DBalloonOperation :
	public CGen2DOperation
{
	GETTYPE(C2D_BALLOON_OPERATION)
public:
	C2DBalloonOperation(CString strName = 0);
	~C2DBalloonOperation();

	virtual CEMMADoc *Add() { return nullptr; }
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc) { return false; }
	virtual bool CanAdd() { return false; }
};

IODEF(C2DBalloonOperation)