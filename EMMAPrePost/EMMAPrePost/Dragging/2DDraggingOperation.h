#pragma once
#include "../Gen2Doperation.h"

class C2DDraggingOperation : public CGen2DOperation {
	GETTYPE(C2D_DRAGGING_OPERATION)
public:
	C2DDraggingOperation(CString strName = _T(""));
	~C2DDraggingOperation(void);

	virtual CEMMADoc *Add() { return nullptr; }	//запрещаем добавлять операции
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc) { return false; }
	virtual bool CanAdd() { return false; }	
};

IODEF(C2DDraggingOperation)