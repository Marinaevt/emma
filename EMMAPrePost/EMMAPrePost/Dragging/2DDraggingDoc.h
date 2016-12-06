#pragma once
#include "../Gen2DDoc.h"

/*
	Заготовка для Волочения (протяжка). НЕ прокатка!
	Официальное название - Drawing, но это плохо, т.к. схоже с черчением\рисованием.
	Альтернативы - dragging, haul, traction, lug.
*/

class C2DDraggingDoc : public CGen2DDoc
{
	GETTYPE(C2D_DRAGGING_DOC)

public:
	
	C2DDraggingDoc(void);
	~C2DDraggingDoc(void);

	virtual CEMMADoc *Add();
	virtual bool CanAdd() { return true; }
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc) { return true; }
	
};

IODEF(C2DDraggingDoc)