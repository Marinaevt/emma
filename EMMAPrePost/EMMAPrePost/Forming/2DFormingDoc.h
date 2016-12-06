#pragma once
#include "../Gen2DDoc.h"

class C2DFormingDoc : public CGen2DDoc
{
	GETTYPE (C2D_FORMING_DOC)
public:
	C2DFormingDoc(void);
	~C2DFormingDoc(void);

	virtual CEMMADoc *Add();
	virtual bool CanAdd() {return true;}
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc){return true;}
};

IODEF (C2DFormingDoc)