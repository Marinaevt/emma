#pragma once
#include "../Gen2DDoc.h"
class C2DBalloonDoc :
	public CGen2DDoc
{

	GETTYPE(C2D_BALLOON_DOC)
public:
	C2DBalloonDoc(void);
	~C2DBalloonDoc(void);

	virtual CEMMADoc *Add();
	virtual bool CanAdd() { return true; }
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc) { return true; }
};

IODEF(C2DBalloonDoc)