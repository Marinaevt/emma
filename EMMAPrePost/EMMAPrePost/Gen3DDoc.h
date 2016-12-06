#pragma once
#include "EmmaDoc.h"
#include "Resource.h"
#include "Child3DView.h"

class CGen3DDoc :	public CEMMADoc
{
	GETTYPE (CGEN3DDOC)
public:
	CGen3DDoc(void);
	CGen3DDoc(int Icon);
	virtual ~CGen3DDoc(void);
};

IODEF (CGen3DDoc)