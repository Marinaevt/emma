#pragma once
#include "../EmmaDoc.h"
#include "../Resource.h"
#include "../2DSpecimen.h"
#include "../DlgAddSubDoc.h"
#include "RRollingOperation.h"

//TODO: Класс для прокатки
class CRollingDoc : public CEMMADoc
{
	GETTYPE (CROLLINGDOC)
public:
	CRollingDoc(void);
	virtual ~CRollingDoc(void);

public:
	size_t InsertRollOperation(const CString &strName, ptrdiff_t nIndex = -1);
	virtual CEMMADoc* Add(HTREEITEM hCurItem);
	//! Добавляет в задачу новую операцию с параметрами, задаваемыми пользователем
	virtual CEMMADoc* Add();
	virtual bool CanRemove(HTREEITEM hCurItem);

	virtual bool CanAdd() { return true; }
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc) { return true; }
};

IODEF (CRollingDoc)
