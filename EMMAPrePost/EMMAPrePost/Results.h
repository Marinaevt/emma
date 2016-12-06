#pragma once
#include "EmmaDoc.h"
#include "..\..\Common\FolderWork.h"
#include "RegisteredMessages.h"
#include "..\..\Common\Logger\Logger.h"


class CResults : public CEMMADoc
{
	GETTYPE (CRESULTS)
public:
	CResults(void);
	virtual ~CResults(void);

	size_t LoadResults();
	//! virtual LoadStage for 2D, 3D etc
	virtual bool LoadStage(const CString &strFileName){return false;}
	void Init2DChildViewDescriptor(CRuntimeClass *pViewClass, UNINT nID);
	void Init3DChildViewDescriptor(CRuntimeClass *pViewClass, UNINT nID);
	CViewDescriptor *GetChildViewDescriptor();

	CString MakeName(const CString &strFileName);

public:
	bool m_b3D;	//флаг 3D вида или 2D

protected:
	//Информация о представлении и фрейме поддокументов
	CViewDescriptor *m_p2DChildViewDescriptor;
	CViewDescriptor *m_p3DChildViewDescriptor;

	void DestroyViewDescriptor(CViewDescriptor *pDesctiptor);
};

IODEF (CResults)