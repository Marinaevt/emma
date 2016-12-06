#pragma once
#include "afxdockablepane.h"
#include "EMMADoc.h"

class CEMMADockablePane : public CDockablePane{
	DECLARE_DYNAMIC(CEMMADockablePane)
public:
	CEMMADockablePane(size_t nID);
	virtual ~CEMMADockablePane(void);

	///////ROMA
	virtual void UpdatePane() =0;
	
	virtual void SetDocument(CEMMADoc *pDoc);
	UNINT GetResID(){return m_nResID;}
protected:
	UNINT m_nResID;
	CEMMADoc *m_pDoc;
};
