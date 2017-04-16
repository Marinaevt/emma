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

	DECLARE_MESSAGE_MAP()
	afx_msg void OnPaint();
	afx_msg void OnSize(UNINT nType, int cx, int cy);
	afx_msg int OnEraseBkgnd(CDC* pDC);

protected:
	UNINT m_nResID;
	CEMMADoc *m_pDoc;
};
