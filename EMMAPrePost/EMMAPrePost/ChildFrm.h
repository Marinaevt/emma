// ChildFrm.h : interface of the CChildFrame class
//

#pragma once
#include "BasicView.h"
#include "EMMAPrePost.h"
#include "RegisteredMessages.h"

class CChildFrame : public CMDIChildWndEx
{
	DECLARE_DYNCREATE(CChildFrame)
public:
	CChildFrame();
	virtual ~CChildFrame();
// Attributes
public:

// Operations
public:
	CEMMADoc *GetDocument(){return m_pView ? m_pView->GetDocument() : nullptr;}
	bool ShowToolBar(UNINT nID);
	bool IsPaneVisible(UNINT nID);

// Overrides
public:
	virtual int PreCreateWindow(CREATESTRUCT& cs);
	virtual int OnCmdMsg(UNINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo);

// Implementation
public:
	// view for the client area of the frame.
	CBasicView *m_pView;

public:
	HTREEITEM GetTreeItem(){return m_pView ? m_pView->GetTreeItem() : nullptr;}
	DECLARE_MESSAGE_MAP()

//message map functions
protected:
	afx_msg void OnSetFocus(CWnd* pOldWnd);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	virtual int OnCreateClient(LPCREATESTRUCT lpcs, CCreateContext* pContext);
public:
	afx_msg void OnMDIActivate(int bActivate, CWnd* pActivateWnd, CWnd* pDeactivateWnd);
	afx_msg void OnClose();



#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	afx_msg void OnPaint();
	afx_msg int OnEraseBkgnd(CDC* pDC);
	afx_msg void OnSize(UNINT nType, int cx, int cy);
};
