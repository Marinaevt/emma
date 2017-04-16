// ChildFrm.cpp : implementation of the CChildFrame class
//

#include "stdafx.h"
#include "ChildFrm.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CChildFrame

IMPLEMENT_DYNCREATE(CChildFrame, CMDIChildWndEx)

BEGIN_MESSAGE_MAP(CChildFrame, CMDIChildWndEx)
	
	ON_WM_SETFOCUS()
	ON_WM_CREATE()
	ON_WM_MDIACTIVATE()
	ON_WM_CLOSE()

	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_ERASEBKGND()

END_MESSAGE_MAP()

// CChildFrame construction/destruction

CChildFrame::CChildFrame()
	: m_pView(nullptr)
{
	// TODO: add member initialization code here
}

CChildFrame::~CChildFrame()
{
	//if(m_pView) 
	delete m_pView;
}


int CChildFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying the CREATESTRUCT cs
	if( !CMDIChildWndEx::PreCreateWindow(cs) )
		return 0;

	cs.dwExStyle &= ~WS_EX_CLIENTEDGE;
	cs.lpszClass = AfxRegisterWndClass(0);
	return 1;
}

// CChildFrame diagnostics

#ifdef _DEBUG
void CChildFrame::AssertValid() const
{
	CMDIChildWndEx::AssertValid();
}

void CChildFrame::Dump(CDumpContext& dc) const
{
	CMDIChildWndEx::Dump(dc);
}
#endif //_DEBUG


int CChildFrame::OnCreate(LPCREATESTRUCT lpCreateStruct){
	if (CMDIChildWndEx::OnCreate(lpCreateStruct) == -1)
		return -1;
	return 0;
}

void CChildFrame::OnSetFocus(CWnd* pOldWnd) 
{
	CMDIChildWndEx::OnSetFocus(pOldWnd);

	if(m_pView)m_pView->SetFocus();
}

int CChildFrame::OnCmdMsg(UNINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo)
{
	if(nID == ID_FILE_NEW)return 0;
	// let the view have first crack at the command
	if (m_pView && m_pView->OnCmdMsg(nID, nCode, pExtra, pHandlerInfo))	return 1;
	
	// otherwise, do default handling
	return CMDIChildWndEx::OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
}

int CChildFrame::OnCreateClient(LPCREATESTRUCT lpcs, CCreateContext* pContext){
	// default create client will create a view if asked for it
	if (pContext != nullptr && pContext->m_pNewViewClass != nullptr)	{
		m_pView = (CBasicView *)CreateView(pContext, AFX_IDW_PANE_FIRST);
		if (m_pView == nullptr)return 0;
	}
	return 1;
}

void CChildFrame::OnMDIActivate(int bActivate, CWnd* pActivateWnd, CWnd* pDeactivateWnd) {
	CMDIChildWndEx::OnMDIActivate(bActivate, pActivateWnd, pDeactivateWnd);
	CMDIFrameWnd* pFrame = GetMDIFrame();
	if (bActivate) {
		pFrame->SendMessage(WMR_CHILDACTIVATE, WPARAM(this));
	}
	else if (pActivateWnd == nullptr) {
		pFrame->SendMessage(WM_ERASEBKGND);
		pFrame->SendMessage(WMR_CHILDACTIVATE);
	}
	m_pView->SendMessage(WMR_ACTIVATE, pActivateWnd == this);
}



void CChildFrame::OnClose(){
	CMDIFrameWnd* pFrame = GetMDIFrame();
	pFrame->SendMessage(WMR_CHILDCLOSE, WPARAM(this));
	CMDIChildWndEx::OnClose();
}

bool CChildFrame::ShowToolBar(UNINT nID){
	if(!m_pView)return false;
	return m_pView->ShowToolBar(nID);
}

bool CChildFrame::IsPaneVisible(UNINT nID){
	if (!m_pView) {
		return false;
	}
	CEMMADoc *pDoc = m_pView->GetDocument();
	if (!pDoc) {
		return false;
	}
	
	bool bPaneVisible = pDoc->IsPaneVisible(nID);	//всегда должны знать, что возвращает функция
	return bPaneVisible;
}

int CChildFrame::OnEraseBkgnd(CDC* pDC)
{
	return 1;
}

void CChildFrame::OnPaint()
{
	/*
	CPaintDC dc(this); // device context for painting
					   //Перекраска фона в белый цвет при изменении размера
	CRect rectClient;
	GetClientRect(&rectClient);

	CBrush bg;
	bg.CreateStockObject(WHITE_BRUSH);
	dc.FillRect(&rectClient, &bg);
	*/
	CPaintDC dc(this); // device context for painting

	CRect rectTree;
	GetClientRect(rectTree);

	dc.Draw3dRect(rectTree, ::GetSysColor(COLOR_3DSHADOW), ::GetSysColor(COLOR_3DSHADOW));
}

//////////////////////////////////////////////////////////////////////////////////////
void CChildFrame::OnSize(UNINT nType, int cx, int cy)
{
	CMDIChildWndEx::OnSize(nType, cx, cy);
	//AdjustLayout();
	Invalidate();
}