// MainFrm.cpp : implementation of the CMainFrame class
//

#include "stdafx.h"
#include "MainFrm.h"
#include "EMMASrv.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

//Регистрируем сообщение для "общения" с иконкой в Tray
UINT WM_NIF_MESSAGE = RegisterWindowMessage(_T("NIFF MESSAGE"));

// CMainFrame

IMPLEMENT_DYNAMIC(CMainFrame, CFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CFrameWnd)
	ON_WM_CREATE()
	ON_WM_SETFOCUS()
	ON_WM_DESTROY()
END_MESSAGE_MAP()

// CMainFrame construction/destruction

CMainFrame::CMainFrame()
{
	// TODO: add member initialization code here
}

CMainFrame::~CMainFrame()
{
	DestroyWindow();
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CFrameWnd::OnCreate(lpCreateStruct) == -1) {
	
		CDlgShowError cError(ID_ERROR_MAINFRM_ONCREATEFAIL);
		return -1;
	}
		
	// create a view to occupy the client area of the frame
	if (!m_wndView.Create(nullptr, nullptr, AFX_WS_DEFAULT_VIEW,
		CRect(0, 0, 0, 0), this, AFX_IDW_PANE_FIRST, nullptr))
	{
		
		CDlgShowError cError(ID_ERROR_MAINFRM_VIEW_FAIL);
		return -1;
	}

	SetTrayIcon(NIM_ADD);
	return 0;
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	if( !CFrameWnd::PreCreateWindow(cs) )
		return FALSE;
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	cs.dwExStyle &= ~WS_EX_CLIENTEDGE;
	cs.lpszClass = AfxRegisterWndClass(0);
	return TRUE;
}

// CMainFrame diagnostics

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CFrameWnd::Dump(dc);
}
#endif //_DEBUG


// CMainFrame message handlers

void CMainFrame::OnSetFocus(CWnd* /*pOldWnd*/)
{
	// forward focus to the view window
	m_wndView.SetFocus();
}

void CMainFrame::OnDestroy()
{
	DestroyWindow();
}

BOOL CMainFrame::OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo)
{
	// let the view have first crack at the command
	if (m_wndView.OnCmdMsg(nID, nCode, pExtra, pHandlerInfo))
		return TRUE;

	// otherwise, do default handling
	return CFrameWnd::OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
}

BOOL CMainFrame::SetTrayIcon(DWORD dwMessage){
	
	memset(&m_nid, 0, sizeof(m_nid)); 
	m_nid.cbSize = sizeof(m_nid);
	m_nid.hWnd = GetSafeHwnd(); // хэндл имеющегося окна 
	m_nid.uID = IDC_NOTIFY; // идентификатор иконки
	m_nid.hIcon = theApp.LoadIcon(IDR_MAINFRAME); // иконка
	lstrcpyn(m_nid.szTip, _T("EMMA computation server"), sizeof(m_nid.szTip)/sizeof(m_nid.szTip[0]));// текст подсказки копируется в структуру
	m_nid.uFlags = NIF_ICON | NIF_TIP | NIF_INFO | NIF_MESSAGE;
	m_nid.uCallbackMessage = WM_NIF_MESSAGE;//Сообщение для обработки
	return Shell_NotifyIcon(dwMessage, &m_nid);
}

//If the function succeeds, the return value is nonzero.
BOOL CMainFrame::DestroyWindow(){
	int ret = Shell_NotifyIcon(NIM_DELETE, &m_nid);	//удаление иконки //Returns 1 if successful, or 0 otherwise.
	m_wndView.DestroyWindow();
	return ret; //CFrameWnd::DestroyWindow();
}

LRESULT CMainFrame::WindowProc(UINT message, WPARAM wParam, LPARAM lParam){
	
	if(message == WM_NIF_MESSAGE && wParam == IDC_NOTIFY){
		switch(static_cast<UINT>(lParam)){
		case WM_LBUTTONDBLCLK: 
			ShowWindow(SW_SHOW);
			UpdateWindow();
			return 0;
		}
	}

	if (message == WM_DESTROY) {
		int ret = Shell_NotifyIcon(NIM_DELETE, &m_nid);	//удаление иконки
	}
	return CFrameWnd::WindowProc(message, wParam, lParam);
}
