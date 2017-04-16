#include "StdAfx.h"
#include "EMMADockablePane.h"

IMPLEMENT_DYNAMIC(CEMMADockablePane, CDockablePane)

BEGIN_MESSAGE_MAP(CEMMADockablePane, CDockablePane)
	
	ON_WM_SIZE()
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()

END_MESSAGE_MAP()

CEMMADockablePane::CEMMADockablePane(size_t nID) : m_nResID(nID), m_pDoc(nullptr)
{
}

CEMMADockablePane::~CEMMADockablePane(void)
{
}

void CEMMADockablePane::SetDocument(CEMMADoc *pDoc){
	m_pDoc = pDoc;
}


int CEMMADockablePane::OnEraseBkgnd(CDC* pDC)
{
	int nEraseRightPaneBackground = CWnd::OnEraseBkgnd(pDC);
	return nEraseRightPaneBackground;
}

void CEMMADockablePane::OnPaint()
{
	CPaintDC dc(this); // device context for painting

	CRect rectTree;
	GetClientRect(&rectTree);	//получаем область панели
	dc.Draw3dRect(rectTree, ::GetSysColor(COLOR_3DSHADOW), ::GetSysColor(COLOR_3DSHADOW));	//закрашиваем серым цветом + объём

}

void CEMMADockablePane::OnSize(UNINT nType, int cx, int cy)
{
	CDockablePane::OnSize(nType, cx, cy);
	AdjustLayout();	//новые размеры элементов
	Invalidate();	//перерисовка вида с удалением фона
}