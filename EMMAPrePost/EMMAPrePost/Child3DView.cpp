#include "StdAfx.h"
#include "Child3DView.h"

BEGIN_MESSAGE_MAP(CChild3DView, CGL3DView)
	ON_WM_LBUTTONUP()
END_MESSAGE_MAP()

CChild3DView::CChild3DView(void)
{
}

void CChild3DView::OnLButtonUp(UNINT nFlags, CPoint point)
{
	CGL3DView::OnLButtonUp(nFlags,point);
	CChildFrame *pChildFrame = dynamic_cast<CChildFrame *>(GetParent());
		CMDIFrameWnd* pFrame = pChildFrame->GetMDIFrame();
		pFrame->SendMessage(WMR_UPDATEPANE,(WPARAM) nFlags);
}
CChild3DView::~CChild3DView(void)
{
}


IMPLEMENT_DYNCREATE(CChild3DView, CGL3DView);