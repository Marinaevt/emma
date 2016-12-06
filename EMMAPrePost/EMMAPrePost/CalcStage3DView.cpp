#include "stdafx.h"
#include "CalcStage3DView.h"
#include "3DCalcStage.h"
#include "2DCalcStage.h"
#include "..\..\Common\3DPoint.h"
#include "..\..\Common\2DPoint.h"
#include "RegisteredMessages.h"
#include "ChildFrm.h"


IMPLEMENT_DYNCREATE(CCalcStage3DView, CGL3DView);

BEGIN_MESSAGE_MAP(CCalcStage3DView, CGL3DView)
	ON_REGISTERED_MESSAGE(WMR_DRAFTRECT, OnDraftRect)
	ON_COMMAND(ID_2DVIEW, &CCalcStage3DView::On2D)
	ON_REGISTERED_MESSAGE(WMR_ACTIVATE, OnActivate)
END_MESSAGE_MAP()

CCalcStage3DView::CCalcStage3DView(void)
{

	m_Parameter = new GLParam();
}

CCalcStage3DView::~CCalcStage3DView(void)
{
	delete m_Parameter;
	m_Parameter = nullptr;
}

void CCalcStage3DView::SetDocument(CEMMADoc *pDoc)
{
	

	CChildFrame *pChildFrame = dynamic_cast<CChildFrame *>(GetParent());
	CMDIFrameWnd* pFrame = pChildFrame->GetMDIFrame();
	pFrame->SendMessage(WMR_UPDATEPANE);


	CBasicView::SetDocument(pDoc);
}


//TODO: Установить прямоугольник
LRESULT CCalcStage3DView::OnDraftRect(WPARAM wParam, LPARAM lParam){
	CRect3D *Rect = (CRect3D *)wParam; // Rect3D or 2DRect ?
	//SetDraftRect( );       //   Какой параметр установить тут?                        (*m_DraftRect, lParam != 0);
	SetProjection();
	Invalidate();
	UpdateWindow();	

	return 0;
}


void CCalcStage3DView::On2D(){
	CCalcStage *pDoc = dynamic_cast<CCalcStage *>(GetDocument());
	if(!pDoc){
		return;
	}
	bool bFlagResult = pDoc->Set3DFlag(false);
	//GetParentFrame()->SendMessage(WMR_EP_SELCHANGE, (WPARAM)pDoc);
	GetTopLevelFrame()->SendMessage(WMR_EP_SELCHANGE, (WPARAM)pDoc);
}

LRESULT CCalcStage3DView::OnActivate(WPARAM wParam, LPARAM lParam){
	if(wParam){
		CCalcStage *pDoc = dynamic_cast<CCalcStage *>(GetDocument());
		if(!pDoc)return 0;
		pDoc->Set3DFlag(true);
	}
	return 0;
}

