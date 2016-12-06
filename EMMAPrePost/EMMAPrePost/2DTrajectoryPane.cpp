#include "StdAfx.h"
#include "2DTrajectoryPane.h"


C2DTrajectoryPane::~C2DTrajectoryPane(void)
{
}

BEGIN_MESSAGE_MAP(C2DTrajectoryPane, CEMMARightPane)
	ON_WM_CREATE()
	ON_WM_SIZE()

	////Roma
	ON_REGISTERED_MESSAGE(AFX_WM_PROPERTY_CHANGED,C2DTrajectoryPane::OnPropertyChanged)
	
	
	////TOOLBAR
	ON_COMMAND(ID_CLEAR_SKETCH, &C2DTrajectoryPane::OnClearSketch)
	ON_UPDATE_COMMAND_UI(ID_CLEAR_SKETCH, &C2DTrajectoryPane::OnUpdateClearSketch)

	ON_COMMAND(ID_NEW_NODE, &C2DTrajectoryPane::OnNewNode)
	ON_UPDATE_COMMAND_UI(ID_NEW_NODE, &C2DTrajectoryPane::OnUpdateNewNode)
	
	ON_COMMAND(ID_NEW_CURVE, &C2DTrajectoryPane::OnNewCurve)
	ON_UPDATE_COMMAND_UI(ID_NEW_CURVE, &C2DTrajectoryPane::OnUpdateNewCurve)

	ON_COMMAND(ID_NEW_CIRCLE_ARC, &C2DTrajectoryPane::OnNewCircleArc)
	ON_UPDATE_COMMAND_UI(ID_NEW_CIRCLE_ARC, &C2DTrajectoryPane::OnUpdateNewCircleArc)

	ON_COMMAND(ID_FACET, &C2DTrajectoryPane::OnNewFacet)
	ON_UPDATE_COMMAND_UI(ID_FACET, &C2DTrajectoryPane::OnUpdateNewFacet)

	//Включение анимации движения
	ON_COMMAND(ID_PLAY, &C2DTrajectoryPane::OnPlay)
	ON_UPDATE_COMMAND_UI(ID_PLAY, &C2DTrajectoryPane::OnUpdatePlay)

END_MESSAGE_MAP()

int C2DTrajectoryPane::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CEMMADockablePane::OnCreate(lpCreateStruct) == -1)
		return -1;

	// Создание Тулбара {{
	m_wndToolBar.Create(this, AFX_DEFAULT_TOOLBAR_STYLE, IDR_2DTRAJECTORY_TB);
	m_wndToolBar.LoadToolBar(IDR_2DTRAJECTORY_TB, 0, 0, 1 /* Is locked */);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() | CBRS_TOOLTIPS | CBRS_FLYBY);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() & ~(CBRS_GRIPPER | CBRS_SIZE_DYNAMIC | CBRS_BORDER_TOP | CBRS_BORDER_BOTTOM | CBRS_BORDER_LEFT | CBRS_BORDER_RIGHT));
	
	m_wndToolBar.SetOwner(this);
	// All commands will be routed via this control , not via the parent frame:
	m_wndToolBar.SetRouteCommandsViaFrame(0);
	// }}

	// Создание Тулбара редактирования чертежа {{
	m_wndToolBar_sketch.Create(this, AFX_DEFAULT_TOOLBAR_STYLE, IDR_2DSKETCH_TB);
	m_wndToolBar_sketch.LoadToolBar(IDR_2DSKETCH_TB, 0, 0, 1 /* Is locked */);
	m_wndToolBar_sketch.SetPaneStyle(m_wndToolBar_sketch.GetPaneStyle() | CBRS_TOOLTIPS | CBRS_FLYBY);
	m_wndToolBar_sketch.SetPaneStyle(m_wndToolBar_sketch.GetPaneStyle() & ~(CBRS_GRIPPER | CBRS_SIZE_DYNAMIC | CBRS_BORDER_TOP | CBRS_BORDER_BOTTOM | CBRS_BORDER_LEFT | CBRS_BORDER_RIGHT));
	
	m_wndToolBar_sketch.SetOwner(this);
	// All commands will be routed via this control , not via the parent frame:
	m_wndToolBar_sketch.SetRouteCommandsViaFrame(0);
	// }}

	CRect rectDummy;
	rectDummy.SetRectEmpty();

	if (!m_wndPropList.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, 1))
	{
		CDlgShowError cError(ID_ERROR_2DTRAJECTPANE_FAIL_PROPGRID); //_T("Failed to create Properties Grid \n"));
		return -1;      // fail to create
	}

	SetColumnNames();
	AdjustLayout();

	m_cursor = 0;
	
	return 0;
}

void C2DTrajectoryPane::AdjustLayout(void){
	if (GetSafeHwnd() == nullptr)		return;

	CRect rectClient;
	GetClientRect(rectClient);

	int cyTlb = m_wndToolBar.CalcFixedLayout(0, 1).cy;
	int cyTlb_sketch = m_wndToolBar_sketch.CalcFixedLayout(0, 1).cy;

	m_wndToolBar.SetWindowPos(nullptr, rectClient.left, rectClient.top, rectClient.Width(), cyTlb, SWP_NOACTIVATE | SWP_NOZORDER);
	m_wndToolBar_sketch.SetWindowPos(nullptr, rectClient.left, rectClient.top + cyTlb, rectClient.Width(), cyTlb_sketch, SWP_NOACTIVATE | SWP_NOZORDER);
	m_wndPropList.SetWindowPos(nullptr, rectClient.left, rectClient.top + cyTlb + cyTlb_sketch, rectClient.Width(), rectClient.Height() - cyTlb - cyTlb_sketch, SWP_NOACTIVATE | SWP_NOZORDER);
}

////////TOOLBAR
void C2DTrajectoryPane::OnClearSketch(){

	GetTrajectory()->SetClearOutline();
	Math::C2DPoint *p = new Math::C2DPoint(0,0);
	GetTrajectory()->AddNode(p);
	
	m_pDoc->GetView()->Invalidate();
	UpdatePane();
}

void C2DTrajectoryPane::OnUpdateClearSketch(CCmdUI *pCmdUI){
	// TODO: Add your command update UI handler code here
}

void C2DTrajectoryPane::OnNewNode(){
	
	if(m_cursor != 100){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 100);
		m_cursor = 100;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
	m_pDoc->GetView()->Invalidate();
	UpdatePane();
}

void C2DTrajectoryPane::OnUpdateNewNode(CCmdUI *pCmdUI){

	if(m_cursor == 100) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);

}

void C2DTrajectoryPane::OnNewCurve(){
	
	if(m_cursor != 101){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 101);
		m_cursor = 101;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
	m_pDoc->GetView()->Invalidate();
	UpdatePane();
}
void C2DTrajectoryPane::OnUpdateNewCurve(CCmdUI *pCmdUI){
	
	if(m_cursor == 101) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);

}

void C2DTrajectoryPane::OnNewCircleArc(){
	
	if(m_cursor != 102){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 102);
		m_cursor = 102;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
	m_pDoc->GetView()->Invalidate();
	UpdatePane();

}
void C2DTrajectoryPane::OnUpdateNewCircleArc(CCmdUI *pCmdUI){
	
	if(m_cursor == 102) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);

}

void C2DTrajectoryPane::OnNewFacet(){
	
	if(m_cursor != 103){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 103);
		m_cursor = 103;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
}

void C2DTrajectoryPane::OnUpdateNewFacet(CCmdUI *pCmdUI){

	if(m_cursor == 103) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DTrajectoryPane::OnPlay(){
	
	C2DTrajectory* pTraj = GetTrajectory();
	if(!pTraj) return;

	pTraj->m_dCurTime = 0.0;

	if(!pTraj->m_fPlay){
		pTraj->m_fPlay = true;
		m_pDoc->GetView()->SetTimer(51,4,nullptr);
	}else{
		pTraj->m_fPlay = false;
		m_pDoc->GetView()->KillTimer(51);
	}
	//m_pDoc->GetView()->Invalidate(0);
	UpdatePane();
}

void C2DTrajectoryPane::OnUpdatePlay(CCmdUI *pCmdUI){
	
	C2DTrajectory* pTraj = GetTrajectory();

	if (pTraj) {
		if (GetTrajectory()->m_fPlay) {
			pCmdUI->SetCheck(1);
		}
	}
	else {
		pCmdUI->SetCheck(0);
	}
	//*/
}