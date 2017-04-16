#include "StdAfx.h"
#include "2DToolPane.h"


C2DToolPane::~C2DToolPane(void)
{
}

BEGIN_MESSAGE_MAP(C2DToolPane, CEMMARightPane)
	ON_WM_CREATE()
	ON_WM_SIZE()

	//ON_REGISTERED_MESSAGE(AFX_WM_PROPERTY_CHANGED,C2DToolPane::OnPropertyChanged)
	
	// TOOLBAR от чертежа
	ON_COMMAND(ID_CLEAR_SKETCH, &C2DToolPane::OnClearSketch)
	ON_UPDATE_COMMAND_UI(ID_CLEAR_SKETCH, &C2DToolPane::OnUpdateClearSketch)

	ON_COMMAND(ID_NEW_NODE, &C2DToolPane::OnNewNode)
	ON_UPDATE_COMMAND_UI(ID_NEW_NODE, &C2DToolPane::OnUpdateNewNode)

	ON_COMMAND(ID_NEW_CURVE, &C2DToolPane::OnNewCurve)
	ON_UPDATE_COMMAND_UI(ID_NEW_CURVE, &C2DToolPane::OnUpdateNewCurve)

	ON_COMMAND(ID_NEW_CIRCLE_ARC, &C2DToolPane::OnNewCircleArc)
	ON_UPDATE_COMMAND_UI(ID_NEW_CIRCLE_ARC, &C2DToolPane::OnUpdateNewCircleArc)

	ON_COMMAND(ID_FACET, &C2DToolPane::OnNewFacet)
	ON_UPDATE_COMMAND_UI(ID_FACET, &C2DToolPane::OnUpdateNewFacet)

	//Точка привязки
	ON_COMMAND(ID_NEW_TOOL_NODE, &C2DToolPane::OnNewToolNode)
	ON_UPDATE_COMMAND_UI(ID_NEW_TOOL_NODE, &C2DToolPane::OnUpdateNewToolNode)

	//Buttons
	ON_BN_CLICKED(ID_RIGHTPANE_BUTTON_APPLY, C2DToolPane::ButtonApply)

END_MESSAGE_MAP()

int C2DToolPane::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CEMMADockablePane::OnCreate(lpCreateStruct) == -1)
		return -1;

	// Создание Тулбара {{
	m_wndToolBar.Create(this, AFX_DEFAULT_TOOLBAR_STYLE, IDR_2DTOOL_TB);
	m_wndToolBar.LoadToolBar(IDR_2DTOOL_TB, 0, 0, 1 /* Is locked */);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() | CBRS_TOOLTIPS | CBRS_FLYBY);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() & ~(CBRS_GRIPPER | CBRS_SIZE_DYNAMIC | CBRS_BORDER_TOP | CBRS_BORDER_BOTTOM | CBRS_BORDER_LEFT | CBRS_BORDER_RIGHT));
	
	m_wndToolBar.SetOwner(this);
	// All commands will be routed via this control , not via the parent frame:
	m_wndToolBar.SetRouteCommandsViaFrame(0);
	// }}

	// Создание Тулбара редактирования чертежа{{
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
		CDlgShowError cError(ID_ERROR_2DTOOLPANE_FAIL_PROPGRID); //_T("Failed to create Properties Grid \n"));
		return -1;      // fail to create
	}

	//creating Apply button for PropertyGrid
	int nButton = CreateApplyButton(rectDummy);
	if (nButton == -1) {
		return -1;	//тут же выходим
	}

	SetColumnNames();
	AdjustLayout();

	m_cursor = 0;

	return 0;
}

void C2DToolPane::AdjustLayout(void){
	if (GetSafeHwnd() == nullptr)		return;

	CRect rectClient;
	GetClientRect(rectClient);

	int cyTlb = m_wndToolBar.CalcFixedLayout(0, 1).cy;
	int cyTlb_sketch = m_wndToolBar_sketch.CalcFixedLayout(0, 1).cy;

	m_wndToolBar.SetWindowPos(nullptr, rectClient.left, rectClient.top, rectClient.Width(), cyTlb, SWP_NOACTIVATE | SWP_NOZORDER);
	m_wndToolBar_sketch.SetWindowPos(nullptr, rectClient.left, rectClient.top + cyTlb, rectClient.Width(), cyTlb_sketch, SWP_NOACTIVATE | SWP_NOZORDER);
	m_wndPropList.SetWindowPos(nullptr, rectClient.left, rectClient.top + cyTlb + cyTlb_sketch, rectClient.Width(), (rectClient.Height() - cyTlb - cyTlb_sketch) / 2, SWP_NOACTIVATE | SWP_NOZORDER);

	//устанавливаем размеры Кнопки
	int cyPropList = rectClient.top + cyTlb + cyTlb_sketch + (rectClient.Height() - cyTlb - cyTlb_sketch) / 2 + 10;	//чтобы не плотно +10
	int nFromLeft = rectClient.left + int(rectClient.Width() / 5.0);
	int nFromRight = rectClient.Width() - int(rectClient.Width() / 5.0) * 2;	//обязательное умножение, чтобы был отступ справа
	int nSetButton = m_buttonApply.SetWindowPos(nullptr, nFromLeft, cyPropList, nFromRight, 20, SWP_NOOWNERZORDER | SWP_SHOWWINDOW);
}

////////TOOLBAR

void C2DToolPane::OnNewToolNode(){
	
	if(m_cursor != 200){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 200);
		m_cursor = 200;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
	UpdatePane();
}

void C2DToolPane::OnUpdateNewToolNode(CCmdUI *pCmdUI){

	if(m_cursor == 200) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}


void C2DToolPane::OnClearSketch(){
	
	GetTool()->SetClearOutline();
	m_pDoc->GetView()->Invalidate();
	UpdatePane();
}

void C2DToolPane::OnUpdateClearSketch(CCmdUI *pCmdUI){
	// TODO: Add your command update UI handler code here
}

void C2DToolPane::OnNewNode(){
	
	if(m_cursor != 100){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 100);
		m_cursor = 100;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}

}
void C2DToolPane::OnUpdateNewNode(CCmdUI *pCmdUI){

	if(m_cursor == 100) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DToolPane::OnNewCurve(){
	
	if(m_cursor != 101){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 101);
		m_cursor = 101;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
}
void C2DToolPane::OnUpdateNewCurve(CCmdUI *pCmdUI){

	if(m_cursor == 101) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DToolPane::OnNewCircleArc(){
	
	if(m_cursor != 102){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 102);
		m_cursor = 102;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
}
void C2DToolPane::OnUpdateNewCircleArc(CCmdUI *pCmdUI){

	if(m_cursor == 102) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DToolPane::OnNewFacet(){
	
	if(m_cursor != 103){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 103);
		m_cursor = 103;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
}
void C2DToolPane::OnUpdateNewFacet(CCmdUI *pCmdUI){

	if(m_cursor == 103) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DToolPane::ButtonApply()
{
	m_pDoc->UpdatePropList(&m_wndPropList);
}