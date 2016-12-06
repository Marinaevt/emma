#include "StdAfx.h"
#include "2DBuilderPane.h"


C2DBuilderPane::~C2DBuilderPane(void)
{
}

BEGIN_MESSAGE_MAP(C2DBuilderPane, CEMMARightPane)
	ON_WM_CREATE()
	ON_WM_SIZE()

	ON_REGISTERED_MESSAGE(AFX_WM_PROPERTY_CHANGED,C2DBuilderPane::OnPropertyChanged)
	
	////TOOLBAR

	ON_COMMAND(ID_BUILD_SAVE, &C2DBuilderPane::OnSave)
	ON_UPDATE_COMMAND_UI(ID_BUILD_SAVE, &C2DBuilderPane::OnUpdateSave)


END_MESSAGE_MAP()

int C2DBuilderPane::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CEMMADockablePane::OnCreate(lpCreateStruct) == -1)
		return -1;

	// Создание Тулбара {{
	m_wndToolBar.Create(this, AFX_DEFAULT_TOOLBAR_STYLE, IDR_2DBUILDER_TB);
	m_wndToolBar.LoadToolBar(IDR_2DBUILDER_TB, 0, 0, 1 /* Is locked */);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() | CBRS_TOOLTIPS | CBRS_FLYBY);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() & ~(CBRS_GRIPPER | CBRS_SIZE_DYNAMIC | CBRS_BORDER_TOP | CBRS_BORDER_BOTTOM | CBRS_BORDER_LEFT | CBRS_BORDER_RIGHT));
	
	m_wndToolBar.SetOwner(this);
	// All commands will be routed via this control , not via the parent frame:
	m_wndToolBar.SetRouteCommandsViaFrame(0);
	// }}


	CRect rectDummy;
	rectDummy.SetRectEmpty();

	if (!m_wndPropList.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, 1))
	{
		CDlgShowError cError(ID_ERROR_2DBUILDERPANE_FAIL_PROP_GRID); //_T("Failed to create Properties Grid \n"));
		return -1;      // fail to create
	}

	SetColumnNames();
	AdjustLayout();

	//m_cursor = 0;

	return 0;
}

void C2DBuilderPane::AdjustLayout(void){
	if (GetSafeHwnd() == nullptr)		return;

	CRect rectClient;
	GetClientRect(rectClient);

	int cyTlb = m_wndToolBar.CalcFixedLayout(0, 1).cy;

	m_wndToolBar.SetWindowPos(nullptr, rectClient.left, rectClient.top, rectClient.Width(), cyTlb, SWP_NOACTIVATE | SWP_NOZORDER);
	m_wndPropList.SetWindowPos(nullptr, rectClient.left, rectClient.top + cyTlb, rectClient.Width(), rectClient.Height() - cyTlb, SWP_NOACTIVATE | SWP_NOZORDER);
}

void C2DBuilderPane::OnSave(){

	C2DBuilder *pBuilder = dynamic_cast<C2DBuilder *>(m_pDoc);	//получаем ссылку на документ (Сборка)
	if (pBuilder){
		pBuilder->Build();
	}
}

void C2DBuilderPane::OnUpdateSave(CCmdUI *pCmdUI){

	//if(t_cursor == 9) pCmdUI->SetCheck(1);
	//else 
	pCmdUI->SetCheck(0);
}
