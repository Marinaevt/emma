#include "stdafx.h"
#include "ExplorerPane.h"


// CExplorerPane

IMPLEMENT_DYNAMIC(CExplorerPane, CEMMADockablePane)

CExplorerPane::CExplorerPane(size_t nID) :
	CEMMADockablePane(nID),
	m_bDontSendWM(false)
{
	if (!m_imageList.Create(16, 16, ILC_COLOR32, 1, 0)) {
		CDlgShowError cError(ID_ERROR_EXPLORERPANE_FAIL_IMAGELIST); //_T("Failed to create ImageList to Tree \n"));	// fail to create
	}
	m_icons_name = 0;
}

CExplorerPane::~CExplorerPane()
{
}


BEGIN_MESSAGE_MAP(CExplorerPane, CDockablePane)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()

	ON_NOTIFY(TVN_SELCHANGED, IDC_EP_TREE, &CExplorerPane::OnTvnSelchangedEPTree)

	//Добавляем подэлемент
	ON_COMMAND(ID_EXPLORER_ADD, &CExplorerPane::OnAdd)
	ON_UPDATE_COMMAND_UI(ID_EXPLORER_ADD, &CExplorerPane::OnUpdateAdd)

	//Удаляем подэлемент
	ON_COMMAND(ID_EXPLORER_REMOVE, &CExplorerPane::OnRemove)
	ON_UPDATE_COMMAND_UI(ID_EXPLORER_REMOVE, &CExplorerPane::OnUpdateRemove)

	//Смотрим результаты
	ON_COMMAND(ID_EXPLORER_RESULTS, &CExplorerPane::OnLoadResults)
	ON_UPDATE_COMMAND_UI(ID_EXPLORER_RESULTS, &CExplorerPane::OnUpdateResults)

END_MESSAGE_MAP()


//////////////////////////////////////////////////////////////////////////////////////
int CExplorerPane::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CDockablePane::OnCreate(lpCreateStruct) == -1)
		return -1;

	// Создание Тулбара {{
	m_wndToolBar.Create(this, AFX_DEFAULT_TOOLBAR_STYLE, IDR_EXPLORER_TB);
	m_wndToolBar.LoadToolBar(IDR_EXPLORER_TB, 0, 0, 1 /* Is locked */);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() | CBRS_TOOLTIPS | CBRS_FLYBY);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() & ~(CBRS_GRIPPER | CBRS_SIZE_DYNAMIC | CBRS_BORDER_TOP | CBRS_BORDER_BOTTOM | CBRS_BORDER_LEFT | CBRS_BORDER_RIGHT));
	m_wndToolBar.SetOwner(this);
	// If this parameter is TRUE, the parent frame sends commands to the toolbar.
	// Otherwise, the owner sends commands to the toolbar.
	m_wndToolBar.SetRouteCommandsViaFrame(0);
	// }}

	CRect rectDummy;
	rectDummy.SetRectEmpty();

	// Создание дерева {{
	const UNLONG dwTreeStyle = WS_CHILD | WS_VISIBLE | TVS_HASLINES | TVS_HASBUTTONS | TVS_SHOWSELALWAYS | TVS_FULLROWSELECT;

	if (!m_wndTree.Create(dwTreeStyle, rectDummy, this, IDC_EP_TREE)) {
		CDlgShowError cError(ID_ERROR_EXPLORERPANE_FAIL_TREE); //_T("Failed to create Explorer Tree \n"));
		return -1;      // fail to create
	}
	// }}

	AdjustLayout();

	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////
void CExplorerPane::AdjustLayout() {
	if (GetSafeHwnd() == nullptr) {
		return;
	}

	CRect rectClient;
	GetClientRect(rectClient);

	int cyTlb = m_wndToolBar.CalcFixedLayout(0, 1).cy;

	m_wndToolBar.SetWindowPos(nullptr, rectClient.left, rectClient.top, rectClient.Width(), cyTlb, SWP_NOACTIVATE | SWP_NOZORDER);
	m_wndTree.SetWindowPos(nullptr, rectClient.left + 1, rectClient.top + cyTlb + 1, rectClient.Width() - 2, rectClient.Height() - cyTlb - 2, SWP_NOACTIVATE | SWP_NOZORDER);
}

//////////////////////////////////////////////////////////////////////////////////////
CEMMADoc *CExplorerPane::GetCurDoc() {
	if (!m_wndTree) {
		return nullptr;	//Проверяем демонов наших, перед вызовом их...
	}

	HTREEITEM hItem = m_wndTree.GetSelectedItem();
	if (!hItem) {
		//TODO: show error
		//CDlgShowError cExplorerError(ID_???);
		return nullptr;
	}
	return (CEMMADoc *)m_wndTree.GetItemData(hItem);
}

//////////////////////////////////////////////////////////////////////////////////////
void CExplorerPane::OnSize(UNINT nType, int cx, int cy)
{
	CDockablePane::OnSize(nType, cx, cy);
	AdjustLayout();
	Invalidate();
}

//////////////////////////////////////////////////////////////////////////////////////
void CExplorerPane::OnPaint()
{
	CPaintDC dc(this); // device context for painting

	CRect rectTree;
	m_wndTree.GetWindowRect(rectTree);
	ScreenToClient(rectTree);

	dc.Draw3dRect(rectTree, ::GetSysColor(COLOR_3DSHADOW), ::GetSysColor(COLOR_3DSHADOW));
}

int CExplorerPane::OnEraseBkgnd(CDC * pDC)
{
	int nEraseRightPaneBackground = CWnd::OnEraseBkgnd(pDC);
	return nEraseRightPaneBackground;
}

//////////////////////////////////////////////////////////////////////////////////////
void CExplorerPane::FillTree() {
	if (!m_wndTree.GetSafeHwnd() || !m_pDoc)return;



	HTREEITEM hHeadItem = m_wndTree.GetNextItem(TVI_ROOT, TVGN_CHILD);

	//Очистка (сверху) ветвей дерева не соответствующих документу 
	while (hHeadItem != nullptr && hHeadItem != m_pDoc->GetTreeItem()) {
		m_wndTree.DeleteItem(hHeadItem);
		hHeadItem = m_wndTree.GetNextItem(TVI_ROOT, TVGN_CHILD);
	}

	if (hHeadItem == nullptr) {
		//добавление соответствующей документу ветки
		hHeadItem = m_pDoc->AddToTree(m_wndTree, TVI_ROOT, TVI_FIRST, m_icons_name);
	}
	else {

		//модификация соответствующей документу ветки
		HTREEITEM hTempItem = m_pDoc->ModifyTree(m_wndTree, m_icons_name);

		//Очистка (снизу) ветвей дерева не соответствующих документу 
		HTREEITEM hNextItem = m_wndTree.GetNextItem(hHeadItem, TVGN_NEXT);
		while (hNextItem != nullptr) {
			int nDelItem = m_wndTree.DeleteItem(hNextItem);
			hNextItem = m_wndTree.GetNextItem(hHeadItem, TVGN_NEXT);
		}

	}
	//m_wndTree.Expand(hHeadItem, TVE_EXPAND);
	m_wndTree.UpdateWindow();
	//m_wndTree.ShowWindow(SW_SHOW);
}

/*
по имени документа CString находим набор иконок и загружаем их
*/
void CExplorerPane::SetIconsToTree(CString str) {
	CString strTaskName;

	//список документов-задач (главных)
	std::vector<int> docIDs =
	{ IDS_GEN2D, IDS_FORMING,	IDS_ROLLING, IDS_DRAGGING, IDS_GEN3D, IDS_BALLOON };

	//список наборов иконок для этих задач (1:1)
	std::vector<int> imageListIDs =
	{ IDS_TREE_ICONS_GEN2D,	IDS_TREE_ICONS_FORMING,	IDS_TREE_ICONS_ROLLING,	IDS_TREE_ICONS_DRAGGING, IDS_TREE_ICONS_GEN3D, IDS_TREE_ICONS_FORMING };

	//по имени получаем набор иконок
	for (size_t i = 0; i < docIDs.size(); i++) {
		strTaskName.LoadString(docIDs[i]);
		if (str == strTaskName) { m_icons_name = imageListIDs[i]; }
	}

	//устанавливаем выбранный набор в дерево
	int resDelBitmap = m_bitmap.DeleteObject();
	int resLoadBitmap = m_bitmap.LoadBitmap(m_icons_name);

	//удаляем предыдущий загруженный набор иконок и пересоздаём
	m_imageList.DeleteImageList();
	if (!m_imageList.Create(16, 16, ILC_COLOR32, 1, 0)) {
		CDlgShowError cError(ID_ERROR_EXPLORERPANE_FAIL_IMAGELIST); //_T("Failed to create ImageList to Tree \n"));	// fail to create
	}
	int resAddToImageList = m_imageList.Add(&m_bitmap, RGB(0, 0, 0));	//Zero-based index of the first new image if successful; otherwise – 1.

	CImageList* pList = m_wndTree.SetImageList(&m_imageList, TVSIL_NORMAL);
}


void CExplorerPane::SetDocument(CEMMADoc *pDoc) {
	//устанавливаем документ
	CEMMADockablePane::SetDocument(pDoc);

	CString str = pDoc->GetName();
	//устанавливаем наборы иконок для дерева
	SetIconsToTree(str);

	//заполняем все поддокументы (вроде Чертежа, Материалов и т.д.)
	FillTree();
}

////////////////////////////////////////////////////////////////////////////////////
void CExplorerPane::SetCurSel(HTREEITEM hTreeItem) {
	if (!m_wndTree.GetSafeHwnd()) {
		return;
	}
	m_bDontSendWM = true;
	m_wndTree.SelectItem(hTreeItem);
	m_bDontSendWM = false;
}

////////////////////////////////////////////////////////////////////////////////////
void CExplorerPane::OnTvnSelchangedEPTree(NMHDR *pNMHDR, LRESULT *pResult) {
	if (m_bDontSendWM)return;

	LPNMTREEVIEW pNMTreeView = reinterpret_cast<LPNMTREEVIEW>(pNMHDR);
	*pResult = 0;

	CWnd* pFrame = GetParentFrame();
	HTREEITEM hItem = m_wndTree.GetSelectedItem();
	if (hItem) {
		CEMMADoc *pDoc = (CEMMADoc *)m_wndTree.GetItemData(hItem);
		pFrame->SendMessage(WMR_EP_SELCHANGE, (WPARAM)m_wndTree.GetItemData(hItem));
	}
}

void CExplorerPane::OnAdd() {
	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		CEMMADoc *pNewDoc = pCurDoc->Add();
		if (pNewDoc) {
			FillTree();
			m_wndTree.Expand(pNewDoc->GetTreeItem(), TVE_EXPAND);	//"раскрытие" у нового эл-та
			m_wndTree.UpdateWindow();
			// Активация созданного итема.
			m_wndTree.SelectItem(pNewDoc->GetTreeItem());
		}
	}
}

void CExplorerPane::OnUpdateAdd(CCmdUI *pCmdUI) {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCmdUI->Enable(GetCurDoc()->CanAdd());
	}
	else {
		pCmdUI->Enable(false);
	}
}

void CExplorerPane::OnRemove() {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCurDoc->Remove();
		FillTree();
		//m_wndTree.UpdateWindow();
	}
}

void CExplorerPane::OnRemoveAllSubDoc() {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		size_t nSub = pCurDoc->GetSubDocNum();

		while (pCurDoc->GetSubDocNum() > 0) {
			pCurDoc->SubDoc(0)->Remove();
		}

		FillTree();
		//m_wndTree.UpdateWindow();
	}
}

void CExplorerPane::OnUpdateRemove(CCmdUI *pCmdUI) {
	if (GetCurDoc()) {
		pCmdUI->Enable(GetCurDoc()->IsRemovable());
	}
	else {
		pCmdUI->Enable(false);
	}
}


void CExplorerPane::OnLoadResults() {

	if (GetCurDoc()) {

		CEMMADoc *pResDoc = GetCurDoc();
		C2DResults *pRes = dynamic_cast<C2DResults*>(GetCurDoc());
		CEMMADoc *pOperDoc = GetCurDoc()->GetParent();
		CGen2DOperation *pOper = dynamic_cast<CGen2DOperation*>(GetCurDoc()->GetParent());
		/*//Call MainFrame, Destroy all ChildFrames of Results
		if(pRes->GetChildViewDescriptor()){
		CWnd *pView = pRes->GetChildViewDescriptor()->GetView();
		if(pView){
		pView->GetTopLevelParent()->SendMessage(WMR_REOPEN_SAVES, (WPARAM)pRes);
		}
		}*/

		if (GetCurDoc()->GetSubDocNum() > 0) {
			pOperDoc->Remove();
		}
		else {

			size_t nNum = pRes->LoadResults();		//nNum for debug

		}
		/*C2DResults *pResults = new C2DResults;
		pResults->SetName(IDS_RESULTS);
		pOper->InsertSubDoc(pResults);
		m_wndTree.SelectItem(pResults->GetTreeItem());	//фокус на эл-нт

		int nNum = pResults->LoadResults();		//nNum for debug
		//*/
		FillTree();
	}



	/*
	C2DResults *pRes = dynamic_cast<C2DResults*>(GetCurDoc());
	CEMMADoc *pResDoc = GetCurDoc();
	CGen2Operation *pOper = dynamic_cast<CGen2Operation*>(GetCurDoc()->GetParent());

	if(pRes && pResDoc){

	//Call MainFrame, Destroy all ChildFrames of Results
	if(pRes->GetChildViewDescriptor()){
	CWnd *pView = pRes->GetChildViewDescriptor()->GetView();
	if(pView){
	pView->GetTopLevelParent()->SendMessage(WMR_REOPEN_SAVES, (WPARAM)pRes);
	}
	}

	//TODO: Разобраться с Result, C2DResult, CalcStage, C2DCalcStage
	// & Deactivate();

	C2DResults *pResults = new C2DResults;

	int nNum = pResults->LoadResults();		//nNum for debug
	if (nNum > 0){

	//pRes->Deactivate();	//если не удалять старое, то просто добавится новый раздел снизу
	//pRes->Remove();	// краш (выясняю)
	pResDoc->Remove();

	pResults->SetName(IDS_RESULTS);
	pOper->InsertSubDoc(pResults);
	m_wndTree.SelectItem(pResults->GetTreeItem());	//фокус на эл-нт

	}else{
	delete pResults;
	}

	FillTree();
	}
	*/
}

void CExplorerPane::OnUpdateResults(CCmdUI *pCmdUI) {
	C2DResults* pRes = dynamic_cast<C2DResults*>(GetCurDoc());
	if (pRes) {
		pCmdUI->Enable(true); //Включение или отключение элемента интерфейса для данной команды.
	}
	else {
		pCmdUI->Enable(false);
	}
}


void CExplorerPane::OnPressCloseButton() {
	CWnd* pFrame = GetParentFrame();

	//GetSafeHwnd() Returns the window handle for a window.
	//GetSafeHwnd() Returns NULL if the CWnd is not attached to a window or if it is used with a NULL CWnd pointer.
	if (pFrame->GetSafeHwnd() == nullptr) {
		CDlgShowError cError(ID_ERROR_EXPLORERPANE_GETSAFEWND_NULL); //_T("pFrame->GetSafeHwnd() is null"));
	}
	pFrame->SendMessage(WMR_EP_CLOSE);
}

void CExplorerPane::UpdatePane() {

}

//! Раскрываем все подэлементы в Обозревателе проекта
void CExplorerPane::ExpandAll() {

	HTREEITEM hItem = m_wndTree.GetFirstVisibleItem();

	while (hItem != nullptr)
	{
		m_wndTree.Expand(hItem, TVE_EXPAND);
		hItem = m_wndTree.GetNextItem(hItem, TVGN_NEXTVISIBLE);
	}
}
