#include "StdAfx.h"
#include "EMMADoc.h"

IOIMPL (CEMMADoc, CEMMADOC)

CEMMADoc::CEMMADoc(CString strName) 
	: m_hTreeItem(nullptr),
	  m_pParentDoc(nullptr),
	  m_pViewDescriptor(nullptr)
{
	RegisterMembers();
	if(strName)SetName(strName);
	m_icon = 0;	//по-умолчанию иконки нет
}

CEMMADoc::CEMMADoc(CRuntimeClass *pViewClass, size_t nID, CString strName)
	: m_hTreeItem(nullptr),
	  m_pParentDoc(nullptr),
	  m_pViewDescriptor(nullptr)
{
	RegisterMembers();

	if (pViewClass == nullptr) {
		CDlgShowError cError(ID_ERROR_EMMADOC_VIEWCLASS_NULL); //_T("pViewClass is null"));
		return;
	}

	InitViewDescriptor(pViewClass, nID);
	if(strName)SetName(strName);
	m_icon = 0;	//по-умолчанию иконки нет
}

CEMMADoc::~CEMMADoc(void){
	if(m_pViewDescriptor){
		CWnd *pView = m_pViewDescriptor->GetView();
		if(pView){
			//pView->GetTopLevelParent()->SendMessage(WMR_SUBDOCDESTROY, (WPARAM)this);
			pView->GetTopLevelFrame()->SendMessage(WMR_SUBDOCDESTROY, (WPARAM)this);
		}
		delete m_pViewDescriptor;
	}
	m_arSubDocs.DeleteObjects();
	//for(int i=0; i<GetSubDocNum(); i++){
	//	delete m_arSubDocs()[i];
	//}
}

void CEMMADoc::RegisterMembers(){
	RegisterMember(&m_strName);
	RegisterMember(&m_strFilePath);
	RegisterMember(&m_arSubDocs);
}

bool CEMMADoc::IsRemovable() {
	if (!m_pParentDoc) return false;
	return m_pParentDoc->IsSubDocRemovable(this);
}

bool CEMMADoc::LoadBody(CStorage& file){

	if (!IO::CClass::LoadBody(file)){
		return false;
	}

	size_t nSub = GetSubDocNum();

	for(size_t i=0; i<nSub; i++){
		CEMMADoc *pDoc = SubDoc(i);
		pDoc->SetParent(this);
	}

	return true;
}

void CEMMADoc::RemoveBranch(HTREEITEM hItem){
	ptrdiff_t nBranch = FindBranch(hItem);
	if (nBranch >= 0) {
		RemoveAt(nBranch);
	}
}

size_t CEMMADoc::FindBranch(HTREEITEM hItem){
	if(hItem == m_hTreeItem) return -2;

	for(size_t i=0; i<GetSubDocNum(); i++){
		if (SubDoc(i)->FindBranch(hItem) != -1) {
			return i;
		}
	}
	return -1;
}

CEMMADoc *CEMMADoc::FindDoc(HTREEITEM hItem){
	if(hItem == m_hTreeItem) return this;

	size_t nBranch = FindBranch(hItem);
	//if(nBranch >=0 && nBranch < GetSubDocNum()){
	if (nBranch < GetSubDocNum()) {
		return SubDoc(nBranch)->FindDoc(hItem);
	}

	return nullptr;
}

void CEMMADoc::RemoveAt(size_t nItem){

	//if (nItem < 0 || nItem >= GetSubDocNum()) {
	if (nItem >= GetSubDocNum()) {
		CDlgShowError cError(ID_ERROR_EMMADOC_NITEM_WRONG); //_T("EMMADoc: nItem is wrong"));
		return;
	}

	delete m_arSubDocs()[nItem];
	//m_arSubDocs().RemoveAt(nItem);
	m_arSubDocs().erase(m_arSubDocs().begin() + nItem);
}

size_t CEMMADoc::InsertSubDoc(CEMMADoc *pDoc, ptrdiff_t nIndex/* = -1*/ ){
	pDoc->SetParent(this);
	if (nIndex == -1) {
		m_arSubDocs().push_back(pDoc);
		return 0;	//нулевой индекс у первого элемента
	}

	m_arSubDocs().insert(m_arSubDocs().begin() + nIndex, pDoc);
	return nIndex;
}

int CEMMADoc::ChooseIconByID(int Icons_name, int Icon_name) {
	switch (Icons_name) {
	case IDS_TREE_ICONS_GEN2D:
		return ShowIconNumber2D(Icon_name);
		break;
	case IDS_TREE_ICONS_FORMING:
		return ShowIconNumberForming(Icon_name);
		break;
	case IDS_TREE_ICONS_ROLLING:
		return ShowIconNumberRolling(Icon_name);
		break;
	case IDS_TREE_ICONS_DRAGGING:
		return ShowIconNumberDragging(Icon_name);
		break;
	case IDS_TREE_ICONS_GEN3D:
		return ShowIconNumber3D(Icon_name);
		break;
	default:
		return 0;
		break;
	}
}

int CEMMADoc::ShowIconNumber2D( int IDS_IconName) {
	switch (IDS_IconName) {
	case IDS_ICON_GEN2D:
		return 0;
		break;
	case IDS_ICON_SPECIMEN:
		return 1;
		break;
	case IDS_ICON_SKETCH:
		return 2;
		break;
	case IDS_ICON_MESH:
		return 3;
		break;
	case IDS_ICON_MATERIALS:
		return 4;
		break;
	case IDS_ICON_OPERATION_2D:
		return 5;
		break;
	case IDS_ICON_TOOL:
		return 6;
		break;
	case IDS_ICON_TRAJECTORY:
		return 7;
		break;
	case IDS_ICON_BUILDER:
		return 8;
		break;
	case IDS_ICON_RESULT:
		return 9;
		break;
	case IDS_ICON_MATERIALDOC:
		return 10;
		break;
	case IDS_ICON_RESULTDOC:
		return 11;
		break;
	default:
		return 0;
		break;
	}
}

int CEMMADoc::ShowIconNumberForming( int IDS_IconName) {
	switch (IDS_IconName) {
	case IDS_ICON_FORMING:
		return 0;
		break;
	case IDS_ICON_SPECIMEN:
		return 1;
		break;
	case IDS_ICON_SKETCH:
		return 2;
		break;
	case IDS_ICON_MESH:
		return 3;
		break;
	case IDS_ICON_MATERIALS:
		return 4;
		break;
	case IDS_ICON_OPERATION_FORMING:
		return 5;
		break;
	case IDS_ICON_TOOL:
		return 6;
		break;
	case IDS_ICON_TRAJECTORY:
		return 7;
		break;
	case IDS_ICON_BUILDER:
		return 8;
		break;
	case IDS_ICON_RESULT:
		return 9;
		break;
	case IDS_ICON_MATERIALDOC:
		return 10;
		break;
	case IDS_ICON_RESULTDOC:
		return 11;
		break;
	default:
		return 0;
		break;
	}
}

int CEMMADoc::ShowIconNumberRolling( int IDS_IconName) {
	switch (IDS_IconName) {
	case IDS_ICON_ROLLING:
		return 0;
		break;
	case IDS_ICON_PODKAT:
		return 1;
		break;
	case IDS_ICON_SKETCH:
		return 2;
		break;
	case IDS_ICON_MESH:
		return 3;
		break;
	case IDS_ICON_MATERIALS:
		return 4;
		break;
	case IDS_ICON_ROLLINGOPERATION:
		return 5;
		break;
	case IDS_ICON_CALIBERS:
		return 6;
		break;
	case IDS_ICON_TOPROLLER:
		return 7;
		break;
	case IDS_ICON_BOTTOMROLLER:
		return 8;
		break;
	case IDS_ICON_BUILDER:
		return 9;
		break;
	case IDS_ICON_RESULT:
		return 10;
		break;
	case IDS_ICON_MATERIALDOC:
		return 11;
		break;
	case IDS_ICON_RESULTDOC:
		return 12;
		break;
	default:
		return 0;
		break;
	}
}

int CEMMADoc::ShowIconNumberDragging( int IDS_IconName) {
	switch (IDS_IconName) {
	case IDS_ICON_DRAGGING:
		return 0;
		break;
	case IDS_ICON_SPECIMEN:
		return 1;
		break;
	case IDS_ICON_SKETCH:
		return 2;
		break;
	case IDS_ICON_MESH:
		return 3;
		break;
	case IDS_ICON_MATERIALS:
		return 4;
		break;
	case IDS_ICON_OPERATION_DRAGGING:
		return 5;
		break;
	case IDS_ICON_TOOL:
		return 6;
		break;
	case IDS_ICON_TRAJECTORY:
		return 7;
		break;
	case IDS_ICON_BUILDER:
		return 8;
		break;
	case IDS_ICON_RESULT:
		return 9;
		break;
	case IDS_ICON_MATERIALDOC:
		return 10;
		break;
	case IDS_ICON_RESULTDOC:
		return 11;
		break;
	default:
		return 0;
		break;
	}
}

int CEMMADoc::ShowIconNumber3D( int IDS_IconName) {
	switch (IDS_IconName) {
	case IDS_ICON_GEN3D:
		return 0;
		break;
	default:
		return 0;
		break;
	}
}


void CEMMADoc::SetIconName(int name) {
	m_icon = name;
}

int CEMMADoc::GetIconName() {
	return m_icon;
}

HTREEITEM CEMMADoc::AddToTree(CTreeCtrl &TreeCtrl, HTREEITEM hParent, HTREEITEM hInsertAfter, int Icons_name){
	
	int nIconNumber = ChooseIconByID(Icons_name, GetIconName());

	m_hTreeItem = TreeCtrl.InsertItem(m_strName(), nIconNumber, nIconNumber, hParent, hInsertAfter);
	int res = TreeCtrl.SetItemData(m_hTreeItem, (size_t)this);	//для x32 (DWORD u_long) == (Size_t u_int), для x64 будет __int64 у обоих

	size_t nSubDocNum = GetSubDocNum();
	for(size_t i=0; i<nSubDocNum; i++){
		SubDoc(i)->AddToTree(TreeCtrl, m_hTreeItem, TVI_LAST, Icons_name);
	}

	return m_hTreeItem;
}

HTREEITEM CEMMADoc::ModifyTree(CTreeCtrl &TreeCtrl, int Icons_name){
	
	HTREEITEM hItem = TreeCtrl.GetNextItem(m_hTreeItem, TVGN_CHILD);
	ptrdiff_t nStart = 0;

	while(hItem != nullptr){
		HTREEITEM hNext = TreeCtrl.GetNextSiblingItem(hItem);
		HTREEITEM hPrev = TreeCtrl.GetPrevSiblingItem(hItem);
		
		if (hPrev == nullptr) {
			hPrev = TVI_FIRST;
		}
		
		ptrdiff_t n = -1; // Номер элемента массива поддокументов, соответствующего hItem
		
		for(size_t i=nStart; i<GetSubDocNum(); i++){
			if(SubDoc(i)->GetTreeItem() == hItem){
				n=i; 
				break;
			}
		}
		if(n==-1){// Если такого элемента не найдено, удаляем ветку дерева
			int resDel = TreeCtrl.DeleteItem(hItem);
		}else{
			// Добавляем в дерево ветки, соответствующие поддокументам, которые лежат перед n
			for(ptrdiff_t i=n-1; i>=nStart; i--){
				SubDoc(i)->AddToTree(TreeCtrl, m_hTreeItem, hPrev, Icons_name);
			}
			// Модифицируем ветку hItem, соответствующую n-ному элементу массива
			SubDoc(n)->ModifyTree(TreeCtrl, Icons_name);
			nStart = n+1;
		}
		hItem = hNext;
	}

	// Добавляем в дерево оставшиеся элементы массива
	size_t nSubDocNum = GetSubDocNum();
	for(size_t i=nStart; i<nSubDocNum; i++){
		SubDoc(i)->AddToTree(TreeCtrl, m_hTreeItem, TVI_LAST, Icons_name);
	}

	int nIconNumber = ChooseIconByID(Icons_name, GetIconName());
	TreeCtrl.SetItemImage(m_hTreeItem, nIconNumber, nIconNumber);
	int resSetText = TreeCtrl.SetItemText(m_hTreeItem, m_strName());	//не 0, если успешно

	return m_hTreeItem;
}

void CEMMADoc::InitialUpdateView() {
	if (GetView()) {
		GetView()->SendMessage(WMR_INITIALUPDATE);
	}
}

void CEMMADoc::InitialUpdatePane() {
	if (m_pViewDescriptor != nullptr) {

		CWnd *pView = m_pViewDescriptor->GetView();
		if (pView != nullptr) {
			pView->GetTopLevelParent()->SendMessage(WMR_UPDATEPANE, (WPARAM)this);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////
//CEMMADoc CViewDescriptor functions
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
void CEMMADoc::InitViewDescriptor(CRuntimeClass *pViewClass, size_t nID){
	m_pViewDescriptor = new CViewDescriptor(pViewClass, nID);
}

void CEMMADoc::SetView(CWnd *pView){
	if (GetViewDescriptor()) {
		GetViewDescriptor()->SetView(pView);
	}
	InitialUpdateView();
}

CWnd* CEMMADoc::GetView(){
	CViewDescriptor *tViewDescriptor = GetViewDescriptor();
	if (!tViewDescriptor) {
		return nullptr;
	}
	CWnd* pView = tViewDescriptor->GetView();
	return pView;
}

CRuntimeClass *CEMMADoc::GetViewClass(){
	if(!GetViewDescriptor()) return nullptr;
	return GetViewDescriptor()->GetViewClass();
}

HMENU CEMMADoc::GetHMenu(){
	if(!GetViewDescriptor()) return nullptr;
	return GetViewDescriptor()->GetHMenu();
}

HACCEL CEMMADoc::GetHAccel(){
	if(!GetViewDescriptor()) return nullptr;
	return GetViewDescriptor()->GetHAccel();
}

UNINT CEMMADoc::GetResID(){
	if(!GetViewDescriptor()) return 0;
	return GetViewDescriptor()->GetResID();
}


CEMMADoc *CEMMADoc::GetRoot(){
	CEMMADoc *pRoot = this;
	while (pRoot->m_pParentDoc) {
		pRoot = pRoot->m_pParentDoc;
	}
	return pRoot;
}

CEMMADoc *CEMMADoc::GetParent(){
	CEMMADoc *pRoot = this;
	pRoot = pRoot->m_pParentDoc;
	return pRoot;
}