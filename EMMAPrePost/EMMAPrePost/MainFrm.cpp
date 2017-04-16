// MainFrm.cpp : implementation of the CMainFrame class
//

#include "stdafx.h"
#include "MainFrm.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CMainFrame


IMPLEMENT_DYNAMIC(CMainFrame, CMDIFrameWndEx)

BEGIN_MESSAGE_MAP(CMainFrame, CMDIFrameWndEx)
	ON_WM_CREATE()
	// Global help commands
	ON_COMMAND(ID_HELP, &CMDIFrameWndEx::OnHelp)
	ON_COMMAND(ID_HELP_FINDER, &CMDIFrameWndEx::OnHelpFinder)
	ON_COMMAND(ID_CONTEXT_HELP, &CMDIFrameWndEx::OnContextHelp)
	ON_COMMAND(ID_DEFAULT_HELP, &CMDIFrameWndEx::OnHelpFinder)
	ON_COMMAND_RANGE(ID_VIEW_APPLOOK_WIN_2000, ID_VIEW_APPLOOK_WINDOWS_7, &CMainFrame::OnApplicationLook)
	ON_UPDATE_COMMAND_UI_RANGE(ID_VIEW_APPLOOK_WIN_2000, ID_VIEW_APPLOOK_WINDOWS_7, &CMainFrame::OnUpdateApplicationLook)
	// Новый проект\задача
	ON_COMMAND(ID_FILE_NEW, &CMainFrame::OnFileNew)
	// Сохранить
	ON_COMMAND(ID_FILE_SAVE, &CMainFrame::OnFileSave)
	ON_UPDATE_COMMAND_UI(ID_FILE_SAVE, &CMainFrame::OnUpdateFileSave)
	
	// Сохранить
	ON_COMMAND(ID_FILE_SAVE_AS, &CMainFrame::OnFileSaveAs)
	ON_UPDATE_COMMAND_UI(ID_FILE_SAVE_AS, &CMainFrame::OnUpdateFileSaveAs)
	
	// Открыть
	ON_COMMAND(ID_FILE_OPEN, &CMainFrame::OnFileOpen)
	// Закрыть проект
	ON_COMMAND(ID_FILE_CLOSE, &CMainFrame::OnFileClose)
	ON_UPDATE_COMMAND_UI(ID_FILE_CLOSE, &CMainFrame::OnUpdateFileClose)
	// Сменился активный дочерний фрейм
	ON_REGISTERED_MESSAGE(WMR_CHILDACTIVATE, OnChildFrm)
	// Закрывается дочернее окно
	ON_REGISTERED_MESSAGE(WMR_CHILDCLOSE, OnChildClose)
	// Сменился активный элемент дерева в ExplorerPane
	ON_REGISTERED_MESSAGE(WMR_EP_SELCHANGE, OnEPSelChange)
	// Пользователь хочет закрыть ExplorerPane
	ON_REGISTERED_MESSAGE(WMR_EP_CLOSE, OnEPClose)
	// Уничтожен элемент документа
	ON_REGISTERED_MESSAGE(WMR_SUBDOCDESTROY, OnSubDocDestroy)
	/////Roma
	//Обновление всех видимых панелей
	ON_REGISTERED_MESSAGE(WMR_UPDATEPANE, OnUpdatePane)
	//Новая загрузка файлов с результатами
	ON_REGISTERED_MESSAGE(WMR_REOPEN_SAVES, OnReopenSaves)

	//RIBBON ExplorerPane messages
	ON_COMMAND(ID_RIBBON_EXPLORER_ELEMENT_ADD, &CMainFrame::OnExplorerElementAdd)
	ON_UPDATE_COMMAND_UI(ID_RIBBON_EXPLORER_ELEMENT_ADD, &CMainFrame::OnUpdateExplorerElementAdd)
	ON_COMMAND(ID_RIBBON_EXPLORER_ELEMENT_DELETE, &CMainFrame::OnExplorerElementDelete)
	ON_UPDATE_COMMAND_UI(ID_RIBBON_EXPLORER_ELEMENT_DELETE, &CMainFrame::OnUpdateExplorerElementDelete)

	//Cut, Copy, Paste
	ON_COMMAND(ID_EDIT_CUT, &CMainFrame::OnEditCut)
	ON_UPDATE_COMMAND_UI(ID_EDIT_CUT, &CMainFrame::OnUpdateNeedSel)
	ON_COMMAND(ID_EDIT_COPY, &CMainFrame::OnEditCopy)
	ON_UPDATE_COMMAND_UI(ID_EDIT_COPY, &CMainFrame::OnUpdateNeedSel)
	ON_COMMAND(ID_EDIT_PASTE, &CMainFrame::OnEditPaste)
	ON_UPDATE_COMMAND_UI(ID_EDIT_PASTE, &CMainFrame::OnUpdateNeedClip)

	//Sketch
	ON_COMMAND(ID_NEW_NODE, &CMainFrame::OnNewNode)
	ON_UPDATE_COMMAND_UI(ID_NEW_NODE, &CMainFrame::OnUpdateNewNode)

	ON_COMMAND(ID_NEW_CURVE, &CMainFrame::OnNewCurve)
	ON_UPDATE_COMMAND_UI(ID_NEW_CURVE, &CMainFrame::OnUpdateNewCurve)

	ON_COMMAND(ID_NEW_CIRCLE_ARC, &CMainFrame::OnNewCircleArc)
	ON_UPDATE_COMMAND_UI(ID_NEW_CIRCLE_ARC, &CMainFrame::OnUpdateNewCircleArc)

	ON_COMMAND(ID_FACET, &CMainFrame::OnNewFacet)
	ON_UPDATE_COMMAND_UI(ID_FACET, &CMainFrame::OnUpdateNewFacet)

	ON_COMMAND(ID_CLEAR_SKETCH, &CMainFrame::OnClearSketch)
	ON_UPDATE_COMMAND_UI(ID_CLEAR_SKETCH, &CMainFrame::OnUpdateClearSketch)

	ON_COMMAND(ID_2DSKETCHPANE_BUTTON_SHOWCONTOURS, CMainFrame::OnSetContours)
	ON_UPDATE_COMMAND_UI(ID_2DSKETCHPANE_BUTTON_SHOWCONTOURS, CMainFrame::OnUpdateSketchToolTrajectoryDoc)

	ON_COMMAND(ID_2DSKETCHPANE_BUTTON_SAVE_TO_META, &CMainFrame::OnSaveSketchToMeta)
	ON_UPDATE_COMMAND_UI(ID_2DSKETCHPANE_BUTTON_SAVE_TO_META, &CMainFrame::OnUpdateSketchToolTrajectoryDoc)

	//Tool
	ON_COMMAND(ID_NEW_TOOL_NODE, &CMainFrame::OnNewToolNode)
	ON_UPDATE_COMMAND_UI(ID_NEW_TOOL_NODE, &CMainFrame::OnUpdateNewToolNode)

	//Trajectory
	ON_COMMAND(ID_PLAY, &CMainFrame::OnPlay)
	ON_UPDATE_COMMAND_UI(ID_PLAY, &CMainFrame::OnUpdatePlay)

	//Mesher
	
	ON_COMMAND(ID_CLEAR_MESH, &CMainFrame::OnClearMesh)
	ON_UPDATE_COMMAND_UI(ID_CLEAR_MESH, &CMainFrame::OnUpdateMesherDoc)
	ON_COMMAND_RANGE(ID_MESH_TYPE_FRONT, ID_MESH_TYPE_MSW, &CMainFrame::OnSetMeshType)
	ON_UPDATE_COMMAND_UI_RANGE(ID_MESH_TYPE_FRONT, ID_MESH_TYPE_MSW, &CMainFrame::OnUpdateSetMeshType)
	ON_UPDATE_COMMAND_UI(ID_MESH_TYPE_BUTTON, &CMainFrame::OnUpdateMesherDoc)

	ON_COMMAND(ID_NEW_MESH, &CMainFrame::OnNewMesher)
	ON_UPDATE_COMMAND_UI(ID_NEW_MESH, &CMainFrame::OnUpdateMesherDoc)
	ON_COMMAND(ID_IMPORT_SKETCH, &CMainFrame::OnImportSketch)
	ON_UPDATE_COMMAND_UI(ID_IMPORT_SKETCH, &CMainFrame::OnUpdateMesherDoc)

	ON_COMMAND(ID_NEW_MESHPOINT, &CMainFrame::OnNewMeshPoint)
	ON_UPDATE_COMMAND_UI(ID_NEW_MESHPOINT, &CMainFrame::OnUpdateMesherDoc)
	//ON_BN_CLICKED(ID_UPDATE_CONTOUR, &CMainFrame::OnShowContours)
	//ON_UPDATE_COMMAND_UI(ID_UPDATE_CONTOUR, &CMainFrame::OnUpdateMesherDoc)

	//Builder
	ON_COMMAND(ID_BUILD_SAVE, &CMainFrame::OnBuilderSave)
	ON_UPDATE_COMMAND_UI(ID_BUILD_SAVE, &CMainFrame::OnUpdateBuilderDoc)

	//Results
	ON_COMMAND(ID_EXPLORER_RESULTS, &CMainFrame::OnLoadResults)
	ON_UPDATE_COMMAND_UI(ID_EXPLORER_RESULTS, &CMainFrame::OnUpdateResultsDoc)

	ON_COMMAND(ID_EXPLORER_RESULTS_CLEAR, &CMainFrame::OnClearResults)
	ON_UPDATE_COMMAND_UI(ID_EXPLORER_RESULTS_CLEAR, &CMainFrame::OnUpdateResultsDoc)
		
	ON_WM_MDIACTIVATE()
	ON_WM_CLOSE()
	
	//ON_WM_ERASEBKGND()
	//ON_WM_PAINT()

END_MESSAGE_MAP()

static UNINT g_indicators[] =
{
	ID_SEPARATOR,      // status line indicator
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

//Массив ToolBar
const UNINT g_ToolBarIDs[] =
{
	//IDR_MAIN_TB,
	//IDR_2DVIEW_TB,
	//IDR_3DVIEW_TB,
	// ---- ROMA
	IDR_2DTRAJECTORY_TB,
	// YURA
	IDR_2DSKETCH_TB,
	IDR_2DMESH_TB,
	IDR_2DTOOL_TB,
	IDR_2DBUILDER_TB,

	// TIMUR
	IDR_2DCALCSTAGE_TB,
	// ----
	IDR_3DCALCSTAGE_TB
};

//////////////////////////////////////////////////////////////////////////
// CMainFrame construction/destruction
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
CMainFrame::CMainFrame() :m_wndExplorerPane(IDR_EXPLORER_VIEW) {
	CString str_AppLook;
	int n = str_AppLook.LoadString(IDS_APP_LOOK);
	theApp.m_nAppLook = theApp.GetInt(str_AppLook, ID_VIEW_APPLOOK_VS_2005);
	m_pDoc = nullptr;

	//m_DockBars.push_back(new CPropertiesPane(IDR_PROPERTIES_PANE));
	m_DockBars.push_back(new C2DSketchPane(IDR_2DSKETCH_PANE));
	m_DockBars.push_back(new C2DMeshPane(IDR_2DMESH_PANE));
	m_DockBars.push_back(new CMaterialPane(IDR_MATERIALS_PANE));
	m_DockBars.push_back(new C2DTrajectoryPane(IDR_2DTRAJECTORY_PANE));
	m_DockBars.push_back(new C2DToolPane(IDR_2DTOOL_PANE));
	m_DockBars.push_back(new C2DBuilderPane(IDR_2DBUILDER_PANE));
	m_DockBars.push_back(new C2DCalcStagePane(IDR_2DCALCSTAGE_PANE));
	m_DockBars.push_back(new C3DCalcStagePane(IDR_3DCALCSTAGE_PANE));
}

CMainFrame::~CMainFrame() {
	for (size_t i = 0; i<m_ToolBars.size(); i++) {
		//if(m_ToolBars[i])
		delete m_ToolBars[i];
	}
	m_ToolBars.resize(0);

	for (size_t i = 0; i<m_DockBars.size(); i++) {
		delete m_DockBars[i];
	}/*/
	 m_DockBars.resize(0); //*/
	DestroyDocument();
}


//////////////////////////////////////////////////////////////////////////
// CMainFrame implementation
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
int CMainFrame::LoadFrame(UNINT nIDResource, UNLONG dwDefaultStyle, CWnd* pParentWnd, CCreateContext* pContext) {
	if (!CMDIFrameWndEx::LoadFrame(nIDResource, dwDefaultStyle, pParentWnd, pContext)) {
		return 0;
	}

	int nSize = sizeof(g_ToolBarIDs) / sizeof(g_ToolBarIDs[0]);
	m_ToolBars.resize(nSize);

	for (int i = 0; i<nSize; i++) {
		m_ToolBars[i] = new CMFCToolBar;
		if (!m_ToolBars[i]->CreateEx(this, TBSTYLE_FLAT, WS_CHILD | CBRS_TOP | CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC, CRect(1, 1, 1, 1), 1) || !m_ToolBars[i]->LoadToolBar(g_ToolBarIDs[i])) {
			CDlgShowError cError(ID_ERROR_MAINFRM_TOOLBAR_FAIL); //_T("Failed to create toolbar \n"));
			return -1;   // fail to create
		}
		CString strToolBarName;
		int bNameValid = strToolBarName.LoadString(g_ToolBarIDs[i]);

		if (!bNameValid) {
			CDlgShowError cError(ID_ERROR_MAINFRM_WRONG_TOOLBAR_NAME); //_T("Wrong Toolbar name"));
			return -1;
		}		

		m_ToolBars[i]->SetWindowText(strToolBarName);//*/
	}

	//m_wndMenuBar.EnableDocking(CBRS_ALIGN_ANY);
	for (int i = 0; i < nSize; i++) {
		m_ToolBars[i]->EnableDocking(CBRS_ALIGN_ANY);
	}
	EnableDocking(CBRS_ALIGN_ANY);
	//DockPane(&m_wndMenuBar);
	for (int i = 0; i < nSize; i++) {
		DockPane(m_ToolBars[i]);
	}

	CreateDockPanes();
	ShowControlBars();

	return 1;
}

/*
/// Ribbon commands
*/

//determine if it is a CEdit control that has the focus :
int CMainFrame::IsEdit(CWnd* pWnd) {
	if (pWnd == nullptr) {
		return 0;
	};

	HWND hWnd = pWnd->GetSafeHwnd();
	if (hWnd == 0) {
		return 0;
	}

	// Use the SDK ::GetClassName() function because MFC IsKindOf()
	// fails if no CEdit variable has been created for the control you're trying to test.

	TCHAR szClassName[6];
	int n1 = ::GetClassName(hWnd, szClassName, 6);
	int n2 = (_tcsicmp(szClassName, _T("Edit")) == 0);
	return  n1 && n2;
}

// UPDATE_COMMAND_UI handler for Edit Copy and Edit Cut which both
// require that the current focus is on an edit control that has text selected.

void CMainFrame::OnUpdateNeedSel(CCmdUI* pCmdUI)
{
	// get the current focus & determine if its on a CEdit control
	CWnd* pWnd = GetFocus();
	if (nullptr == pWnd || !IsEdit(pWnd)){
		pCmdUI->Enable(0);
	}else{
		CEdit* pEdit = (CEdit*)pWnd;
		int nBeg, nEnd;

		pEdit->GetSel(nBeg, nEnd);
		pCmdUI->Enable(nBeg != nEnd);
	}
}

// UPDATE_COMMAND_UI handlers for Edit Paste requires that focus be
// on an edit control and that the clipboard contains text to be
// pasted into the control.

void CMainFrame::OnUpdateNeedClip(CCmdUI* pCmdUI)
{
	// get the current focus & determine if its on a CEdit control
	// also check to see if the control is read-only.
	CWnd* pWnd = GetFocus();
	if (nullptr == pWnd ||	!IsEdit(pWnd) || (pWnd->GetStyle() & ES_READONLY) != 0){
		pCmdUI->Enable(0);
	}else{
		pCmdUI->Enable(::IsClipboardFormatAvailable(CF_TEXT));
	}
}


void CMainFrame::OnEditCopy()
{
	CEdit* pEdit = (CEdit*)GetFocus();
	ASSERT(IsEdit(pEdit));
	pEdit->Copy();
}

void CMainFrame::OnEditCut()
{
	CEdit* pEdit = (CEdit*)GetFocus();
	ASSERT(IsEdit(pEdit));
	pEdit->Cut();
}

void CMainFrame::OnEditPaste()
{
	CEdit* pEdit = (CEdit*)GetFocus();
	ASSERT(IsEdit(pEdit));
	ASSERT(::IsClipboardFormatAvailable(CF_TEXT));
	pEdit->Paste();
}

// ExplorerPane
void CMainFrame::OnExplorerElementAdd() {
	// Добавление нового подэлемента (Операции, Инструмента и прочее)
	m_wndExplorerPane.OnAdd();
}

void CMainFrame::OnUpdateExplorerElementAdd(CCmdUI* pCmdUI) {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCmdUI->Enable(pCurDoc->CanAdd());
	}
	else {
		pCmdUI->Enable(false);
	}
}

void CMainFrame::OnExplorerElementDelete() {
	// Добавление нового подэлемента (Операции, Инструмента и прочее)
	m_wndExplorerPane.OnRemove();
}

void CMainFrame::OnUpdateExplorerElementDelete(CCmdUI* pCmdUI) {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCmdUI->Enable(pCurDoc->IsRemovable());
	}
	else {
		pCmdUI->Enable(false);
	}
}

// Sketch

// Метод для включения всех кнопок Чертежа
void CMainFrame::OnUpdateSketchToolTrajectoryDoc(CCmdUI* pCmdUI) {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		bool bType1 = (pCurDoc->GetType() == C2DSKETCH);
		bool bType2 = (pCurDoc->GetType() == C2DTOOL);
		bool bType3 = (pCurDoc->GetType() == C2DTRAJECTORY);
		pCmdUI->Enable(bType1 || bType2 || bType3);
	}
	else {
		pCmdUI->Enable(false);
	}
}

void CMainFrame::OnUpdateSetCursorCheck(CCmdUI* pCmdUI, int nType){
	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		CWnd *pView = GetCurDoc()->GetView();
		if (pView) {
			CGLBasicView *pBasicView = dynamic_cast<CGLBasicView*>(GetCurDoc()->GetView());
			if ((pBasicView != nullptr) && pBasicView->GetCursorType() == nType) {
				pCmdUI->SetCheck(1);
			}
			else {
				pCmdUI->SetCheck(0);
			}
		}
	}
}

void CMainFrame::OnNewNode() {
	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCurDoc->GetView()->SendMessage(WMR_SETCURSOR, 100);
	}
	
}

void CMainFrame::OnUpdateNewNode(CCmdUI *pCmdUI) {
	OnUpdateSketchToolTrajectoryDoc(pCmdUI);
	OnUpdateSetCursorCheck(pCmdUI, 100);
}

void CMainFrame::OnNewCurve() {
	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCurDoc->GetView()->SendMessage(WMR_SETCURSOR, 101);
	}

}

void CMainFrame::OnUpdateNewCurve(CCmdUI *pCmdUI) {
	OnUpdateSketchToolTrajectoryDoc(pCmdUI);
	OnUpdateSetCursorCheck(pCmdUI, 101);
}

void CMainFrame::OnNewCircleArc() {
	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCurDoc->GetView()->SendMessage(WMR_SETCURSOR, 102);
	}

}

void CMainFrame::OnUpdateNewCircleArc(CCmdUI *pCmdUI) {
	OnUpdateSketchToolTrajectoryDoc(pCmdUI);
	OnUpdateSetCursorCheck(pCmdUI, 102);
}

void CMainFrame::OnNewFacet() {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCurDoc->GetView()->SendMessage(WMR_SETCURSOR, 103);
	}
}

void CMainFrame::OnUpdateNewFacet(CCmdUI *pCmdUI) {

	OnUpdateSketchToolTrajectoryDoc(pCmdUI);
	OnUpdateSetCursorCheck(pCmdUI, 103);
}

void CMainFrame::OnClearSketch() {

	I2DSketch* pSketch = GetCurDoc() ? dynamic_cast<I2DSketch*> (GetCurDoc()) : nullptr;
	if (pSketch) {
		pSketch->SetClearOutline();
	}
}

void CMainFrame::OnUpdateClearSketch(CCmdUI *pCmdUI) {
	OnUpdateSketchToolTrajectoryDoc(pCmdUI);
}

//Сохранение всего чертежа в meta-файл
void CMainFrame::OnSaveSketchToMeta() {
	
	I2DSketch* pSketch = GetCurDoc() ? dynamic_cast<I2DSketch*> (GetCurDoc()) : nullptr;
	if (pSketch) {
		C2DOutline *pOutline = pSketch->GetOutline();
		//Две точки могут образовать контур (две дуги из двух точек с разными радиусами), меньше - нет
		if (pOutline->GetNodeCount() > 1) {

			//Формирование строки для записи в файл
			CString sMetaToSave = pOutline->SaveMetaToString();

			//Показываем окно для сохранения .meta файла
			CFileDialog fileDlg(0, _T("meta"), 0, 4 | 2, _T("meta file (*.meta)|*.meta|"));
			if (fileDlg.DoModal() == IDOK) {
				CString fname = fileDlg.GetPathName();
				CString fpath = fname;

				//обрабатываем путь
				int nPosSlash = fpath.ReverseFind('\\');
				fname = fname.Right(fname.GetLength() - nPosSlash - 1);
				fpath.Truncate(nPosSlash + 1);

				CFileException e;
				CStorage file;
				if (file.Open(fpath + fname, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeText | CStorage::shareDenyWrite, &e)) {
					file.WriteString(sMetaToSave);
				}
				else {
					CString str_error_fnf;
					int nload = str_error_fnf.LoadStringW(ID_ERROR_2DSKETCHPANE_NO_FILE);
					CDlgShowError cError(str_error_fnf + AllToString(e.m_cause));	//Показываем окно ошибки
				}
				file.Close();
			}
		}
		else {
			CDlgShowError Diag_err(ID_ERROR_2DSKETCHPANE_NO_NODES);	//Показываем окно ошибки
		}
	}
	else {
		CDlgShowError Diag_err(ID_ERROR_2DSKETCHPANE_NO_NODES);	//Показываем окно ошибки
	}
}


void CMainFrame::OnSetContours() {
	I2DSketch* pI2DSketch = GetCurDoc() ? dynamic_cast<I2DSketch*> (GetCurDoc()) : nullptr;
	if (pI2DSketch) {
		pI2DSketch->ClearActiveObjects();
		pI2DSketch->SetContours();
	}
	else {
		CDlgShowError Diag_err(ID_ERROR_2DSPECIMEN_SKETCH_NULL);	//Показываем окно ошибки
	}
}

/////////////////
// Mesh

void CMainFrame::OnClearMesh() {
	
	CEMMADoc *pCurDoc = GetCurDoc();
	C2DMesher *pMesher = pCurDoc ? dynamic_cast<C2DMesher*>(pCurDoc) : nullptr;
	if (pMesher) {
		CMesh* pMesh = pMesher->GetMesh();
		if (pMesh) {
			pMesh->DeleteMesh();
			pMesh->DeleteMeshPoint();
		}
		pMesher->GetView()->Invalidate();

	}
}

void CMainFrame::OnSetMeshType(UNINT id) {

	CEMMADoc *pCurDoc = GetCurDoc();
	C2DMesher *pMesher = pCurDoc ? dynamic_cast<C2DMesher*>(pCurDoc) : nullptr;
	if (pMesher) {

		pMesher->SetMeshType(id);
	}
}

void CMainFrame::OnUpdateSetMeshTypeMainButton(CCmdUI* pCmdUI) {
	CEMMADoc *pCurDoc = GetCurDoc();
	C2DMesher *pMesher = pCurDoc ? dynamic_cast<C2DMesher*>(pCurDoc) : nullptr;
	if (pMesher) {
		pCmdUI->Enable(true);
	}
	else {
		pCmdUI->Enable(false);
	}
}


void CMainFrame::OnUpdateSetMeshType(CCmdUI* pCmdUI) {

	bool bAppLook = false;
	CEMMADoc *pCurDoc = GetCurDoc();
	C2DMesher *pMesher = pCurDoc ? dynamic_cast<C2DMesher*>(pCurDoc) : nullptr;
	if (pMesher) {
	UNINT nType = pMesher->GetMeshType();
		bAppLook = (nType == pCmdUI->m_nID);
	}
	
	pCmdUI->SetRadio(bAppLook);
}

void CMainFrame::OnNewMesher() {

	CEMMADoc *pCurDoc = GetCurDoc();
	C2DMesher* pMesher = pCurDoc ? dynamic_cast<C2DMesher*>(pCurDoc) : nullptr;
	//TODO: продумать проверки
	if (pMesher) {
		if (pMesher->GetOutline()->GetCurveCount()>0) {
			CMesh* pMesh = pMesher->GetMesh();
			if (pMesh != nullptr) {
				pMesh->DeleteMesh();
			}

			//TODO: проверка на чертёж (возможно не нужна, если Женин вариант будет работать)
			pMesh->GenerateMesh(pMesh->m_dRibLen, pMesher->GetOutline());
			pMesh->ConstructMap();

		}
		else {
			//Нету чертежа или не загружен
			CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_CONTOURS);
		}
		pCurDoc->GetView()->Invalidate();
	}
	else {
		//Нету документа (Mesher-a) или не загружен
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_SKETCH);
	}
	
}

/*
void CMainFrame::OnShowContours() {

	CEMMADoc *pCurDoc = GetCurDoc();
	C2DMesher* pMesher = pCurDoc ? dynamic_cast<C2DMesher*>(pCurDoc) : nullptr;
	C2DOutline *pOutline = pMesher->GetOutline();

	if (!pOutline) {
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_CONTOURS); //_T("Contours not found."));	//Показываем окно ошибки
		return;
	}

	if (pOutline->FindContours(pMesher->GetMesh()->m_dRibLen))
		pMesher->GetView()->SendMessage(WMR_INITIALUPDATE); // обновляем представление
	else
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_CONTOURS); //_T("Contours not found."));	//Показываем окно ошибки
}
*/

void CMainFrame::OnNewMeshPoint() {

	if (GetCurDoc()) {
		GetCurDoc()->GetView()->SendMessage(WMR_SETCURSOR, 300);
	}
	
	//UpdatePane();
}

void CMainFrame::OnImportSketch() {

	CEMMADoc *pCurDoc = GetCurDoc();
	C2DMesher* pMesher = pCurDoc ? dynamic_cast<C2DMesher*>(pCurDoc) : nullptr;

	C2DSketch* parent = dynamic_cast<C2DSketch*>(pCurDoc->GetParent()->SubDoc(0));
	C2DOutline *pOutline = nullptr;
	if (parent) {
		pOutline = parent->GetOutline();
	}
	

	if (pOutline) {
		//TODO: проверка должна быть полнее - например, на незамкнутость контура(ов)
		if (pOutline->GetNodeCount() < 3 || pOutline->GetCurveCount() < 2) {
			pOutline = nullptr;

			//вывод сообщения об ошибке
			CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_SKETCH);	//Показываем окно ошибки
			return;
		}
		pMesher->SetOutline(pOutline);
	}
	else {
		//вывод сообщения об ошибке
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_SKETCH);	//Показываем окно ошибки
		return;
	}
	//*/
}

// Метод для включения всех кнопок КЭ сетки 
void CMainFrame::OnUpdateMesherDoc(CCmdUI* pCmdUI) {

	RibbonButtonOnOff(pCmdUI, C2DMESHER);
}

/////////////////
//Tool

void CMainFrame::OnNewToolNode() {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		pCurDoc->GetView()->SendMessage(WMR_SETCURSOR, 200);
	}
}

void CMainFrame::OnUpdateNewToolNode(CCmdUI *pCmdUI) {

	RibbonButtonOnOff(pCmdUI, C2DTOOL);
	OnUpdateSetCursorCheck(pCmdUI, 200);

}

//////////////////
//Trajectory


void CMainFrame::OnPlay() {

	C2DTrajectory* pTraj = GetCurDoc() ? dynamic_cast<C2DTrajectory*>(GetCurDoc()) : nullptr;
	if (!pTraj) return;

	pTraj->m_dCurTime = 0.0;

	if (!pTraj->m_fPlay) {
		pTraj->m_fPlay = true;
		pTraj->GetView()->SetTimer(51, 4, nullptr);
	}
	else {
		pTraj->m_fPlay = false;
		pTraj->GetView()->KillTimer(51);
	}
}

void CMainFrame::OnUpdatePlay(CCmdUI *pCmdUI) {

	RibbonButtonOnOff(pCmdUI, C2DTRAJECTORY);
	C2DTrajectory* pTraj = GetCurDoc() ? dynamic_cast<C2DTrajectory*>(GetCurDoc()) : nullptr;
	if (!pTraj) return;
	if (pTraj->m_fPlay) {
		pCmdUI->SetCheck(1);
	}
	else {
		pCmdUI->SetCheck(0);
	}
	//*/
}

//Builder

void CMainFrame::OnBuilderSave() {

	CEMMADoc *pCurDoc = GetCurDoc();
	C2DBuilder *pBuilder = dynamic_cast<C2DBuilder *>(pCurDoc);	//получаем ссылку на документ (Сборка)
	if (pBuilder) {
		pBuilder->Build();
	}
}

//Активируем кнопки Сборки
void CMainFrame::OnUpdateBuilderDoc(CCmdUI *pCmdUI) {
	RibbonButtonOnOff(pCmdUI, C2DBUILDER);
}

void CMainFrame::OnLoadResults() {

	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {

		C2DResults *pRes = dynamic_cast<C2DResults*>(pCurDoc);
		CEMMADoc *pOperDoc = pCurDoc->GetParent();
		CGen2DOperation *pOper = dynamic_cast<CGen2DOperation*>(pOperDoc);

		size_t nSize = pRes->GetSubDocNum();
		if (nSize > 0) {
			
			CloseAllChild();
			//m_wndExplorerPane.SetDocument(pRes);
			m_wndExplorerPane.OnRemoveAllSubDoc();

		}
		else {
			size_t nNum = pRes->LoadResults();		//nNum for debug
		}
		
		m_wndExplorerPane.FillTree();
	}

}

void CMainFrame::OnClearResults() {
	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {

		CResults *pRes = dynamic_cast<CResults*>(pCurDoc);

		size_t nSize = pRes->GetSubDocNum();
		if (nSize > 0) {

			CloseAllChild();
			for (size_t i = 0; i < m_DockBars.size(); i++) {
				m_DockBars[i]->ShowPane(0, 0, 0);
			}
			
			// Удаляем все загруженные результаты из дерева
			m_wndExplorerPane.OnRemoveAllSubDoc();

		}
	}
}

void CMainFrame::OnUpdateResultsDoc(CCmdUI *pCmdUI) {
	RibbonButtonOnOff(pCmdUI, C2DRESULTS);
}

//Активирует или деактивирует кнопки Ribbon в зависимости от типа документа
void CMainFrame::RibbonButtonOnOff(CCmdUI *pCmdUI, TYPEID nId) {
	CEMMADoc *pCurDoc = GetCurDoc();
	if (pCurDoc) {
		bool bType = pCurDoc->GetType() == nId;
		pCmdUI->Enable(bType);
	}
	else {
		pCmdUI->Enable(false);
	} 
	//*/
}

//////////////////////////////////////////////////////////////////////////
void CMainFrame::CloseDocument() {
	CloseAllChild();
	if (m_wndExplorerPane.GetSafeHwnd()) {
		m_wndExplorerPane.DestroyWindow();
	}
	DestroyDocument();
	RecalcLayout();	//Called by the framework when the standard control bars are toggled on or off or when the frame window is resized
}

//////////////////////////////////////////////////////////////////////////
void CMainFrame::DestroyDocument() {
	if (m_pDoc == nullptr) return;
	delete m_pDoc;
	m_pDoc = nullptr;
}

// Close all CChildFrame
void CMainFrame::CloseAllChild() {
	while (m_CFList.GetHeadPosition() != nullptr) {
		CChildFrame *pFrame = m_CFList.GetTail();
		if (pFrame) {
			pFrame->DestroyWindow();
			m_CFList.RemoveAt(m_CFList.GetTailPosition());
		}
	}

	/*
	while (m_CFList.front() != nullptr) {
	CChildFrame *pFrame = m_CFList.back();
	if (pFrame) {
	pFrame->DestroyWindow();
	m_CFList.erase(m_CFList.end());
	}
	}
	*/

	ShowControlBars();
}

//////////////////////////////////////////////////////////////////////////
void CMainFrame::ShowControlBars() {
	
	CChildFrame *pChild = dynamic_cast<CChildFrame *>(GetActiveFrame());
	if (!pChild) {
		return;
	}
	bool bNoView = !m_CFList.Find(pChild);	//если не нашли окно CChildFrame
	size_t ToolBarsSize = m_ToolBars.size();
	
	for (size_t i = 0; i<ToolBarsSize; i++) {
		/*	//no more MAIN toolbar (2016.09)
		if (g_ToolBarIDs[i] == IDR_MAIN_TB) {
			ShowPane(m_ToolBars[i], 1, 1, 0);
		}
		//*/
		//else 
		if (!bNoView) {
			ShowPane(m_ToolBars[i], pChild->ShowToolBar(g_ToolBarIDs[i]), 1, 0);
		}
		else {
			ShowPane(m_ToolBars[i], 0, 1, 0);
		}
	}

	HideAllPanels(); //Hide all panels\dockbars

	size_t DockBarsSize = m_DockBars.size();
	for (size_t i = 0; i < DockBarsSize; i++) {
		if (!bNoView) {
			bool bVisiblePane = pChild->IsPaneVisible(m_DockBars[i]->GetResID());
			ShowPane(m_DockBars[i], bVisiblePane, 1, 0);
			//m_DockBars[i]->Invalidate();
			m_DockBars[i]->SetDocument(pChild->GetDocument());
			
			//если панель видна, то обновляем\заполняем её контролы значениями
			if (bVisiblePane) {
				m_DockBars[i]->UpdatePane();
				m_DockBars[i]->Invalidate();
			}
		} else {
			ShowPane(m_DockBars[i], 0, 1, 0);
			m_DockBars[i]->SetDocument(nullptr);
		}
	}
	
	if (pChild->m_pView) {
		pChild->m_pView->Invalidate();
	}
	/*
	if (pChild->GetActiveView()) {
		pChild->GetActiveView()->Invalidate();
	} */

	//Invalidate();
	//RedrawWindow(nullptr, nullptr, RDW_ALLCHILDREN | RDW_INVALIDATE | RDW_UPDATENOW | RDW_FRAME | RDW_ERASE);
}

void CMainFrame::HideAllPanels() {
	size_t DockBarsSize = m_DockBars.size();
	for (size_t i = 0; i < DockBarsSize; i++) {
		m_DockBars[i]->Invalidate();
		ShowPane(m_DockBars[i], 0, 1, 0);	//hide all, visibility == 0
		m_DockBars[i]->Invalidate();
	}
}

//////////////////////////////////////////////////////////////////////////
void CMainFrame::CreateExplorerPane() {
	if (!m_pDoc) {
		return;		//Если нет документа, то нет и эксплорера
	}

	if (m_wndExplorerPane.GetSafeHwnd()) {
		m_wndExplorerPane.DestroyWindow();
	}

	CString strExplorer;
	int bNameValid = strExplorer.LoadString(IDS_EXPLORER_VIEW);

	if (!bNameValid) {
		CDlgShowError cError(ID_ERROR_MAINFRM_WRONG_EXPLORER_NAME); //_T("Wrong name IDS_EXPLORER_VIEW"));
		return;
	}

	if (!m_wndExplorerPane.Create(strExplorer, this, CRect(0, 0, 200, 200), 1 , m_wndExplorerPane.GetResID(), WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_CLIPCHILDREN | CBRS_LEFT | CBRS_FLOAT_MULTI))
	{
		CDlgShowError cError(ID_ERROR_MAINFRM_EXPLORERPANE_FAIL); //_T("Failed to create Class View window\n"));
		return; // failed to create
	}
	m_wndExplorerPane.SetDocument(m_pDoc);
	m_wndExplorerPane.EnableDocking(CBRS_ALIGN_LEFT);
	m_wndExplorerPane.ExpandAll();	//Раскрываем все подэлементы

	DockPane(&m_wndExplorerPane);
}

//////////////////////////////////////////////////////////////////////////
CMDIChildWnd* CMainFrame::CreateNewChild(CEMMADoc *pDoc) {
	CCreateContext context;
	context.m_pCurrentFrame = this;
	context.m_pNewViewClass = pDoc->GetViewClass();

	if (context.m_pNewViewClass == 0) {
		return nullptr;
	}

	CChildFrame* pFrame = new CChildFrame;
	m_CFList.AddTail(pFrame);
	// load the frame
	pFrame->SetHandles(pDoc->GetHMenu(), pDoc->GetHAccel());
	if (!pFrame->LoadFrame(pDoc->GetResID(), WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE, nullptr, &context)) {
		
		CDlgShowError cError(ID_ERROR_MAINFRM_FRAME_LOAD_FAIL); //_T("Couldn't load frame window. \n"));
		m_CFList.RemoveTail();
		return nullptr;
	}

	if (pFrame->m_pView) {
		pFrame->m_pView->SetDocument(pDoc);
	}

	// redraw the frame and parent
	pFrame->SetTitle(pDoc->GetName());
	pFrame->InitialUpdateFrame(nullptr, 1);
	pFrame->MDIActivate();//Активируем соответствующий фрейм
	return pFrame;
}

//////////////////////////////////////////////////////////////////////////
void CMainFrame::CreateDockPanes() {

	size_t nDocksSize = m_DockBars.size();
	for (size_t i = 0; i<nDocksSize; i++) {
		CString strName;
		UNINT uResID = m_DockBars[i]->GetResID();	//для наглядности
		int bNameValid = strName.LoadString(uResID);
		
		if (!bNameValid) {
			CDlgShowError cError(ID_ERROR_MAINFRM_DOCBAR_NAME_WRONG); //_T("Wrong strName in GetResID()"));
			return;
		}
		
		if (!m_DockBars[i]->Create(strName,
			this,
			CRect(0, 0, 200, 200),
			1,
			uResID,
			WS_CHILD | WS_CLIPSIBLINGS | WS_CLIPCHILDREN | CBRS_RIGHT | CBRS_FLOAT_MULTI))
		{
			CDlgShowError cError(ID_ERROR_MAINFRM_DOCBARS_PANES_VIEW_FAIL); //_T("Failed to create Class View window. \n"));
			return; // failed to create
		}
		m_DockBars[i]->EnableDocking(CBRS_ALIGN_ANY);
		DockPane(m_DockBars[i]);
	}

}

////////////////////////////////////////////////////////////////////////////////////////
// CMainFrame message handlers
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CMDIFrameWndEx::OnCreate(lpCreateStruct) == -1)
		return -1;


	// set the visual manager and style based on persisted value
	OnApplicationLook(theApp.m_nAppLook);

	CMDITabInfo mdiTabParams;
	mdiTabParams.m_style = CMFCTabCtrl::STYLE_3D_ONENOTE; // другие доступные стили...
	mdiTabParams.m_bActiveTabCloseButton = 1;   // установите значение 0, чтобы расположить кнопку \"Закрыть\" в правой части области вкладки
	mdiTabParams.m_bTabIcons = 0;  // установите значение 1, чтобы включить значки документов на вкладках MDI
	mdiTabParams.m_bAutoColor = 1;  // установите значение 0, чтобы отключить автоматическое выделение цветом вкладок MDI
	mdiTabParams.m_bDocumentMenu = 1; // включить меню документа на правой границе области вкладки
	EnableMDITabbedGroups(1, mdiTabParams);

	/*
	if (!m_wndMenuBar.Create(this))
	{
		//CDlgShowError cError(_T("Failed to create menubar. \n"));
		return -1;   // fail to create
	}

	m_wndMenuBar.SetPaneStyle(m_wndMenuBar.GetPaneStyle() | CBRS_SIZE_DYNAMIC | CBRS_TOOLTIPS | CBRS_FLYBY);
	*/

	if (!m_wndStatusBar.Create(this)) {
		CDlgShowError cError(ID_ERROR_MAINFRM_STATUSBAR_FAIL); //_T("Failed to create status bar. \n"));
		return -1;   // fail to create
	}
	m_wndStatusBar.SetIndicators(g_indicators, sizeof(g_indicators) / sizeof(UNINT));

	// prevent the menu bar from taking the focus on activation
	CMFCPopupMenu::SetForceMenuFocus(0);

	// enable Visual Studio 2005 style docking window behavior
	CDockingManager::SetDockingMode(DT_SMART);
	// enable Visual Studio 2005 style docking window auto-hide behavior
	EnableAutoHidePanes(CBRS_ALIGN_ANY);

	m_wndRibbonBar.Create(this);
	m_wndRibbonBar.LoadFromResource(IDR_RIBBON);

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
void CMainFrame::OnApplicationLook(UNINT id)
{
	CWaitCursor wait;

	theApp.m_nAppLook = id;

	switch (theApp.m_nAppLook)
	{
	case ID_VIEW_APPLOOK_WIN_2000:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManager));
		m_wndRibbonBar.SetWindows7Look(0);
		break;

	case ID_VIEW_APPLOOK_OFF_XP:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerOfficeXP));
		m_wndRibbonBar.SetWindows7Look(0);
		break;

	case ID_VIEW_APPLOOK_WIN_XP:
		CMFCVisualManagerWindows::m_b3DTabsXPTheme = 1;
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerWindows));
		m_wndRibbonBar.SetWindows7Look(0);
		break;

	case ID_VIEW_APPLOOK_OFF_2003:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerOffice2003));
		CDockingManager::SetDockingMode(DT_SMART);
		m_wndRibbonBar.SetWindows7Look(0);
		break;

	case ID_VIEW_APPLOOK_VS_2005:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerVS2005));
		CDockingManager::SetDockingMode(DT_SMART);
		m_wndRibbonBar.SetWindows7Look(0);
		break;

	case ID_VIEW_APPLOOK_WINDOWS_7:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerWindows7));
		CDockingManager::SetDockingMode(DT_SMART);
		m_wndRibbonBar.SetWindows7Look(1);
		break;

	default:
		switch (theApp.m_nAppLook)
		{
		case ID_VIEW_APPLOOK_OFF_2007_BLUE:
			CMFCVisualManagerOffice2007::SetStyle(CMFCVisualManagerOffice2007::Office2007_LunaBlue);
			break;

		case ID_VIEW_APPLOOK_OFF_2007_BLACK:
			CMFCVisualManagerOffice2007::SetStyle(CMFCVisualManagerOffice2007::Office2007_ObsidianBlack);
			break;

		case ID_VIEW_APPLOOK_OFF_2007_SILVER:
			CMFCVisualManagerOffice2007::SetStyle(CMFCVisualManagerOffice2007::Office2007_Silver);
			break;

		case ID_VIEW_APPLOOK_OFF_2007_AQUA:
			CMFCVisualManagerOffice2007::SetStyle(CMFCVisualManagerOffice2007::Office2007_Aqua);
			break;
		}

		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerOffice2007));
		CDockingManager::SetDockingMode(DT_SMART);
		m_wndRibbonBar.SetWindows7Look(0);
	}

	//RedrawWindow(nullptr, nullptr, RDW_ALLCHILDREN | RDW_INVALIDATE | RDW_UPDATENOW | RDW_FRAME | RDW_ERASE);

	CString str_AppLook;
	int nLoadStr = str_AppLook.LoadString(IDS_APP_LOOK);
	theApp.WriteInt(str_AppLook, theApp.m_nAppLook);
}

/////////////////////////////////////////////////////////////////////////////////////////////
void CMainFrame::OnUpdateApplicationLook(CCmdUI* pCmdUI) {
	bool bAppLook = (theApp.m_nAppLook == pCmdUI->m_nID);	//для отладки
	pCmdUI->SetRadio(bAppLook);
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Новый документ
void CMainFrame::OnFileNew() {

	//Вызвать диалог выбора типа задачи
	CDlgNewTask DlgNewTask;
	if (DlgNewTask.DoModal() != IDOK) return;

	//Уничтожаем старый документ
	CloseDocument();

	//Создаем объект документа, соответствующий выбранному типу задачи
	//и сохранить адрес созданного объекта в m_pDoc
	CString strTaskName;	//название задач берём из ресурсов
	int nLoadStr = -1;

	size_t nSel = DlgNewTask.GetSel();
	switch (nSel) {
	case IDS_GEN2D:// 2D
	{
		m_pDoc = new CGen2DDoc;
		nLoadStr = strTaskName.LoadString(IDS_GEN2D);
		m_pDoc->SetName(strTaskName);
		dynamic_cast<CGen2DDoc*>(m_pDoc)->AddOperationOnConstructor();
	}
	break;
	case IDS_FORMING:// формовка
	{
		m_pDoc = new C2DFormingDoc;
		nLoadStr = strTaskName.LoadString(IDS_FORMING);
		m_pDoc->SetName(strTaskName);
	}
	break;
	case IDS_ROLLING:// Прокатка
	{
		m_pDoc = new CRollingDoc;
		nLoadStr = strTaskName.LoadString(IDS_ROLLING);
		m_pDoc->SetName(strTaskName);
	}
	break;
	case IDS_DRAGGING:// волочение
	{
		m_pDoc = new C2DDraggingDoc;
		nLoadStr = strTaskName.LoadString(IDS_DRAGGING);
		m_pDoc->SetName(strTaskName);
	}
	break;
	case IDS_GEN3D:// 3D
	{
		m_pDoc = new CGen3DDoc;
		nLoadStr = strTaskName.LoadString(IDS_GEN3D);
		m_pDoc->SetName(strTaskName);
	}
	break;
	case IDS_BALLOON:
	{
		m_pDoc = new C2DBalloonDoc;
		nLoadStr = strTaskName.LoadString(IDS_BALLOON);
		m_pDoc->SetName(strTaskName);
	}
	break;
	default:
		CDlgShowError cError(ID_ERROR_MAINFRM_TASKTYPE_WRONG); //_T("Wrong Task type"));
		return;
	}
	//Создаем левую боковую панель (и дерево)
	CreateExplorerPane();
}

//////////////////////////////////////////////////////////////////////////
void CMainFrame::OnFileOpen() {

	//LOGGER.Init(CString("..\\..\\Logs\\CMainFrame.OnFileOpen.txt"));
	CFileDialog fileDlg( 1 , _T("emma"), 0, 4 | 2, _T("emma files (*.emma)|*.emma|"));
	if (fileDlg.DoModal() == IDOK) {
		CloseDocument();	//Уничтожаем старый документ после выбора нового
		CString strFile = fileDlg.GetPathName();
		CStorage file;
		if (!file.Open(strFile, CStorage::modeRead | CStorage::typeBinary)) {
			CDlgShowError(ID_ERROR_MAINFRM_LOADFILE_WRONG); //_T("Невозможно загрузить файл"));
			return;
		};
		m_pDoc = dynamic_cast<CEMMADoc*>(IO::LoadObject(file));
	}
	else {
		return;
	}

	if (m_pDoc) {
		//Создаем левую боковую панель, в ней - древовидную структуру
		CreateExplorerPane();
	}
	else {
		CDlgShowError(ID_ERROR_MAINFRM_LOADFILE_WRONG); //_T("Невозможно загрузить файл"));
	}
	
}

//////////////////////////////////////////////////////////////////////////
void CMainFrame::OnFileSave() {
	if (!m_pDoc)return;

	CString strFile = m_pDoc->GetFilePath();
	if (!strFile.IsEmpty()) {
		CStorage file;
		int bOpen = file.Open(strFile, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeBinary);
		bool bSave = m_pDoc->Save(file);
	}else {
		//Если сохраняем первый раз, то необходимо указать путь до файла
		OnFileSaveAs();
	}

}

void CMainFrame::OnUpdateFileSave(CCmdUI* pCmdUI) {
	if (m_pDoc == nullptr) {
		pCmdUI->Enable(false);
	}else{
		pCmdUI->Enable(true);
	}
}

void CMainFrame::OnFileSaveAs() {
	if (!m_pDoc)return;
	LOGGER.Init(CString("..\\..\\Logs\\CMainFrame.OnFileSave.txt"));
	CFileDialog fileDlg(0, _T("emma"), 0, 4 | 2, _T("emma files (*.emma)|*.emma|"));
	CString strFile;
	bool bSave = false;
	if (fileDlg.DoModal() == IDOK) {
		strFile = fileDlg.GetPathName();
		CStorage file;
		int bOpen = file.Open(strFile, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeBinary);
		bSave = m_pDoc->Save(file);
	}
	//Если успешно сохранили в новый файл, записываем путь до него, чтобы сохранять потом в него же
	if (bSave) {
		m_pDoc->SetFilePath(strFile);
	}
}

void CMainFrame::OnUpdateFileSaveAs(CCmdUI* pCmdUI) {
	if (m_pDoc == nullptr) {
		pCmdUI->Enable(false);
	}
	else {
		pCmdUI->Enable(true);
	}
}

void CMainFrame::OnFileClose() {
	CloseDocument();
}

void CMainFrame::OnUpdateFileClose(CCmdUI* pCmdUI) {
	if (m_pDoc == nullptr) {
		pCmdUI->Enable(false);
	}
	else {
		pCmdUI->Enable(true);
	}
}

//////////////////////////////////////////////////////////////////////////
LRESULT CMainFrame::OnChildFrm(WPARAM wParam, LPARAM lParam) {
	HTREEITEM hTreeItem;
	if (wParam) {
		CChildFrame *pChild = (CChildFrame *)wParam;
		hTreeItem = pChild->GetTreeItem();
		/*if (pChild->GetActiveView()) {
			pChild->GetActiveView()->Invalidate();
		}*/
	}
	else {
		//hTreeItem = m_pDoc ? m_pDoc->GetTreeItem() : nullptr;
		hTreeItem = nullptr;
		//HideAllPanels();	//hide previous panels
		//UpdateAllViews();	//
	}
	m_wndExplorerPane.SetCurSel(hTreeItem);
	//Меняем активную категорию
	SetActiveRibbonCategory();
	ShowControlBars();
	return 0;
}

//////////////////////////////////////////////////////////////////////////
LRESULT CMainFrame::OnChildClose(WPARAM wParam, LPARAM lParam) {
	CChildFrame *pChild = (CChildFrame *)wParam;
	POSITION pos = m_CFList.Find(pChild);

	if (pos == nullptr) {
		CDlgShowError cError(ID_ERROR_MAINFRM_POSITION_NOT_FOUND); //_T("Позиция не найдена m_CFList.Find(pChild)"));
		return 0;
	}
	pChild->CloseWindow();
	//если нашли, то удаляем
	m_CFList.RemoveAt(pos);
	UpdateAllViews();
	return 0;
}

//////////////////////////////////////////////////////////////////////////
LRESULT CMainFrame::OnEPSelChange(WPARAM wParam, LPARAM lParam) {
	CEMMADoc *pDoc = (CEMMADoc *)wParam;
	CWnd *pView = pDoc->GetView();	//для отладки

	if (pView != nullptr) {
		// Если представление у документа есть, ищем его
		POSITION pos = m_CFList.GetHeadPosition();
		bool bFound = false;
		while (pos != nullptr) {
			CChildFrame *pFrame = m_CFList.GetNext(pos);
			if (pFrame->m_pView == pView) {//Представление найдено

				// Эти строчки необходимы, если у нескольких документов одно представление на всех {{
				CEMMADoc *pOldDoc = pFrame->m_pView->GetDocument();

				if (pOldDoc && pOldDoc != pDoc)
				{
					pOldDoc->Deactivate();

					//Проверяем панели Результатов (2D и 3D вид)
					if (pDoc->IsPaneVisible(IDR_2DCALCSTAGE_TB) || pDoc->IsPaneVisible(IDR_3DCALCSTAGE_TB)) {
						for (size_t i = 0; i < m_DockBars.size(); i++) {
							if (m_DockBars[i]->GetResID() == IDR_2DCALCSTAGE_PANE || m_DockBars[i]->GetResID() == IDR_3DCALCSTAGE_PANE) {
								m_DockBars[i]->SetDocument(pDoc);
								break;
							}
						}
					}
				}

				(dynamic_cast<CBasicView *>(pView))->SetDocument(pDoc);
				pFrame->SetTitle(pDoc->GetName());
				pFrame->InitialUpdateFrame(nullptr, 1);	// }}

				pFrame->MDIActivate();//Активируем соответствующий фрейм
				bFound = true;
				break;
			}
		}

		if (!bFound) {
			CDlgShowError cError(ID_ERROR_MAINFRM_VIEW_NOT_FOUND); //_T("Представление не найдено"));
			return 0;
		}
	}
	else {					 // Если представления у документа нет
		CreateNewChild(pDoc);// Создаем для него новый фрейм и представление
	}

	//Меняем активную категорию
	SetActiveRibbonCategory();

	return 0;
}

//! Меняем активные категории\вкладки у Ribbon панели
void CMainFrame::SetActiveRibbonCategory() {
	CEMMADoc *pDoc = GetCurDoc();

	if (!pDoc) {
		//CDlgShowError cError();
		return;
	}

	TYPEID type = pDoc->GetType();
	int nCategoryIndex = 1;	//первая вкладка всегда показывается
	
	//Результаты - категория\вкладка 2 (Вид)
	if (type == C2DCALCSTAGE || type == C3DCALCSTAGE) {
		nCategoryIndex = 2;
	}

	//Чертёж, Инструмент, Траектория - категория\вкладка 3 (Редактор)
	if (type == C2DSKETCH || type == C2DTOOL || type == C2DTRAJECTORY) {
		nCategoryIndex = 3;
	}

	//КЭ сетка - категория\вкладка 4 (КЭ сетка)
	if (type == C2DMESHER) {
		nCategoryIndex = 4;
	}

	//Cборка проекта и результаты - категория\вкладка 5 (Проект)
	if (type == C2DBUILDER || type == C2DRESULTS) {
		nCategoryIndex = 5;
	}

	CMFCRibbonCategory* pCategory = m_wndRibbonBar.GetCategory(nCategoryIndex);
	if (!pCategory) {
		//CDlgShowError cError();
		return;
	}
	int nRes = m_wndRibbonBar.SetActiveCategory(pCategory, 1);
}

//////////////////////////////////////////////////////////////////////////
void CMainFrame::UpdateAllViews() {
	POSITION pos = m_CFList.GetHeadPosition();
	
	while (pos != nullptr) {
		CChildFrame *pFrame = m_CFList.GetNext(pos);
		pFrame->m_pView->Invalidate();
		pFrame->m_pView->UpdateWindow();
		//pFrame->InitialUpdateFrame(nullptr, 1);	// }}
	}
}

//////////////////////////////////////////////////////////////////////////
LRESULT CMainFrame::OnEPClose(WPARAM wParam, LPARAM lParam) {
	CloseDocument();
	return 0;
}

// 
LRESULT CMainFrame::OnReopenSaves(WPARAM wParam, LPARAM lParam) {

	//CloseAllChild();
	CEMMADoc *pDoc = (CEMMADoc *)wParam;
	//pDoc->m_arSubDocs

	//TODO: закрыть только избранное (CResults) окно

	/*while(m_CFList.GetHeadPosition() != nullptr){
	CChildFrame *pFrame = m_CFList.GetTail();
	if (pFrame->m_pView == pDoc->GetView()){ //pFrame
	pFrame->DestroyWindow();
	m_CFList.RemoveAt(m_CFList.GetTailPosition());
	}
	}*/
	POSITION pos = m_CFList.GetHeadPosition();
	
	while (pos != nullptr) {
		CChildFrame *pFrame = m_CFList.GetNext(pos);
		if (pFrame->m_pView == pDoc->GetView()) {//Представление найдено
			pFrame->DestroyWindow();
			//m_CFList.RemoveAt(pos);
			break;
		}
	}
	ShowControlBars();

	return 0;
}

//////////////////////////////////////////////////////////////////////////
LRESULT CMainFrame::OnSubDocDestroy(WPARAM wParam, LPARAM lParam) {
	CEMMADoc *pDoc = (CEMMADoc *)wParam;

	if (pDoc->GetView()) {// Если представление у документа есть, ищем его
		POSITION pos = m_CFList.GetHeadPosition();
		
		while (pos != nullptr) {
			CChildFrame *pFrame = m_CFList.GetNext(pos);
			if (pFrame->m_pView == pDoc->GetView()) {//Представление найдено
				pFrame->SendMessage(WM_CLOSE);//Закрываем соответствующий фрейм
				break;
			}
		}
	}

	return 0;
}



//////////////////ROMA
LRESULT CMainFrame::OnUpdatePane(WPARAM wParam, LPARAM lParam) {
	CChildFrame *pChild = dynamic_cast<CChildFrame *>(GetActiveFrame());
	if (!pChild) return 0;

	size_t nSize = m_DockBars.size();
	for (size_t i = 0; i<nSize; i++) {
		UNINT nResID = m_DockBars[i]->GetResID();
		if (pChild->IsPaneVisible(nResID)) {
			m_DockBars[i]->UpdatePane();
		}
	}

	return 0;
}

/*
int CMainFrame::OnEraseBkgnd(CDC* pDC)
{
	int nErase = CMDIFrameWndEx::OnEraseBkgnd(pDC);
	return nErase;
}

void CMainFrame::OnPaint()
{
	CPaintDC dc(this); // device context for painting
					   //Перекраска фона в белый цвет при изменении размера
	CRect rectClient;
	GetClientRect(&rectClient);

	CBrush bg;
	bg.CreateStockObject(WHITE_BRUSH);
	dc.FillRect(&rectClient, &bg);
}
*/
////////////////////////

//////////////////////////////////////////////////////////////////////////
// CMainFrame diagnostics
//////////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CMDIFrameWndEx::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CMDIFrameWndEx::Dump(dc);
}
#endif //_DEBUG


void CMainFrame::OnClose() {
	//TODO: Если запущена в режиме отладки, выходить по кнопке Stop!!!
	CMDIFrameWndEx::OnClose();
}
