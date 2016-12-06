#include "StdAfx.h"
#include "2DSketchPane.h"


C2DSketchPane::~C2DSketchPane(void)
{
}

BEGIN_MESSAGE_MAP(C2DSketchPane, CEMMARightPane)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_WM_PAINT()

	////Roma
	ON_REGISTERED_MESSAGE(AFX_WM_PROPERTY_CHANGED, C2DSketchPane::OnPropertyChanged)

	////TOOLBAR
	ON_COMMAND(ID_CLEAR_SKETCH, &C2DSketchPane::OnClearSketch)
	ON_UPDATE_COMMAND_UI(ID_CLEAR_SKETCH, &C2DSketchPane::OnUpdateClearSketch)

	ON_COMMAND(ID_NEW_NODE, &C2DSketchPane::OnNewNode)
	ON_UPDATE_COMMAND_UI(ID_NEW_NODE, &C2DSketchPane::OnUpdateNewNode)

	ON_COMMAND(ID_NEW_CURVE, &C2DSketchPane::OnNewCurve)
	ON_UPDATE_COMMAND_UI(ID_NEW_CURVE, &C2DSketchPane::OnUpdateNewCurve)

	ON_COMMAND(ID_NEW_CIRCLE_ARC, &C2DSketchPane::OnNewCircleArc)
	ON_UPDATE_COMMAND_UI(ID_NEW_CIRCLE_ARC, &C2DSketchPane::OnUpdateNewCircleArc)

	ON_COMMAND(ID_FACET, &C2DSketchPane::OnNewFacet)
	ON_UPDATE_COMMAND_UI(ID_FACET, &C2DSketchPane::OnUpdateNewFacet)

	ON_BN_CLICKED(ID_2DSKETCHPANE_BUTTON_SHOWCONTOURS, C2DSketchPane::ShowContours)
	ON_UPDATE_COMMAND_UI(ID_2DSKETCHPANE_BUTTON_SHOWCONTOURS, C2DSketchPane::OnEnableButton)

	ON_BN_CLICKED(ID_2DSKETCHPANE_BUTTON_SAVE_TO_META, C2DSketchPane::SaveSketchToMeta)
	ON_UPDATE_COMMAND_UI(ID_2DSKETCHPANE_BUTTON_SAVE_TO_META, C2DSketchPane::OnEnableButtonSave)

END_MESSAGE_MAP()

int C2DSketchPane::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CEMMADockablePane::OnCreate(lpCreateStruct) == -1)
		return -1;

	// Создание Тулбара {{
	m_wndToolBar.Create(this, AFX_DEFAULT_TOOLBAR_STYLE, IDR_2DSKETCH_TB);
	m_wndToolBar.LoadToolBar(IDR_2DSKETCH_TB, 0, 0, 1 /* Is locked */);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() | CBRS_TOOLTIPS | CBRS_FLYBY);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() & ~(CBRS_GRIPPER | CBRS_SIZE_DYNAMIC | CBRS_BORDER_TOP | CBRS_BORDER_BOTTOM | CBRS_BORDER_LEFT | CBRS_BORDER_RIGHT));
	
	m_wndToolBar.SetOwner(this);
	// All commands will be routed via this control , not via the parent frame:
	m_wndToolBar.SetRouteCommandsViaFrame(0);
	// }}


	CRect rectDummy;
	rectDummy.SetRectEmpty();

	if (!m_wndPropList.Create(WS_CHILD | WS_VISIBLE, rectDummy, this, 1))
	{
		CDlgShowError cError(ID_ERROR_2DSKETCHPANE_FAIL_PROPGRID); //_T("Failed to create Properties Grid \n"));
		return -1;      // fail to create
	}

	CString str4button_showcontours;
	str4button_showcontours.LoadStringW(ID_2DSKETCHPANE_BUTTON_SHOWCONTOURS);
	if (!m_button_contour.Create(str4button_showcontours, WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON, rectDummy, this, ID_2DSKETCHPANE_BUTTON_SHOWCONTOURS)) {
		CDlgShowError cError(ID_ERROR_2DSKETCHPANE_FAIL_BUTTON); //_T("Failed to create Button \n"));
		return -1;      // fail to create
	}
	//устанавливаем такой же шрифт, как в таблице свойств
	m_button_contour.SetFont(m_wndPropList.GetFont());

	CString str4button_save2meta;
	str4button_save2meta.LoadStringW(ID_2DSKETCHPANE_BUTTON_SAVE_TO_META);
	if (!m_button_save_to_meta.Create(str4button_save2meta, WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON, rectDummy, this, ID_2DSKETCHPANE_BUTTON_SAVE_TO_META)) {
		CDlgShowError cError(ID_ERROR_2DSKETCHPANE_FAIL_BUTTON); //_T("Failed to create Button \n"));
		return -1;      // fail to create
	}
	//устанавливаем такой же шрифт, как в таблице свойств
	m_button_save_to_meta.SetFont(m_wndPropList.GetFont());

	//! Устанавливаем имена колонок Таблицы
	SetColumnNames();
	//! Устанавливаем размеры Тулбара и Таблицы и Кнопок
	AdjustLayout();

	m_cursor = 0;

	return 0;
}

//! Устанавливаем размеры Тулбара и Таблицы и Кнопок
void C2DSketchPane::AdjustLayout(void){
	if (GetSafeHwnd() == nullptr) {
		return;
	}

	CRect rectClient;
	GetClientRect(rectClient);

	int cyTlb = m_wndToolBar.CalcFixedLayout(0, 1).cy;

	//устанавливает размеры тулбара (справа)
	m_wndToolBar.SetWindowPos(nullptr, rectClient.left, rectClient.top, rectClient.Width(), cyTlb, SWP_NOACTIVATE | SWP_NOZORDER);
	//устанавливает размеры Таблицы свойств (справа)
	m_wndPropList.SetWindowPos(nullptr, rectClient.left, rectClient.top + cyTlb, rectClient.Width(), (rectClient.Height() - cyTlb)/2, SWP_NOACTIVATE | SWP_NOZORDER);
	
	//устанавливаем размеры Кнопки
	int cyPropList = rectClient.top + cyTlb + (rectClient.Height() - cyTlb) / 2 + 10;
	int nFromLeft = rectClient.left + int(rectClient.Width() / 5.0);
	int nFromRight = rectClient.Width() - int(rectClient.Width() / 5.0)*2;	//обязательное умножение, чтобы был отступ справа
	m_button_contour.SetWindowPos(nullptr, nFromLeft, cyPropList, nFromRight, 20, SWP_NOOWNERZORDER | SWP_SHOWWINDOW);
	cyPropList += 20;
	m_button_save_to_meta.SetWindowPos(nullptr, nFromLeft, cyPropList, nFromRight, 20, SWP_NOOWNERZORDER | SWP_SHOWWINDOW);
}


////////TOOLBAR
void C2DSketchPane::OnClearSketch(){

	//GetSketch()->GetOutline()->ClearAll();	//включено в SetClearOutline
	GetSketch()->SetClearOutline();
	UpdatePane();
}

void C2DSketchPane::OnUpdateClearSketch(CCmdUI *pCmdUI){
	// TODO: Add your command update UI handler code here
}


void C2DSketchPane::OnNewNode(){
	
	if(m_cursor != 100){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 100);
		m_cursor = 100;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}

}
void C2DSketchPane::OnUpdateNewNode(CCmdUI *pCmdUI){

	if(m_cursor == 100) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DSketchPane::OnNewCurve(){
	
	if(m_cursor != 101){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 101);
		m_cursor = 101;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
}
void C2DSketchPane::OnUpdateNewCurve(CCmdUI *pCmdUI){

	if(m_cursor == 101) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DSketchPane::OnNewCircleArc(){
	
	if(m_cursor != 102){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 102);
		m_cursor = 102;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
}
void C2DSketchPane::OnUpdateNewCircleArc(CCmdUI *pCmdUI){

	if(m_cursor == 102) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DSketchPane::OnNewFacet(){
	
	if(m_cursor != 103){
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR, 103);
		m_cursor = 103;
	}
	else{
		m_pDoc->GetView()->SendMessage(WMR_SETCURSOR,0);
		m_cursor = 0;
	}
}

void C2DSketchPane::OnUpdateNewFacet(CCmdUI *pCmdUI){

	if(m_cursor == 103) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void C2DSketchPane::ShowContours() {
	C2DOutline *pOutline = GetSketch()->GetOutline();
	//Две точки могут образовать контур (две дуги из двух точек с разными радиусами), меньше - нет
	if (pOutline->GetNodeCount() > 1) {
		//TODO: выделение всех контуров

	}
	else {
		CDlgShowError Diag_err(ID_ERROR_2DSKETCHPANE_NO_CONTOURS); //_T("Contours not found."));	//Показываем окно ошибки
	}
}

void C2DSketchPane::OnEnableButton(CCmdUI *pCmdUI) {
	//Приходится активировать кнопку для нажатий (MFC style)
	pCmdUI->Enable();
}

//Сохранение всего чертежа в meta-файл
void C2DSketchPane::SaveSketchToMeta() {
	C2DOutline *pOutline = GetSketch()->GetOutline();
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

void C2DSketchPane::OnEnableButtonSave(CCmdUI *pCmdUI) {
	//Приходится активировать кнопку для нажатий (MFC style)
	pCmdUI->Enable();
}
