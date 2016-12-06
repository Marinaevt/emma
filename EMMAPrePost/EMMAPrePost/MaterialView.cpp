#include "StdAfx.h"
#include "MaterialView.h"

BEGIN_MESSAGE_MAP(CMaterialView, CMaterialView_Parent)
	
	ON_WM_CREATE()
	ON_WM_PAINT()
	ON_WM_SIZE()

	ON_REGISTERED_MESSAGE(WMR_INITIALUPDATE, OnInitialUpdate)
	ON_WM_MOUSEWHEEL()
	ON_WM_MOUSEMOVE()
	
END_MESSAGE_MAP()

IMPLEMENT_DYNCREATE(CMaterialView, CMaterialView_Parent);

CMaterialView::CMaterialView(void) : CMaterialView_Parent()
{
}

CMaterialView::~CMaterialView(void)
{
}

int CMaterialView::OnCreate(LPCREATESTRUCT lpCreateStruct) {

	if (CMaterialView_Parent::OnCreate(lpCreateStruct) == -1)
		return -1;

	//Устанавливаем размер окна редактора в 0
	CRect rectDummy;
	rectDummy.SetRectEmpty();

	//создаём поле редактора
	if (!m_MatEditArea.Create(WS_CHILD | WS_VISIBLE | WS_TABSTOP | ES_AUTOVSCROLL | ES_AUTOHSCROLL | ES_MULTILINE | WS_HSCROLL | WS_VSCROLL,
		rectDummy, this, IDS_MATERIAL_EDIT)) {
		CDlgShowError cError(ID_ERROR_MATERIALVIEW_EDITAREA_FAIL); //_T("Failed to create CEdit \n"));
		return -1;      // fail to create
	}

	m_MatEditArea.SetOwner(this);

	//устанавливаем больший объём текста (2147483646 (0x7FFFFFFE))
	m_MatEditArea.SetLimitText(0);

	return 0;
}

void CMaterialView::OnSize(UNINT nType, int cx, int cy)
{
	//устанавливаем размер текстовой области
	if (::IsWindow(m_MatEditArea.m_hWnd)) {
		m_MatEditArea.MoveWindow(0, 0, cx, cy, 1);
	}
	CMaterialView_Parent::OnSize(nType, cx, cy);
}

void CMaterialView::SetDocument(CEMMADoc *pDoc) {
	CMaterialView_Parent::SetDocument(pDoc);

	FillEditArea(pDoc);
	Invalidate(0);
}

void CMaterialView::FillEditArea(CEMMADoc *pDoc) {

	//устанавливаем текст в окно области редактирования
	CMaterialDoc *pMatDoc = dynamic_cast<CMaterialDoc *>(pDoc);
	if (!pMatDoc) {
		CDlgShowError cError(ID_ERROR_MATERIALVIEW_MATERIALDOC_NULL); //_T("Failed to get CMaterialDoc \n"));
		return;
	}
	//если всё хорошо, то устанавливаем содержимое
	CString sArea = pMatDoc->GetMaterialFileContent();

	sArea.Replace(_T("\r"), _T(""));
	sArea.Replace(_T("\n"), _T("\r\n"));
	m_MatEditArea.SetWindowText(sArea);
	m_MatEditArea.SetFocus();
}

LRESULT CMaterialView::OnInitialUpdate(WPARAM wParam, LPARAM lParam) {
	if (!m_pDoc)return 0;
	//пусто, пока нечего обновлять
	return 0;	
}

void CMaterialView::OnPaint() {
	CPaintDC dc(this); // device context for painting
					   //Перекраска фона в белый цвет при изменении размера
	CRect rectClient;
	GetClientRect(&rectClient);

	CBrush bg;
	bg.CreateStockObject(WHITE_BRUSH);
	dc.FillRect(&rectClient, &bg);
}