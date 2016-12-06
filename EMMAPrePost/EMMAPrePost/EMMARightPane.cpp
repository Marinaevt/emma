#include "stdafx.h"
#include "EMMARightPane.h"

CEMMARightPane::~CEMMARightPane(void)
{
}

BEGIN_MESSAGE_MAP(CEMMARightPane, CEMMADockablePane)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_WM_PAINT()
	//Обновление всех видимых панелей
	//ON_REGISTERED_MESSAGE(AFX_WM_PROPERTY_CHANGED, OnPropertyChanged)

END_MESSAGE_MAP()

//! Устанавливаем заголовки для колонок таблицы
void CEMMARightPane::SetColumnNames(void) {

	CString sParam, sValue;
	int nNumPar = sParam.LoadString(IDS_PARAM);
	int nNumVal = sValue.LoadString(IDS_PARAM_VALUE);
	m_wndPropList.EnableHeaderCtrl( 1 , sParam, sValue);
}

int CEMMARightPane::OnEraseBkgnd(CDC* pDC) {
	return CWnd::OnEraseBkgnd(pDC);
}

void CEMMARightPane::OnPaint()
{
	CPaintDC dc(this); // device context for painting

	CRect rect;
	GetClientRect(&rect);	//Получаем размеры
							//Перекраска фона в белый цвет при изменении размера
	CBrush bg;
	bg.CreateStockObject(WHITE_BRUSH);
	dc.FillRect(&rect, &bg);

}

void CEMMARightPane::UpdatePane() {

	//Заполняем таблицу параметров из документа
	m_pDoc->FillPropList(&m_wndPropList);

}

void CEMMARightPane::OnSize(UNINT nType, int cx, int cy)
{
	CEMMADockablePane::OnSize(nType, cx, cy);
	AdjustLayout();
}

LRESULT CEMMARightPane::OnPropertyChanged(__in WPARAM wparam, __in LPARAM lparam)
{

	CMFCPropertyGridProperty * pProperty = (CMFCPropertyGridProperty *)lparam;
	COleVariant v = pProperty->GetValue();
	UNLONG id = pProperty->GetData();
	v.ChangeType(VT_R8);

	m_pDoc->UpdateProp(&v.dblVal, &id);
	//m_pDoc->UpdatePropList(&m_wndPropList);
	//UpdatePane();
	//Заполняем таблицу параметров из документа
	m_pDoc->FillPropList(&m_wndPropList);
	return 0;
}