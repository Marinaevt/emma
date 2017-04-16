// BasicView.cpp : implementation file
//
#include "stdafx.h"
#include "BasicView.h"
#include "EMMAPrePost.h"

// CBasicView
BEGIN_MESSAGE_MAP(CBasicView, CBasicView_Parent)
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_MOUSEWHEEL()
	ON_WM_MOUSELEAVE()
	ON_WM_MBUTTONDOWN()
	ON_WM_MBUTTONUP()
	//ON_WM_ERASEBKGND()
	ON_WM_DESTROY()
END_MESSAGE_MAP()
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
CBasicView::CBasicView(){
		m_pDoc = nullptr;
		m_bLButtonPressed = false;
		m_bRButtonPressed = false;
		m_bMButtonPressed = false;
		m_bMouseIsInside = false;
		m_MousePos = CPoint(0, 0);
		m_Parameter = nullptr;
}
////////////////////////////////////////////////////////////////////////////////
CBasicView::~CBasicView()
{
}

////////////////////////////////////////////////////////////////////////////////
void CBasicView::SetDocument(CEMMADoc *pDoc){
	m_pDoc = pDoc;

	if(m_pDoc!=nullptr){
		m_pDoc->SetView(this);
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int CBasicView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	//
	if (CBasicView_Parent::OnCreate(lpCreateStruct) == -1)return -1;
//	TrackMouseEvent(TME_HOVER | TME_LEAVE);// пока не знаем где мышь - ждем оба сообщения
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CBasicView::OnLButtonDown(UNINT nFlags, CPoint point)
{// Нажатие левой кнопки
	CBasicView_Parent::OnLButtonDown(nFlags, point);// Базовая реализация
	SetCapture();
	SetFocus();
	m_bLButtonPressed = true;
}
////////////////////////////////////////////////////////////////////////////////
void CBasicView::OnLButtonUp(UNINT nFlags, CPoint point)
{// Отпускание левой кнопки
	CBasicView_Parent::OnLButtonUp(nFlags, point);// Базовая реализация
	if(!m_bRButtonPressed && !m_bMButtonPressed) ReleaseCapture();
	m_bLButtonPressed = false;
}
////////////////////////////////////////////////////////////////////////////////
void CBasicView::OnRButtonDown(UNINT nFlags, CPoint point)
{// Нажатие правой кнопки
	CBasicView_Parent::OnRButtonDown(nFlags, point);// Базовая реализация
	SetCapture();
	m_bRButtonPressed = true;
}
////////////////////////////////////////////////////////////////////////////////
void CBasicView::OnRButtonUp(UNINT nFlags, CPoint point)
{// отпускание правой кнопки
	CBasicView_Parent::OnRButtonUp(nFlags, point);// Базовая реализация
	if(!m_bLButtonPressed && !m_bMButtonPressed) ReleaseCapture();
	m_bRButtonPressed = false;
}
////////////////////////////////////////////////////////////////////////////////
void CBasicView::OnMButtonDown(UNINT nFlags, CPoint point)
{// Нажатие средней кнопки
	CBasicView_Parent::OnMButtonDown(nFlags, point);// Базовая реализация
	if(!m_bLButtonPressed) SetCapture();
	m_bMButtonPressed = true;
}
////////////////////////////////////////////////////////////////////////////////
void CBasicView::OnMButtonUp(UNINT nFlags, CPoint point)
{// отпускание средней кнопки
	CBasicView_Parent::OnMButtonUp(nFlags, point);// Базовая реализация
	if(!m_bLButtonPressed && !m_bRButtonPressed) ReleaseCapture();
	m_bMButtonPressed = false;
}
////////////////////////////////////////////////////////////////////////////////
void CBasicView::OnMouseMove(UNINT nFlags, CPoint point)
{// движение мыши
	CPoint Shift;
	//
	CBasicView_Parent::OnMouseMove(nFlags, point);// Базовая реализация
	if(!m_bMouseIsInside){
		m_bMouseIsInside = true;
		m_MousePos = point;
		SetFocus();
	}// если мышь была за пределами окна, устанавливаем флаг, что она вернулась
	Shift=m_MousePos-point;//Смещение мыши
	m_MousePos = point;// Позиция мыши
//	TrackMouseEvent(TME_LEAVE);// мышь внутри окна, ждем когда уйдет
}
////////////////////////////////////////////////////////////////////////////////
int CBasicView::OnMouseWheel(UNINT nFlags, short zDelta, CPoint pt)
{// прокрутка колесика
	return CBasicView_Parent::OnMouseWheel(nFlags, zDelta, pt);
}
////////////////////////////////////////////////////////////////////////////////
void CBasicView::OnMouseLeave()
{// Мышь покинула область окна
	CBasicView_Parent::OnMouseLeave();
	m_bMouseIsInside = false;
}

////////////////////////////////////////////////////////////////////////////////
void CBasicView::TrackMouseEvent(UNLONG dwFlags)
{// Вспомогательная функция, заставляет объект обрабатывать уход мыши за пределы окна
	TRACKMOUSEEVENT     csTME;
	csTME.cbSize = sizeof(csTME);
	csTME.dwFlags = dwFlags;
	csTME.hwndTrack = GetSafeHwnd();
	::_TrackMouseEvent(&csTME);
}

////////////////////////////////////////////////////////////////////////////////
int CBasicView::OnEraseBkgnd(CDC* pDC){
	return 1;
}

void CBasicView::OnDestroy(){
	if(m_pDoc){
		m_pDoc->SetView(nullptr);
		m_pDoc = nullptr;
	}
	CWnd::OnDestroy();
}

void CBasicView::DrawDoc(GLParam &parameter){
	if (m_pDoc){
		m_pDoc->DrawGL(parameter);
	}
}
