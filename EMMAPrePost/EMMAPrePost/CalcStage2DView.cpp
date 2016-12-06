#include "StdAfx.h"
#include "CalcStage2DView.h"

IMPLEMENT_DYNCREATE(CCalcStage2DView, CGL2DView);

BEGIN_MESSAGE_MAP(CCalcStage2DView, CGL2DView)
	ON_REGISTERED_MESSAGE(WMR_DRAFTRECT, OnDraftRect)
	ON_WM_SIZE()
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_SETCURSOR()
	ON_WM_LBUTTONUP()
	ON_COMMAND(ID_3DVIEW, &CCalcStage2DView::On3D)
	ON_REGISTERED_MESSAGE(WMR_ACTIVATE, OnActivate)
END_MESSAGE_MAP()

CCalcStage2DView::CCalcStage2DView(void)
{
	m_CurrentCursor = 0;
	m_ResizeState = 0;
	m_IsScalePicked = false;
	m_IsScaleResizeMode = false;
	
	m_Parameter = new GLParam();
}

CCalcStage2DView::~CCalcStage2DView(void)
{
	if (m_Parameter){
		delete m_Parameter;
		m_Parameter = nullptr;
	}
}

void CCalcStage2DView::SetDocument(CEMMADoc *pDoc)
{
	m_Parameter->m_PolygonBuilder.ClearDisplayList();
	m_Parameter->m_ScaleVisualizer.ClearDisplayList();
	m_Parameter->m_ScaleVisualizer.SetFont(&m_Font);

	CChildFrame *pChildFrame = dynamic_cast<CChildFrame *>(GetParent());
	CMDIFrameWnd* pFrame = pChildFrame->GetMDIFrame();
	pFrame->SendMessage(WMR_UPDATEPANE);	//TODO: Error & Crash после повторной загрузки результатов


	CBasicView::SetDocument(pDoc);
}

void CCalcStage2DView::UpdateResizeState(const CPoint &point)
{
	int result = 0;
	//CRect scaleRect = m_Parameter.m_ScaleVisualizer.m_ScaleRect;
	CRect scaleRect = m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ToCRect();

	if ((-10 < point.x - scaleRect.left) && (point.x - scaleRect.left < 5) &&
		(scaleRect.top < point.y) && (point.y < scaleRect.bottom)){
		result |= RESIZE_LEFT;
	}

	if ((-5 < point.x - scaleRect.right) && (point.x - scaleRect.right < 10) &&
		(scaleRect.top < point.y) && (point.y < scaleRect.bottom)){
		result |= RESIZE_RIGHT;
	}

	if ((-10 < point.y - scaleRect.top) && (point.y - scaleRect.top < 5) &&
		(scaleRect.left < point.x) && (point.x < scaleRect.right)){
		result |= RESIZE_TOP;
	}

	if ((-5 < point.y - scaleRect.bottom) && (point.y - scaleRect.bottom < 10) &&
		(scaleRect.left < point.x) && (point.x < scaleRect.right)){
		result |= RESIZE_BOTTOM;
	}

	m_ResizeState = result;
}

bool CCalcStage2DView::IsInScale(const CPoint &point)
{
	//CRect scaleRect = m_Parameter.m_ScaleVisualizer.m_ScaleRect;
	CRect scaleRect = m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ToCRect();
	return	(point.x < scaleRect.right && point.y < scaleRect.bottom) &&
			(point.x > scaleRect.left && point.y > scaleRect.top);
}

void CCalcStage2DView::SetCursor(const CPoint &point)
{	
	if (m_nCurType != 0 || !IsInScale(point)){
		m_CurrentCursor = 0;
	}
	else{
		UpdateResizeState(point);
		
		switch(m_ResizeState)
		{
			case RESIZE_LEFT:
			case RESIZE_RIGHT:
				::SetCursor(AfxGetApp()->LoadStandardCursor(IDC_SIZEWE));
				m_CurrentCursor = (ptrdiff_t)IDC_SIZEWE;
				break;

			case RESIZE_TOP:
			case RESIZE_BOTTOM:
				::SetCursor(AfxGetApp()->LoadStandardCursor(IDC_SIZENS));
				m_CurrentCursor = (ptrdiff_t)IDC_SIZENS;
				break;

			case RESIZE_LEFT | RESIZE_BOTTOM:
			case RESIZE_RIGHT | RESIZE_TOP:
				::SetCursor(AfxGetApp()->LoadStandardCursor(IDC_SIZENESW));
				m_CurrentCursor = (ptrdiff_t)IDC_SIZENESW;
				break;

			case RESIZE_LEFT | RESIZE_TOP:
			case RESIZE_RIGHT | RESIZE_BOTTOM:
				::SetCursor(AfxGetApp()->LoadStandardCursor(IDC_SIZENWSE));
				m_CurrentCursor = (ptrdiff_t)IDC_SIZENWSE;
				break;

			default:
				::SetCursor(AfxGetApp()->LoadCursorW(IDC_HAND_O));
				m_CurrentCursor = IDC_HAND_O;
				break;
		}
	}	
}

void CCalcStage2DView::MoveScale(const CPoint &point)
{
	CPoint diff = (point - m_pPrevPos);

	m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ShiftBy(diff);

	m_pPrevPos = point;
	m_Parameter->m_ScaleVisualizer.m_IsMoved = true;
	Invalidate();		
}

void CCalcStage2DView::ResizeScale(const CPoint &point)
{
	CPoint diff = (point - m_pPrevPos);

	switch(m_ResizeState)
	{
		case RESIZE_LEFT:
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeLeftBottom(diff.x, 0);
			break;

		case RESIZE_RIGHT:
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeRightTop(diff.x, 0);
			break;

		case RESIZE_TOP:
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeRightTop(0, diff.y);
			break;

		case RESIZE_BOTTOM:
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeLeftBottom(0, diff.y);
			break;

		case RESIZE_LEFT | RESIZE_BOTTOM:
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeLeftBottom(diff);
			break;

		case RESIZE_RIGHT | RESIZE_TOP:
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeRightTop(diff);

			break;

		case RESIZE_LEFT | RESIZE_TOP:
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeLeftBottom(diff.x, 0);
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeRightTop(0, diff.y);
			break;

		case RESIZE_RIGHT | RESIZE_BOTTOM:
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeRightTop(diff.x, 0);
			m_Parameter->m_ScaleVisualizer.m_ScaleViewport.ResizeLeftBottom(0, diff.y);
			break;
	}

	m_pPrevPos = point;
	m_Parameter->m_ScaleVisualizer.m_IsMoved = true;
	Invalidate();
}

LRESULT CCalcStage2DView::OnDraftRect(WPARAM wParam, LPARAM lParam){
	CRect2D *rect = (CRect2D *)wParam;
	SetDraftRect(*rect, lParam != 0);
	SetProjection();
	Invalidate();
	UpdateWindow();	

	return 0;
}

void CCalcStage2DView::OnSize(UNINT nType, int cx, int cy)
{
	CGL2DView::OnSize(nType, cx, cy);

	// Если отображение шкалы никто не двигал и не менял ёё размеры, то устанавливаем default значения:
	// Высота = 1/2 высоты окна вывода,
	// Ширина = 1/3 ширина окна,
	// Отступ сверху = 5% высоты экрана,
	// Отступ справа = 5% ширины экрана
	if (!m_Parameter->m_ScaleVisualizer.m_IsMoved)
	{
		CRect scaleRect, clientRect;
		GetClientRect(clientRect);
		
		scaleRect.bottom = clientRect.top + static_cast<long>(clientRect.Height() / 2.0);
		scaleRect.left = clientRect.right - static_cast<long>(clientRect.Width() / 3.0);
		scaleRect.top = clientRect.top + static_cast<long>(clientRect.Height() * 0.05);
		scaleRect.right = clientRect.right - static_cast<long>(clientRect.Width() * 0.05);

		m_Parameter->m_ScaleVisualizer.m_ScaleViewport = GLViewport(scaleRect);
	}
}


void CCalcStage2DView::OnMouseMove(UNINT nFlags, CPoint point)
{
	if (m_Parameter->m_ScaleVisualizer.IsScaleVisible()){	
		if (!m_IsScalePicked && !m_IsScaleResizeMode){
			SetCursor(point);
		}
		else if (!m_IsScaleResizeMode){
			MoveScale(point);
		}
		else{
			ResizeScale(point);
		}
	}

	CGL2DView::OnMouseMove(nFlags, point);
}


void CCalcStage2DView::OnLButtonDown(UNINT nFlags, CPoint point)
{	
	m_pPrevPos = point;

	if (m_CurrentCursor == IDC_HAND_O)
	{
		m_IsScalePicked = true;
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_HAND_P));
		m_CurrentCursor = IDC_HAND_P;
	}
	else if (m_CurrentCursor != 0)
	{
		m_IsScaleResizeMode = true;
	}

	CGL2DView::OnLButtonDown(nFlags, point);
}


int CCalcStage2DView::OnSetCursor(CWnd* pWnd, UNINT nHitTest, UNINT message)
{
	if(m_nCurType == 0 && m_CurrentCursor != 0){		
		return 1;
	}

	return CGL2DView::OnSetCursor(pWnd, nHitTest, message);
}


void CCalcStage2DView::OnLButtonUp(UNINT nFlags, CPoint point)
{
	m_IsScalePicked = false;
	m_IsScaleResizeMode = false;

	if (m_CurrentCursor == IDC_HAND_P)
	{
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_HAND_O));
		m_CurrentCursor = IDC_HAND_O;
	}

	CGL2DView::OnLButtonUp(nFlags, point);
}


void CCalcStage2DView::On3D() {
	CCalcStage *pDoc = dynamic_cast<CCalcStage *>(GetDocument());
	if (!pDoc) {
		return;
	}
	bool bFlagResult = pDoc->Set3DFlag(true);
	//GetParentFrame()->SendMessage(WMR_EP_SELCHANGE, (WPARAM)pDoc);
	GetTopLevelFrame()->SendMessage(WMR_EP_SELCHANGE, (WPARAM)pDoc);
}

LRESULT CCalcStage2DView::OnActivate(WPARAM wParam, LPARAM lParam) {
	if (wParam) {
		CCalcStage *pDoc = dynamic_cast<CCalcStage *>(GetDocument());
		if (!pDoc)return 0;
		pDoc->Set3DFlag(false);
	}
	return 0;
}