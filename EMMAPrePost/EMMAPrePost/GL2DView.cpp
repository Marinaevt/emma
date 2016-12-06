#include "StdAfx.h"
#include "GL2DView.h"
#include "Results.h"

BEGIN_MESSAGE_MAP(CGL2DView, CGLBasicView)
	ON_WM_MOUSEWHEEL()
	ON_WM_SETCURSOR()
	ON_WM_LBUTTONDOWN()
	ON_WM_RBUTTONDOWN()
	ON_WM_MOUSEMOVE()

	ON_REGISTERED_MESSAGE(WMR_INITIALUPDATE, OnInitialUpdate)

	ON_COMMAND(ID_ARROW, &CGL2DView::OnArrow)
	ON_UPDATE_COMMAND_UI(ID_ARROW, &CGL2DView::OnUpdateArrow)
	ON_COMMAND(ID_BUTTON_ZOOM, &CGL2DView::OnZoom)
	ON_UPDATE_COMMAND_UI(ID_BUTTON_ZOOM, &CGL2DView::OnUpdateZoom)
	ON_COMMAND(ID_SHOWALL, &CGL2DView::OnShowAll)
	ON_UPDATE_COMMAND_UI(ID_SHOWALL, &CGL2DView::OnUpdateShowAll)
	ON_COMMAND(ID_HAND, &CGL2DView::OnHand)
	ON_UPDATE_COMMAND_UI(ID_HAND, &CGL2DView::OnUpdateHand)
	
	ON_COMMAND(ID_AXIS, &CGL2DView::OnAxis)
	ON_UPDATE_COMMAND_UI(ID_AXIS, &CGL2DView::OnUpdateAxis)
	ON_COMMAND(ID_GRID, &CGL2DView::OnGrid)
	ON_UPDATE_COMMAND_UI(ID_GRID, &CGL2DView::OnUpdateGrid)
	ON_COMMAND(ID_COLOR, &CGL2DView::OnColor)
	ON_UPDATE_COMMAND_UI(ID_COLOR, &CGL2DView::OnUpdateColor)

	ON_COMMAND(ID_3DVIEW, &CGL2DView::On3DView)
	ON_UPDATE_COMMAND_UI(ID_3DVIEW, &CGL2DView::OnUpdate3DView)

END_MESSAGE_MAP()

///////////////////////////////////////////////////////////////////
CGL2DView::CGL2DView(void)
	:m_DraftRect(-1, 1, -1, 1),
	m_dCx(0), m_dCy(0), m_dZoom(1),
	m_bShowAxis(1), m_bShowGrid(1)
{ 
	{
	double clrGrid1[4] = { 0., 0., .5, .2 };
	double clrGrid5[4] = { 0., 0., 0., .4 };
	double clrGrid10[4] = { 0., 0., 0.2, .6 };

	double clrRuler1[4] = { .4, .4, .6, 1. };
	double clrRuler5[4] = { .2, .2, .4, 1. };
	double clrRuler10[4] = { .0, .0, .2, 1. };
	double clrRulerBG[4] = { .8, .8, 1., .7 };
	double clrRulerB[4] = { .3, .3, .5, 1. };

	m_AxisParam.Stroke1.nType = 1;
	m_AxisParam.Stroke5.nType = 5;
	m_AxisParam.Stroke10.nType = 10;
	m_AxisParam.Stroke50.nType = 50;
	m_AxisParam.Stroke100.nType = 100;

	for (int i = 0; i < 4; i++) {
		m_AxisParam.Stroke1.ColorG[i] = clrGrid1[i];
		m_AxisParam.Stroke5.ColorG[i] = clrGrid5[i];
		m_AxisParam.Stroke10.ColorG[i] = clrGrid10[i];
		m_AxisParam.Stroke50.ColorG[i] = clrGrid10[i];
		m_AxisParam.Stroke100.ColorG[i] = clrGrid10[i];

		m_AxisParam.Stroke1.ColorR[i] = clrRuler1[i];
		m_AxisParam.Stroke5.ColorR[i] = clrRuler5[i];
		m_AxisParam.Stroke10.ColorR[i] = clrRuler10[i];
		m_AxisParam.Stroke50.ColorR[i] = clrRuler10[i];
		m_AxisParam.Stroke100.ColorR[i] = clrRuler10[i];

		m_AxisParam.ColorRuler[i] = clrRulerBG[i];
		m_AxisParam.ColorRulerB[i] = clrRulerB[i];
	}
	m_AxisParam.nScaleWidth = 20;
	m_AxisParam.Stroke1.nLehgth = m_AxisParam.nScaleWidth / 4;
	m_AxisParam.Stroke5.nLehgth = m_AxisParam.nScaleWidth / 3;
	m_AxisParam.Stroke10.nLehgth = m_AxisParam.nScaleWidth / 2;
	m_AxisParam.Stroke50.nLehgth = m_AxisParam.nScaleWidth / 2;
	m_AxisParam.Stroke100.nLehgth = m_AxisParam.nScaleWidth / 2;
}
    Is3DViewShow = false;
    View = false;
}

CGL2DView::~CGL2DView(void){

}

///////////////////////////////////////////////////////////////////
void CGL2DView::SetProjection(void){
	CRect2D rectScreen;
	GetScreenRect(rectScreen);

	if(SetCurrentContext()){
		SetViewport();
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(rectScreen.left, 
				rectScreen.right,
				rectScreen.bottom,
				rectScreen.top,
				-1, 1);
		glMatrixMode(GL_MODELVIEW);
		::wglMakeCurrent(nullptr, nullptr);
	}

	RecalcAxisParam();
/*	SWPARAM	Param ={m_dMinX, m_dMaxX, m_dMinY, m_dMaxY, 
					pDoc->GetMinX() - m_dWidth*m_dKff*0.0499,  pDoc->GetMaxX() + m_dWidth*m_dKff*0.0499, 
					pDoc->GetMinY() - m_dHeight*m_dKff*0.0499, pDoc->GetMaxY() + m_dHeight*m_dKff*0.0499};
	
	GetParent()->SendMessage(WM_SW_SCROLL, (WPARAM)&Param);*/
}

///////////////////////////////////////////////////////////////////
void CGL2DView::DrawScene(void){
	glEnable(GL_ALPHA_TEST);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	if (m_bShowGrid) {
		DrawGrid();
	}
	DrawDoc(*m_Parameter);
	if (m_bShowAxis) {
		DrawAxis();
	}
}

///////////////////////////////////////////////////////////////////
bool CGL2DView::GetScreenRect(CRect2D &rect){
	CRect clientRect;
	GetClientRect(&clientRect);

	if(clientRect.Width()==0 || clientRect.Height()==0)return false;
	
	double dWidth, dHeigth;
	if(m_DraftRect.Width()/m_DraftRect.Heigth() > static_cast<double>(clientRect.Width())/static_cast<double>(clientRect.Height())){
		dWidth = m_DraftRect.Width()/m_dZoom;
		dHeigth = m_DraftRect.Width()/static_cast<double>(clientRect.Width())*static_cast<double>(clientRect.Height())/m_dZoom;
	}else{
		dHeigth = m_DraftRect.Heigth()/m_dZoom;
		dWidth = m_DraftRect.Heigth()/static_cast<double>(clientRect.Height())*static_cast<double>(clientRect.Width())/m_dZoom;
	}

	dHeigth = abs(dHeigth);
	dWidth = abs(dWidth);

	rect.left   = m_dCx - dWidth*0.5;
	rect.right  = m_dCx + dWidth*0.5;
	rect.bottom = m_dCy - dHeigth*0.5;
	rect.top    = m_dCy + dHeigth*0.5;

	return true;
}

///////////////////////////////////////////////////////////////////
void CGL2DView::DrawGrid(){
	DrawRuling(XY_axis, ST_grid);
}

///////////////////////////////////////////////////////////////////
void CGL2DView::DrawAxis(){
	glLineWidth(1);
	CRect2D ScreenRect;
	GetScreenRect(ScreenRect);
	double dW = m_AxisParam.dPCX*m_AxisParam.nScaleWidth;
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColor4dv(m_AxisParam.ColorRuler);
	glBegin(GL_QUADS);
	glVertex2d(ScreenRect.left, ScreenRect.bottom);
	glVertex2d(ScreenRect.left+dW, ScreenRect.bottom);
	glVertex2d(ScreenRect.left+dW, ScreenRect.top);
	glVertex2d(ScreenRect.left, ScreenRect.top);
	glVertex2d(ScreenRect.left, ScreenRect.bottom);
	glVertex2d(ScreenRect.right, ScreenRect.bottom);
	glVertex2d(ScreenRect.right, ScreenRect.bottom+dW);
	glVertex2d(ScreenRect.left, ScreenRect.bottom+dW);
	glEnd();

	glColor4dv(m_AxisParam.ColorRulerB);
	glBegin(GL_LINE_STRIP);
	glVertex2d(ScreenRect.left+dW, ScreenRect.top);
	glVertex2d(ScreenRect.left+dW, ScreenRect.bottom+dW);
	glVertex2d(ScreenRect.right,   ScreenRect.bottom+dW);
	glEnd();
	DrawRuling(XY_axis, ST_axis);
	glColor4d(0,0,0,1);
}

///////////////////////////////////////////////////////////////////
void CGL2DView::GetStrokeParam(size_t n, STROKE_PARAM *pStrokeParam){
	if(n%100 == 0){
		if(pStrokeParam) *pStrokeParam = m_AxisParam.Stroke100;
	}else if(n%50 == 0){
		if(pStrokeParam) *pStrokeParam = m_AxisParam.Stroke50;
	}else if(n%10 == 0){
		if(pStrokeParam) *pStrokeParam = m_AxisParam.Stroke10;
	}else if(n%5 == 0){
		if(pStrokeParam) *pStrokeParam = m_AxisParam.Stroke5;
	}else{
		if(pStrokeParam) *pStrokeParam = m_AxisParam.Stroke1;
	}
}

///////////////////////////////////////////////////////////////////
void CGL2DView::DrawRuling(size_t nDirection, size_t nStrokeType){
	glLineWidth(1);
	if(nDirection == XY_axis){
		DrawRuling(X_axis, nStrokeType);
		DrawRuling(Y_axis, nStrokeType);
		return;
	}

	CRect2D rectScreen;
	GetScreenRect(rectScreen);
	double dMin = (nDirection == X_axis) ? rectScreen.left : rectScreen.bottom;
	double dMax = (nDirection == X_axis) ? rectScreen.right : rectScreen.top;
	double XY0 = floor(dMin/m_AxisParam.dDelta/100)*100*m_AxisParam.dDelta;
	double dCurrXY = XY0;
	size_t i = 0;
	
	STROKE_PARAM Stroke;
	//bool	bDrawText = 1; не используется
	double	dTextPos = dCurrXY;


	while(dCurrXY < dMax){
		GetStrokeParam(i, &Stroke);

		DrawStroke(Stroke, dCurrXY, nDirection, nStrokeType);
		
		if(dTextPos <= dCurrXY && nStrokeType == ST_axis){
			double dDrawnTextPos = DrawText(Stroke, dCurrXY, nDirection)*1.3 + dCurrXY;
			dTextPos = max(dDrawnTextPos, dTextPos);
		}
		i++;
		dCurrXY = static_cast<double>(i)*m_AxisParam.dDelta + XY0;
	}
}

///////////////////////////////////////////////////////////////////
double CGL2DView::DrawText(const STROKE_PARAM &StrokeParam, double dPos, int nDirection){
	//if(StrokeParam.nType == 1)return 0;
	double H = 0.5*m_AxisParam.dPCX*m_AxisParam.nScaleWidth;

	CRect2D ScreenRect;
	GetScreenRect(ScreenRect);

	//dPos*=m_dZoom;
	
	CString str;
	str.Format(_T("%f"), dPos);
	//if(str.Find(static_cast<TCHAR>('.'))>=0){
	if (str.Find(_T('.')) >= 0) {
		//while(str.GetAt(str.GetLength()-1) == static_cast<_TCHAR>('0')){
		while (str.GetAt(str.GetLength() - 1) == _T('0')) {
			str = str.Left(str.GetLength()-1);
		}
	}
	double W = m_Font.GetTextWidth(str.GetBuffer(), H);

	if (StrokeParam.nType == 1 && W*1.5 > m_AxisParam.dDelta) {
		return 0;
	}
	if (StrokeParam.nType == 5 && W*1.5 + dPos > m_AxisParam.dDelta * 10 * (floor(dPos / m_AxisParam.dDelta / 10) + 1)) {
		return 0;
	}
	
	double X, Y, alpha;
	if(nDirection == X_axis){
		if (dPos <= ScreenRect.left + m_AxisParam.dPCX*m_AxisParam.nScaleWidth) {
			return 0;
		}
		Y = ScreenRect.bottom + 0.5*(m_AxisParam.dPCX*m_AxisParam.nScaleWidth - H);
		X = dPos;
		alpha = 0;
	}else if(nDirection == Y_axis){
		if (dPos <= ScreenRect.bottom + m_AxisParam.dPCX*m_AxisParam.nScaleWidth) {
			return 0;
		}
		X = ScreenRect.left + 0.5*(m_AxisParam.dPCX*m_AxisParam.nScaleWidth + H);
		Y = dPos;
		alpha = 90;
	}else{
		CDlgShowError cError(ID_ERROR_GL2DVIEW_DIRECTION_WRONG); //_T("Направление неверно nDirection"));
		return 0;
	}

	glColor4d(0,0,0,1);
	m_Font.OutText2D(str.GetBuffer(), X, Y, H, alpha);
	return W;
}

///////////////////////////////////////////////////////////////////
void CGL2DView::DrawStroke(const STROKE_PARAM &StrokeParam, double dPos, int nDirection, int nStrokeType){
	CRect2D ScreenRect;
	GetScreenRect(ScreenRect);

	glBegin(GL_LINES);
	if(nStrokeType == ST_grid){
		glColor4dv(StrokeParam.ColorG);
		if(nDirection == X_axis){
			glVertex2d(dPos, ScreenRect.bottom);
			glVertex2d(dPos, ScreenRect.top);
		}else{
			glVertex2d(ScreenRect.left,  dPos);
			glVertex2d(ScreenRect.right, dPos);
		}
	}
	if(nStrokeType == ST_axis){
		glColor4dv(StrokeParam.ColorR);
		double dA = m_AxisParam.dPCX*m_AxisParam.nScaleWidth;
		double dB = m_AxisParam.dPCX*(m_AxisParam.nScaleWidth-StrokeParam.nLehgth);
		if(nDirection == X_axis && dPos > ScreenRect.left + dA){
			glVertex2d(dPos, ScreenRect.bottom + dA);
			glVertex2d(dPos, ScreenRect.bottom + dB);
		}
		if(nDirection == Y_axis && dPos > ScreenRect.bottom + dA){
			glVertex2d(ScreenRect.left + dA, dPos);
			glVertex2d(ScreenRect.left + dB, dPos);
		}
	}
	glEnd();
	glColor4d(0,0,0,1);
}

///////////////////////////////////////////////////////////////////
void CGL2DView::RecalcAxisParam(){
	CRect ClientRect;
	CRect2D ScreenRect;
	GetClientRect(&ClientRect);
	GetScreenRect(ScreenRect);

	m_AxisParam.dPCX = ScreenRect.Width()/ClientRect.Width();
	m_AxisParam.dDelta = pow(10, floor(log10(ScreenRect.Width()))-1);
}

// Показать оси
void CGL2DView::ShowAxis(bool bShow, bool bUpdate){
	m_bShowAxis = bShow;
	if(bUpdate){
		Invalidate(0);
		UpdateWindow();
	}
}

// Показать сетку
void CGL2DView::ShowGrid(bool bShow, bool bUpdate){
	m_bShowGrid = bShow;
	if(bUpdate){
		Invalidate(0);
		UpdateWindow();
	}
}

// Смена типа курсора
void CGL2DView::ChangeCursor(int curType, bool bUpdate){
	m_nCurType=curType;
	if(bUpdate){
		Invalidate(0);
		UpdateWindow();
	}
}

// Установка масштаба
void CGL2DView::SetZoom(double dZoom, bool bUpdate){
	if(dZoom <= 0.0){ 
		CDlgShowError cError(ID_ERROR_GL2DVIEW_ZOOM_WRONG); //_T("Масштаб меньше 0"));
		return; 
	}

	m_dZoom = dZoom;
	
	if(bUpdate){
		SetProjection();
		Invalidate(0);
		UpdateWindow();
	}	
}

// Установка центра области чертежа
void CGL2DView::SetCentr(double x, double y){
	
	CRect2D scrRect; 
	GetScreenRect(scrRect);

	if(m_DraftRect.Width() <= scrRect.Width()){
		m_dCx = min(max(x, m_DraftRect.left - scrRect.Width()*0.5), m_DraftRect.right + scrRect.Width()*0.5);
	}else{
		m_dCx = max(min(x, m_DraftRect.right + scrRect.Width()*0.5), m_DraftRect.left - scrRect.Width()*0.5);
	}

	if(m_DraftRect.Heigth() <= scrRect.Heigth()){
		m_dCy = min(max(y, m_DraftRect.bottom - scrRect.Heigth()*0.5), m_DraftRect.top + scrRect.Heigth()*0.5);
	}else{
		m_dCy = max(min(y, m_DraftRect.top + scrRect.Heigth()*0.5), m_DraftRect.bottom - scrRect.Heigth()*0.5);
	}

}

int CGL2DView::OnMouseWheel(UNINT nFlags, short zDelta, CPoint pt)
{
	// TODO: Add your message handler code here and/or call default

	Math::C2DPoint asd;
	asd=FizToLogical(m_MousePos);

	if(zDelta<0) 
	{
		SetZoom(GetZoom()*1.1,false);
		SetCentr(asd.x+(m_dCx-asd.x)/1.1,asd.y+(m_dCy-asd.y)/1.1);
	}
	else 
	{
		SetZoom(GetZoom()/1.1,false);
		SetCentr(asd.x+(m_dCx-asd.x)*1.1,asd.y+(m_dCy-asd.y)*1.1);
	}
	
	SetProjection();
	Invalidate(0);

	return CGLBasicView::OnMouseWheel(nFlags, zDelta, pt);
}

// Установка вида курсора
int CGL2DView::OnSetCursor(CWnd* pWnd, UNINT nHitTest, UNINT message)
{
	//Курсор-лупа
	if(m_nCurType==1){
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_ZOOM));
		return 1;
	}

	//Курсор-рука
	if(m_nCurType==2){
		if(m_bLButtonPressed)	::SetCursor(AfxGetApp()->LoadCursorW(IDC_HAND_P));
		else					::SetCursor(AfxGetApp()->LoadCursorW(IDC_HAND_O));
		return 1;
	}

	//Курсор-добавление узла
	if (m_nCurType == 100) {
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_STICK));
		return 1;
	}

	//Курсор-привязка
	if(m_nCurType == 200){
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_STICK));
		return 1;
	}
	
	//Курсор-КЭ сгущение 
	if(m_nCurType == 300){
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_STICK));
		return 1;
	}

	return CGLBasicView::OnSetCursor(pWnd, nHitTest, message);
}

// Действие по левому клику на области чертежа
void CGL2DView::OnLButtonDown(UNINT nFlags, CPoint point)
{
	//Курсор-лупа
	if(m_nCurType==1){
		Math::C2DPoint asd=FizToLogical(point);
		SetZoom(GetZoom()*2,false);
		SetCentr(asd.x,asd.y);
		SetProjection();
		Invalidate(0);
	}
	
	//Курсор-рука
	if(m_nCurType==2){
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_HAND_P));
		m_pPrevPos=point; 
	}
	
	CGLBasicView::OnLButtonDown(nFlags, point);
}

// Действие по правому клику на области чертежа
void CGL2DView::OnRButtonDown(UNINT nFlags, CPoint point)
{
	if(m_nCurType==1){
		Math::C2DPoint asd=FizToLogical(point);
		SetZoom(GetZoom()*0.5, false);
		SetCentr(asd.x,asd.y);
		SetProjection();
		Invalidate(0);
	}

	CGLBasicView::OnRButtonDown(nFlags, point);
}


Math::C2DPoint CGL2DView::FizToLogical(CPoint FizPoint)
{
	Math::C2DPoint Log;
	CRect clientRect;
	GetClientRect(&clientRect);

	Log.x=m_dCx+(FizPoint.x - static_cast<double>(clientRect.right)/2.0)*m_AxisParam.dPCX;
	Log.y=m_dCy+(static_cast<double>(clientRect.bottom)/2.0 - FizPoint.y)*m_AxisParam.dPCX;

	return Log;
}

void CGL2DView::SetDraftRect(const CRect2D &rect, bool bCenter){
	m_DraftRect = rect;
	if(bCenter){
		m_dCx = (rect.right + rect.left)*0.5;
		m_dCy = (rect.top + rect.bottom)*0.5;
		SetCentr(m_dCx, m_dCy);
	}
}


void CGL2DView::OnMouseMove(UNINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	if(m_nCurType==2 && m_bLButtonPressed)
	{
		Math::C2DPoint asd1,asd2;
		asd1=FizToLogical(m_pPrevPos);
		asd2=FizToLogical(point);
		SetCentr(m_dCx+asd1.x-asd2.x,m_dCy+asd1.y-asd2.y);
		m_pPrevPos=point;
		SetProjection();
		Invalidate(0);
	}

	CGLBasicView::OnMouseMove(nFlags, point);
}

// Действия по клику на иконке туллбара
void CGL2DView::OnArrow(){
	// TODO: Add your command handler code here
	ChangeCursor(0,true);
}
void CGL2DView::OnUpdateArrow(CCmdUI *pCmdUI){
	// TODO: Add your command update UI handler code here
	if (GetCursorType() == 0) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void CGL2DView::OnZoom(){
	// TODO: Add your command handler code here
	if(GetCursorType()==1) ChangeCursor(0,true);
	else ChangeCursor(1,true);

	/*
	Invalidate(0);
	SetProjection();
	UpdateWindow();//*/
}
void CGL2DView::OnUpdateZoom(CCmdUI *pCmdUI){
	
	if(GetCursorType()==1) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void CGL2DView::OnShowAll(){
	
	CRect2D rect;

	if (m_pDoc->GetBoundingBox(rect)){
		SetDraftRect(rect, true);
	}else{
		SetDraftRect(CRect2D(-1.0,1.0,-1.0,1.0), true);	//заглушка, если ничего не загружено
	}
	SetZoom(1/1.1);
}

void CGL2DView::OnUpdateShowAll(CCmdUI *pCmdUI){
}

void CGL2DView::OnHand(){
	// TODO: Add your command handler code here
	if(GetCursorType()==2) ChangeCursor(0,true);
	else ChangeCursor(2,true);

	/* 
	Invalidate(0);
	SetProjection();
	UpdateWindow();//*/
}
void CGL2DView::OnUpdateHand(CCmdUI *pCmdUI){
	// TODO: Add your command update UI handler code here
	if(GetCursorType()==2) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);
}

void CGL2DView::OnAxis(){
	ShowAxis(!IsAxisVisible());
}
void CGL2DView::OnUpdateAxis(CCmdUI *pCmdUI){
	pCmdUI->SetCheck(IsAxisVisible());
}

void CGL2DView::OnGrid(){
	ShowGrid(!IsGridVisible());
}
void CGL2DView::OnUpdateGrid(CCmdUI *pCmdUI){
	pCmdUI->SetCheck(IsGridVisible());
}

void CGL2DView::OnColor(){
	// TODO: Add your command handler code here
}
void CGL2DView::OnUpdateColor(CCmdUI *pCmdUI){
	// TODO: Add your command update UI handler code here
}


LRESULT CGL2DView::OnInitialUpdate(WPARAM wParam, LPARAM lParam){
	if(!m_pDoc)return 0;
	SetZoom(1/1.1);
	CRect2D rect;
	if(m_pDoc->GetBoundingBox(rect)) { 
		SendMessage(WMR_DRAFTRECT, (WPARAM)&rect, (LPARAM)1);
	}else{
		//Если документ ещё не загружен, то ставим вид по-умолчанию
		rect = CRect2D(-1.0,1.0,-1.0,1.0);
		SendMessage(WMR_DRAFTRECT, (WPARAM)&rect, (LPARAM)1); 
	};
	return 0;
}

void CGL2DView::On3DView() {

	Is3DViewShow = !Is3DViewShow;

}

void CGL2DView::OnUpdate3DView(CCmdUI *pCmdUI) {

	pCmdUI->SetCheck(Is3DViewShow);
	if (Is3DViewShow) {


	}
	else
	{


	}

}