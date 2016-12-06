// ChildView.cpp : implementation of the CChild2DView class
//

#include "stdafx.h"
#include "Child2DView.h"
#include "2DMesher.h"
#include "2DSketch.h"
#include "2DTool.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CChild2DView

CChild2DView::CChild2DView(){
	m_pI2DSketch = nullptr;
}

CChild2DView::~CChild2DView(){
}


IMPLEMENT_DYNCREATE(CChild2DView, CGL2DView);


BEGIN_MESSAGE_MAP(CChild2DView, CGL2DView)
	ON_REGISTERED_MESSAGE(WMR_DRAFTRECT, OnDraftRect)
	ON_COMMAND(ID_IMPORT_PARAMETRIC, &CChild2DView::OnImportPrm)
	ON_UPDATE_COMMAND_UI(ID_IMPORT_PARAMETRIC, &CChild2DView::OnUpdateImportPrm)
	ON_COMMAND(ID_IMPORT_DXF, &CChild2DView::OnImportDxf)
	ON_UPDATE_COMMAND_UI(ID_IMPORT_DXF, &CChild2DView::OnUpdateImportDxf)

	ON_REGISTERED_MESSAGE(WMR_SETCURSOR, &CChild2DView::OnSetCursor)
	
	

	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONUP()
	ON_WM_KEYDOWN()
	ON_WM_KEYUP()

	ON_WM_TIMER()		// реагировать на таймер

END_MESSAGE_MAP()

// CChild2DView message handlers

int CChild2DView::PreCreateWindow(CREATESTRUCT& cs){
	if (!CGL2DView::PreCreateWindow(cs))
		return 0;

	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(nullptr, IDC_ARROW), reinterpret_cast<HBRUSH>(COLOR_WINDOW+1), nullptr);

	return 1;
}

////////////////////////////////////////////////////////////////////////////////
void CChild2DView::DrawScene(void){
	CGL2DView::DrawScene();
}

////////////////////////////////////////////////////////////////////////////////
LRESULT CChild2DView::OnDraftRect(WPARAM wParam, LPARAM lParam){
	CRect2D *rect = (CRect2D *)wParam;
	SetDraftRect(*rect, lParam != 0);
	SetProjection();
	Invalidate();
	UpdateWindow();
	return 0;
}


////////////////////////////////////////////////////////////////////////////////
void CChild2DView::OnImportPrm(){
	
	// Запрещаем КЭ сетке загружать чертёж
	if (m_pDoc) {
		if (m_pDoc->GetType() == C2DMESHER) {
			return;
		}
	}

	CDlgParametricPreview dlgParametricPreView;
	if(dlgParametricPreView.DoModal() == IDOK){
		C2DOutline *pOutline = dlgParametricPreView.GetOutline();
		if (pOutline == nullptr || pOutline->GetCurveCount() == 0) {
			CDlgShowError cError(ID_ERROR_CHILD2DVIEW_NO_CURVES); //_T("Чертёж пуст или не имеет кривых"));
			return;
		}
		Set2DOutline(pOutline);
	}
}
///////////////////////////////////////////////////////////////////////////////




void CChild2DView::OnTimer(size_t nIDEvent){

	CIO2DMotion *t_motion = dynamic_cast<CIO2DMotion*>(m_pI2DSketch->GetOutline());
	if (t_motion) {
		if (t_motion->m_params().size() > 0) {
			if (t_motion->GetNodeCount() > 0) {
				size_t nNodes = t_motion->GetNodeCount() - 1;
				DBL dExitTime = t_motion->m_params()[nNodes].GetExitTime();
				if (m_pI2DSketch->m_dCurTime > dExitTime) {
					StopTimerPlay(nIDEvent);
				}
				else {
					m_pI2DSketch->ChangeTime(*m_Parameter);
					Invalidate(0);
				}
			}
		}else {
			//если пуст, то параметров нет и нечего запускать
			StopTimerPlay(nIDEvent);
			CDlgShowError cError(ID_ERROR_CHILD2DVIEW_NO_PARAMS);
		}
	} else {
		//если пуст, то траектории нет и нечего запускать
		StopTimerPlay(nIDEvent);
		CDlgShowError cError(ID_ERROR_CHILD2DVIEW_NO_MOTION);
	}
	
}

//останавливаем таймер и анимацию.
void CChild2DView::StopTimerPlay(size_t nIDEvent) {
	KillTimer(nIDEvent);
	m_pI2DSketch->m_fPlay = false;
	m_pI2DSketch->m_dCurTime = 0;
}

////////////////////////////////////////////////////////////////////////////////
void CChild2DView::OnUpdateImportPrm(CCmdUI *pCmdUI){
	// Запрещаем КЭ сетке загружать чертёж
	if (m_pDoc) {
		if (m_pDoc->GetType() == C2DMESHER) {
			pCmdUI->Enable(0);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void CChild2DView::OnImportDxf(){
	
	// Запрещаем КЭ сетке загружать чертёж
	if (m_pDoc) {
		if (m_pDoc->GetType() == C2DMESHER) {
			return;
		}
	}
	
	CFileDialog fileDlg( 1 , _T("dxf"),0, 4|2, _T("dxf files (*.dxf)|*.dxf|"));
	if(fileDlg.DoModal() == IDOK){
		CString strFile = fileDlg.GetPathName();
		CT2CA pszConvertedAnsiString(strFile);
		DXF::CGraph graph;
		if(!graph.LoadFromDXF(pszConvertedAnsiString)){
			CDlgShowError cError(ID_ERROR_CHILD2DVIEW_FAIL_LOADDXF); //_T("Не загружен LoadFromDXF"));
			return;
		}else{
			C2DOutline *pOutline = graph.MakeOutline();
			Set2DOutline(pOutline);
			delete pOutline;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void CChild2DView::OnUpdateImportDxf(CCmdUI *pCmdUI){
	// Запрещаем КЭ сетке загружать чертёж
	if (m_pDoc) {
		if (m_pDoc->GetType() == C2DMESHER) {
			pCmdUI->Enable(0);
		}
	}
}

// Действия при перемещении мыши
void CChild2DView::OnMouseMove(UNINT nFlags, CPoint point)
{
	Math::C2DPoint cursor;
	cursor = FizToLogical(point);
	
	// получаем координаты мыши и переводим их в координаты чертежа
	if (m_nCurType == 0){
		m_pI2DSketch->SetActiveObject(&cursor);
	}
	
	//Если курсор - добавление отрезка, то рисуем вспомогательные кривые
	if (m_nCurType == 101) {
		m_pI2DSketch->SetTempLine(&cursor);
	}

	//Если курсор - добавление дуги, то рисуем вспомогательную дугу
	if (m_nCurType == 102) {
		m_pI2DSketch->SetTempArc(&cursor);
	}

	//Если курсор - фаска, то строим временное скругление
	if (m_nCurType == 103){
		m_pI2DSketch->SetTempFacet(&cursor,false);
	}
	
	Invalidate(0);

	CGL2DView::OnMouseMove(nFlags, point);
}

//! По типу курсора определяем действие 
void CChild2DView::OnLButtonUp(UNINT nFlags, CPoint point)
{
	// TODO: добавьте свой код обработчика сообщений или вызов стандартного
	// Будем делать выбор выделенной кривой
	
	if (nFlags != 8) {
		m_pI2DSketch->m_fCtrl = false;
	}
	else {
		m_pI2DSketch->m_fCtrl = true;
	}

	//если курсор-стрелка
	if (m_nCurType == 0){
		
		m_pI2DSketch->OnClick(&m_nCurType);
		CMDIChildWndEx *pChildFrame = dynamic_cast<CMDIChildWndEx *>(GetParent());
		CMDIFrameWnd* pFrame = pChildFrame->GetMDIFrame();
		pFrame->SendMessage(WMR_UPDATEPANE,(WPARAM) nFlags);

	}

	//Если курсор - рука
	if (m_nCurType == 2) {
		Math::C2DPoint cursor;
		cursor = FizToLogical(point);
	}

	//Если курсор - фаска, то строим скругление угла
	if (m_nCurType == 103) {
		Math::C2DPoint cursor;
		cursor = FizToLogical(point);
		m_pI2DSketch->SetTempFacet(&cursor,true);
	}

	//если добавляем в чертёж (точка, отрезок, дуга)
	if((m_nCurType >= 100) && (m_nCurType <= 102)){
		
		m_pI2DSketch->OnClick(&m_nCurType);
		CMDIChildWndEx *pChildFrame = dynamic_cast<CMDIChildWndEx *>(GetParent());
		CMDIFrameWnd* pFrame = pChildFrame->GetMDIFrame();
		m_pI2DSketch->AddItem(&FizToLogical(point), m_nCurType);
		pFrame->SendMessage(WMR_UPDATEPANE,(WPARAM) nFlags);
	}

	//если привязка траектории
	if (m_nCurType == 200){
		CMDIChildWndEx *pChildFrame = dynamic_cast<CMDIChildWndEx *>(GetParent());
		CMDIFrameWnd* pFrame = pChildFrame->GetMDIFrame();
		
		C2DTool* pTool = dynamic_cast<C2DTool* >(m_pDoc);
		if (pTool){
			Math::C2DPoint cursor;
			cursor = FizToLogical(point);
			pTool->AddStickPoint(cursor, *m_Parameter);
			pFrame->SendMessage(WMR_UPDATEPANE,(WPARAM) nFlags);
		}
	}

	//если точка сгущения
	if (m_nCurType == 300){
		CMDIChildWndEx *pChildFrame = dynamic_cast<CMDIChildWndEx *>(GetParent());
		CMDIFrameWnd* pFrame = pChildFrame->GetMDIFrame();

		C2DMesher* pMesher = dynamic_cast<C2DMesher* >(m_pDoc);
		if (pMesher){
			Math::C2DPoint cursor;
			cursor = FizToLogical(point);
			pMesher->AddMeshPoint(cursor);
			pFrame->SendMessage(WMR_UPDATEPANE,(WPARAM) nFlags);
		}
	}
	Invalidate(0);

	CGL2DView::OnLButtonUp(nFlags, point);
}


LRESULT CChild2DView::OnSetCursor(WPARAM wParam, LPARAM lParam){
	

	ChangeCursor(wParam,false);

	return 0;
}
