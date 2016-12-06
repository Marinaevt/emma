#include "StdAfx.h"
#include "GLModelPreview.h"

#include "../RegisteredMessages.h"

//////////////////////////////////////////////////////
CGLModelPreview::CGLModelPreview()
{
	m_pModelPreviewOutline = nullptr;
}

CGLModelPreview::~CGLModelPreview(void)
{

}

//////////////////////////////////////////////////////
IMPLEMENT_DYNCREATE(CGLModelPreview,CGL2DView);

//////////////////////////////////////////////////////

BEGIN_MESSAGE_MAP(CGLModelPreview, CGL2DView)
		ON_WM_MOUSEWHEEL()
		ON_REGISTERED_MESSAGE(WMR_DRAFTRECT, OnDraftRect)
END_MESSAGE_MAP()

void CGLModelPreview::ClearModel() {
	delete m_pModelPreviewOutline;
	m_pModelPreviewOutline = nullptr;
	Invalidate();
	UpdateWindow();
}

//////////////////////////////////////////////////////
bool CGLModelPreview::SetModel(C2DOutline *pOutline)
{
	if (!pOutline || pOutline->GetNodeCount() <= 0){
		Invalidate();
		UpdateWindow();
		return false;
	}
	
	m_pModelPreviewOutline = pOutline;
	SetZoom(1.0/1.1);	//Чтобы заготовка краями не упиралась в окно
	SetDraftRect(m_pModelPreviewOutline->GetBoundingBox(), true);
	SetProjection();
	Invalidate();
	UpdateWindow();
	pOutline->DelRepNodes(EPS);

	return true;
}
//////////////////////////////////////////////////////
void CGLModelPreview::DrawScene()
{
	glEnable(GL_ALPHA_TEST);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	if (m_bShowGrid) {
		DrawGrid();
	}

	glColor3f(1.0f,0.5f,0.0f);	//orange color
	//glLineWidth(2);	//закомментировал, т.к. иначе не видно узлов и возникают недопонимания...
	
	if(m_pModelPreviewOutline != nullptr && m_pModelPreviewOutline->GetNodeCount() > 0){
		m_pModelPreviewOutline->DrawGL(*m_Parameter);
	}
}
//////////////////////////////////////////////////////
LRESULT CGLModelPreview::OnDraftRect(WPARAM wParam, LPARAM lParam){
	CRect2D *rect = (CRect2D *)wParam;
	SetDraftRect(*rect, lParam != 0);
	return 0;
}

//////////////////////////////////////////////////////
int CGLModelPreview::PreCreateWindow(CREATESTRUCT& cs){
	if (!CGL2DView::PreCreateWindow(cs))
		return 0;

	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(nullptr, IDC_ARROW), reinterpret_cast<HBRUSH>(COLOR_WINDOW+1), nullptr);

	return 1;
}