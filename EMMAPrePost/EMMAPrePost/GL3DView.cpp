#include "StdAfx.h"
#include "GL3DView.h"


#define rad 57.2957795

#define sq3 0.5773502691

BEGIN_MESSAGE_MAP(CGL3DView, CGLBasicView)
	ON_WM_MOUSEWHEEL()
	ON_WM_SETCURSOR()
	ON_WM_LBUTTONDOWN()
	ON_WM_RBUTTONDOWN()
	ON_WM_MOUSEMOVE()

	ON_COMMAND(ID_ARROW, &CGL3DView::OnArrow)
	ON_UPDATE_COMMAND_UI(ID_ARROW, &CGL3DView::OnUpdateArrow)
	ON_COMMAND(ID_BUTTON_ZOOM, &CGL3DView::OnZoom)
	ON_UPDATE_COMMAND_UI(ID_BUTTON_ZOOM, &CGL3DView::OnUpdateZoom)
	ON_COMMAND(ID_HAND, &CGL3DView::OnHand)
	ON_UPDATE_COMMAND_UI(ID_HAND, &CGL3DView::OnUpdateHand)
	ON_COMMAND(ID_AXIS, &CGL3DView::OnAxis)
	ON_UPDATE_COMMAND_UI(ID_AXIS, &CGL3DView::OnUpdateAxis)
	ON_COMMAND(ID_GRID, &CGL3DView::OnGrid)
	ON_UPDATE_COMMAND_UI(ID_GRID, &CGL3DView::OnUpdateGrid)
	ON_COMMAND(ID_COLOR, &CGL3DView::OnColor)
	ON_UPDATE_COMMAND_UI(ID_COLOR, &CGL3DView::OnUpdateColor)

	ON_COMMAND(ID_LIGHT, &CGL3DView::OnLight)
	ON_UPDATE_COMMAND_UI(ID_LIGHT, &CGL3DView::OnUpdateLight)

	ON_COMMAND(ID_PERSPECTIVE, &CGL3DView::OnPerspective)
	ON_UPDATE_COMMAND_UI(ID_PERSPECTIVE, &CGL3DView::OnUpdatePerspective)

	ON_COMMAND(ID_FLY, &CGL3DView::OnFly)
	ON_UPDATE_COMMAND_UI(ID_FLY, &CGL3DView::OnUpdateFly)

	ON_COMMAND(ID_ROTATE_FREE, &CGL3DView::OnRotateFree)
	ON_UPDATE_COMMAND_UI(ID_ROTATE_FREE, &CGL3DView::OnUpdateRotateFree)

	ON_COMMAND(ID_SHOWALL, &CGL3DView::OnShowAll)
	ON_UPDATE_COMMAND_UI(ID_SHOWALL, &CGL3DView::OnUpdateShowAll)

	ON_WM_KEYDOWN()
	ON_WM_KEYUP()
	ON_WM_TIMER()

END_MESSAGE_MAP()

CGL3DView::CGL3DView(void)
//Борхес ковырялся здесь  (-1, 2, -1, 1, -1, 1)
//плохо - надо забирать нормальные данные
	:m_DraftRect(0, 20, 0, 50, 0, 8), m_dZoom(1), VectorCam(sq3, sq3, sq3), alpha(45), DirectionUp(0, 1, 0),
	m_bShowAxis(1), m_bShowGrid(1), m_bShowPerspective(0), Is2DViewShow(false), 
	IsWPressed(false),	IsAPressed(false), IsSPressed(false), IsDPressed(false), 
	IsShiftPressed(false), IsLeftPressed(false), IsRightPressed(false)
{
	m_ShowLight3Dobj = true;
	IsCentr3DShow = true;
	IsLightOn = true;
	
	CenterScene.x = (m_DraftRect.m_minX + m_DraftRect.m_maxX) / 2;
	CenterScene.y = (m_DraftRect.m_minY + m_DraftRect.m_maxY) / 2;
	CenterScene.z = (m_DraftRect.m_minZ + m_DraftRect.m_maxZ) / 2;

	LookCenter = CenterScene; //помещаем фокус камеры в центр

	RScene = sqrt((m_DraftRect.m_maxX - m_DraftRect.m_minX)*(m_DraftRect.m_maxX - m_DraftRect.m_minX)
		+ (m_DraftRect.m_maxY - m_DraftRect.m_minY)*(m_DraftRect.m_maxY - m_DraftRect.m_minY)
		+ (m_DraftRect.m_maxZ - m_DraftRect.m_minZ)*(m_DraftRect.m_maxZ - m_DraftRect.m_minZ)) / 2;

	DistCam_default = RScene / sin(alpha / (2 * rad));

	DistCam = DistCam_default;

	LookPos = LookCenter + DistCam*VectorCam;

	double x0, y0, z0;
	double x1, y1, z1;
	double u1, u2, u3;
	double v1, v2, v3;
	double A, B, C;
	double NormU, NormV;
	CPoint3D Log_0, Log0, LogU;

	Log_0.x = 0;
	Log_0.y = 0;
	Log_0.z = 0;

	x1 = LookPos.x;
	y1 = LookPos.y;
	z1 = LookPos.z;

	x0 = LookCenter.x;
	y0 = LookCenter.y;
	z0 = LookCenter.z;

	A = x0 - x1;
	B = y0 - y1;
	C = z0 - z1;

	Log0 = ProjToSq(Log_0, A, B, C);
	LogU = ProjToSq(DirectionUp, A, B, C);

	u1 = LogU.x - Log0.x;
	u2 = LogU.y - Log0.y;
	u3 = LogU.z - Log0.z;

	v1 = u2*C - u3*B;
	v2 = u3*A - u1*C;
	v3 = u1*B - u2*A;

	NormU = sqrt(u1*u1 + u2*u2 + u3*u3);
	NormV = sqrt(v1*v1 + v2*v2 + v3*v3);

	u1 /= NormU;
	u2 /= NormU;
	u3 /= NormU;

	v1 /= -NormV;
	v2 /= -NormV;
	v3 /= -NormV;

	DirectionUp.x = u1;
	DirectionUp.y = u2;
	DirectionUp.z = u3;

	DirectionRight.x = v1;
	DirectionRight.y = v2;
	DirectionRight.z = v3;

}

CGL3DView::~CGL3DView(void) {
}

///////////////////////////////////////////////////////////////////

void CGL3DView::SetProjection(void) {

	CRect clientRect;
	GetClientRect(&clientRect);

	LookPos = LookCenter + DistCam*VectorCam;

	if ((LookCenter - CenterScene).Norm2() > RScene*RScene) LimitArea();

	//m_dZoom=DistCam/DistCam_default;


	if (SetCurrentContext()) {

		SetViewport();

		glEnable(GL_ALPHA_TEST);
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(LookPos.x, LookPos.y, LookPos.z, LookCenter.x, LookCenter.y, LookCenter.z, DirectionUp.x, DirectionUp.y, DirectionUp.z);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		if (m_bShowPerspective)
		{
			double a, b;

			if (clientRect.Width() < clientRect.Height())
			{
				a = 2 * rad*atan(static_cast<double>(clientRect.Height())*tan(alpha / (2 * rad)) / clientRect.Width());
				b = alpha / a;
			}

			else
			{
				a = alpha;
				b = 2 * rad*atan(static_cast<double>(clientRect.Width())*tan(a / (2 * rad)) / clientRect.Height()) / a;
			}

			double persp_param = 50, dd = DistCam + 10 * RScene;

			gluPerspective(a, b, dd / persp_param, dd);
		}

		else
		{
			CRect2D rectScreen;
			GetScreenRect(rectScreen);
			glOrtho(rectScreen.left, rectScreen.right, rectScreen.bottom, rectScreen.top, 0, DistCam + 10 * RScene);
		}

		::wglMakeCurrent(nullptr, nullptr);

	}

}

///////////////////////////////////////////////////////////////////
void CGL3DView::DrawScene(void) {

	if (m_bShowGrid)DrawGrid();
	if (m_bShowAxis)DrawAxis();
	//SetDraftRect(0,30,0,30,0,4);
	//SetDraftRect(0,20,0,50,0,4);
	SetDraftRect(0, 20, 0, 30, 0, 8);
	//m_pDoc->DrawGL(*m_Parameter);
	m_pDoc->DrawGL3D(*m_Parameter);        //***///
	
}

///////////////////////////////////////////////////////////////////
bool CGL3DView::GetScreenRect(CRect2D &rect) {
	CRect clientRect;
	GetClientRect(&clientRect);

	double dZoom = DistCam / DistCam_default;


	if (clientRect.Width() == 0 || clientRect.Height() == 0)return false;
	double w = clientRect.Width(), h = clientRect.Height();
	double dWidth, dHeigth, asd;

	asd = max(dZoom, RScene / (DistCam_default + RScene));

	if (w >= h)
	{
		dHeigth = 2 * RScene*asd;
		dWidth = dHeigth*w / h;
	}
	else
	{
		dWidth = 2 * RScene*asd;
		dHeigth = dWidth*h / w;
	}

	rect.left = -dWidth*0.5;
	rect.right = dWidth*0.5;
	rect.bottom = -dHeigth*0.5;
	rect.top = dHeigth*0.5;

	return true;
}

///////////////////////////////////////////////////////////////////
void CGL3DView::DrawGrid() {

	glLineWidth(0.5);
	glBegin(GL_LINES);
	glColor3d(0.7, 0.7, 0.7);

	int i, 
		n_x = 0, n_y = 0, 
		n_z = 0, n = 25;
	double asd;

	if (m_DraftRect.dX() <= m_DraftRect.dY() && m_DraftRect.dX() <= m_DraftRect.dZ())
	{
		n_x = n;
		n_y = static_cast<int>(n*m_DraftRect.dY() / m_DraftRect.dX());
		n_z = static_cast<int>(n*m_DraftRect.dZ() / m_DraftRect.dX());
		asd = m_DraftRect.dX() / n;
	}

	if (m_DraftRect.dY() <= m_DraftRect.dX() && m_DraftRect.dY() <= m_DraftRect.dZ())
	{
		n_y = n;
		n_x = static_cast<int>(n*m_DraftRect.dX() / m_DraftRect.dY());
		n_z = static_cast<int>(n*m_DraftRect.dZ() / m_DraftRect.dY());
		asd = m_DraftRect.dY() / n;
	}

	if (m_DraftRect.dZ() <= m_DraftRect.dX() && m_DraftRect.dZ() <= m_DraftRect.dY())
	{
		n_z = n;
		n_x = static_cast<int>(n*m_DraftRect.dX() / m_DraftRect.dZ());
		n_y = static_cast<int>(n*m_DraftRect.dY() / m_DraftRect.dZ());
		asd = m_DraftRect.dZ() / n;
	}

	//===============================================

	if ((LookPos.x - LookCenter.x) >= 0)
	{
		for (i = 1;i <= n_y;i++)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY + i*asd, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY + i*asd, m_DraftRect.m_maxZ);
		}

		for (i = 1;i <= n_z;i++)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_minZ + i*asd);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_minZ + i*asd);
		}
	}

	else
	{
		for (i = 1;i <= n_y;i++)
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY + i*asd, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY + i*asd, m_DraftRect.m_maxZ);
		}

		for (i = 1;i <= n_z;i++)
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_minZ + i*asd);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_minZ + i*asd);
		}

	}

	//==============================================
	if ((LookPos.y - LookCenter.y) >= 0)
	{
		for (i = 1;i <= n_x;i++)
		{
			glVertex3d(m_DraftRect.m_minX + i*asd, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX + i*asd, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
		}

		for (i = 1;i <= n_z;i++)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_minZ + i*asd);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_minZ + i*asd);
		}
	}

	else
	{
		for (i = 1;i <= n_x;i++)
		{
			glVertex3d(m_DraftRect.m_minX + i*asd, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX + i*asd, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}

		for (i = 1;i <= n_z;i++)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_minZ + i*asd);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_minZ + i*asd);
		}

	}

	//===============================================

	if ((LookPos.z - LookCenter.z) >= 0)
	{
		for (i = 1;i <= n_x;i++)
		{
			glVertex3d(m_DraftRect.m_minX + i*asd, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX + i*asd, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
		}

		for (i = 1;i <= n_y;i++)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY + i*asd, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY + i*asd, m_DraftRect.m_minZ);
		}
	}

	else
	{
		for (i = 1;i <= n_x;i++)
		{
			glVertex3d(m_DraftRect.m_minX + i*asd, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_minX + i*asd, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}

		for (i = 1;i <= n_y;i++)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY + i*asd, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY + i*asd, m_DraftRect.m_maxZ);
		}

	}
	//===============================================

	glEnd();
	glFlush();
}

///////////////////////////////////////////////////////////////////
void CGL3DView::DrawAxis() {

	//X - красная ось, Y - зелёная ось, Z - синяя ось
	glLineWidth(2);

	glBegin(GL_LINES);

	if ((LookPos.z - LookCenter.z) >= 0)
	{
		glColor3f(0, 1, 0);
		if ((LookPos.x - LookCenter.x) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
		}
		else
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
		}

		glColor3f(1, 0, 0);
		if ((LookPos.y - LookCenter.y) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
		}
		else
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
		}

	}

	else
	{
		glColor3f(0, 1, 0);
		if ((LookPos.x - LookCenter.x) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}
		else
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}

		glColor3f(1, 0, 0);
		if ((LookPos.y - LookCenter.y) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
		}
		else
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}

	}

	if ((LookPos.x - LookCenter.x) >= 0)
	{
		glColor3f(0, 1, 0);

		if ((LookPos.z - LookCenter.z) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
		}

		else
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}


		glColor3f(0, 0, 1);
		if ((LookPos.y - LookCenter.y) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
		}
		else
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}
	}

	else
	{
		glColor3f(0, 1, 0);

		if ((LookPos.z - LookCenter.z) >= 0)
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
		}

		else
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}

		glColor3f(0, 0, 1);
		if ((LookPos.y - LookCenter.y) >= 0)
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
		}
		else
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}
	}

	if ((LookPos.y - LookCenter.y) >= 0)
	{
		glColor3f(0, 0, 1);
		if ((LookPos.x - LookCenter.x) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
		}

		else
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
		}

		glColor3f(1, 0, 0);
		if ((LookPos.z - LookCenter.z) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_minZ);
		}

		else
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_minY, m_DraftRect.m_maxZ);
		}

	}

	else
	{
		glColor3f(0, 0, 1);
		if ((LookPos.x - LookCenter.x) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}

		else
		{
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}

		glColor3f(1, 0, 0);
		if ((LookPos.z - LookCenter.z) >= 0)
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_minZ);
		}

		else
		{
			glVertex3d(m_DraftRect.m_minX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
			glVertex3d(m_DraftRect.m_maxX, m_DraftRect.m_maxY, m_DraftRect.m_maxZ);
		}
	}

	glEnd();
	glFlush();
}

///////////////////////////////////////////////////////////////////
void CGL3DView::ShowAxis(bool bShow, bool bUpdate) {
	m_bShowAxis = bShow;
	if (bUpdate) {
		Invalidate(0);
		SetProjection();
		UpdateWindow();
	}
}

///////////////////////////////////////////////////////////////////
void CGL3DView::ShowGrid(bool bShow, bool bUpdate) {
	m_bShowGrid = bShow;
	if (bUpdate) {
		Invalidate(0);
		SetProjection();
		UpdateWindow();
	}
}

void CGL3DView::ShowPerspective(bool bShow, bool bUpdate) {
	m_bShowPerspective = bShow;
	if (bUpdate) {
		Invalidate(0);
		SetProjection();
		UpdateWindow();
	}
}

void CGL3DView::ChangeCursor(int curType, bool bUpdate) {
	m_nCurType = curType;
	if (bUpdate) {
		Invalidate(0);
		UpdateWindow();
	}
}

///////////////////////////////////////////////////////////////////
void CGL3DView::SetZoom(double dZoom, bool bUpdate) {
	if (dZoom <= 0) { 
		CDlgShowError cError(ID_ERROR_GL3DVIEW_ZOOM_WRONG); //_T("Масштаб меньше нуля"));
		return; 
	}

	m_dZoom = dZoom;

	DistCam = m_dZoom*DistCam_default;

	if (bUpdate) {
		SetProjection();
		Invalidate(0);
		UpdateWindow();
	}
}


int CGL3DView::OnMouseWheel(UNINT nFlags, short zDelta, CPoint pt)
{
	//CPoint3D asd;	//не используется
	if (zDelta < 0)
	{
		SetZoom(GetZoom()*1.1, false);
		SetProjection();
		Invalidate(0);
	}
	else
	{
		SetZoom(GetZoom() / 1.1, false);
		SetProjection();
		Invalidate(0);
	}

	return CGLBasicView::OnMouseWheel(nFlags, zDelta, pt);
}

int CGL3DView::OnSetCursor(CWnd* pWnd, UNINT nHitTest, UNINT message)
{
	switch (m_nCurType)
	{
	case 1:
	{
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_ZOOM));
		return 1;
	}
	break;

	case 2:
	{
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_HAND_O));
		return 1;
	}
	break;

	case 3:
	{
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_FLY));
		return 1;
	}
	break;

	case 4:
	{
		::SetCursor(AfxGetApp()->LoadCursorW(IDC_ROTATE_FREE));
		return 1;
	}
	break;

	}

	return CGLBasicView::OnSetCursor(pWnd, nHitTest, message);
}

void CGL3DView::OnLButtonDown(UNINT nFlags, CPoint point)
{
	switch (m_nCurType)
	{
	case 1:
	{
		//C3DPoint asd;
		SetZoom(GetZoom()*0.5, false);
		SetProjection();
		Invalidate(0);
	}
	break;

	case 2:
	{

	}
	break;

	case 3:
	{
	}
	break;

	case 4:
	{
		m_pPrevPos = point;
		m_pPrevPosSp = FizToLogical3D(point);

	}
	break;

	}

	CGLBasicView::OnLButtonDown(nFlags, point);
}

void CGL3DView::OnRButtonDown(UNINT nFlags, CPoint point)
{
	if (m_nCurType == 1)
	{
		//C3DPoint asd;
		SetZoom(GetZoom() * 2, false);

		SetProjection();
		Invalidate(0);
	}

	CPoint3D Log;
	Log = FizToLogical3D(point);

	CGLBasicView::OnRButtonDown(nFlags, point);
}



CPoint3D CGL3DView::ProjToSq(CPoint3D point, double A, double B, double C)
{
	CPoint3D proj;

	double x0 = point.x, y0 = point.y, z0 = point.z;
	double x1 = LookPos.x, y1 = LookPos.y, z1 = LookPos.z;
	double p;

	p = (A*(x1 - x0) + B*(y1 - y0) + C*(z1 - z0)) / (A*A + B*B + C*C);

	proj.x = x0 + p*A;
	proj.y = y0 + p*B;
	proj.z = z0 + p*C;

	return proj;
}

CPoint3D CGL3DView::FizToLogical3D(CPoint FizPoint)
{
	CPoint3D Log;
	CRect clientRect;
	GetClientRect(&clientRect);

	CRect2D rectScreen;
	GetScreenRect(rectScreen);

	double x0, y0, z0;
	double x1, y1, z1;
	double u1, u2, u3;
	double v1, v2, v3;
	double A, B, C;
	double NormU, NormV;
	double Scr_Cli;
	CPoint3D LogC, Log_0, Log0, Log0_, LogU, LogV;

	Log_0.x = 0;
	Log_0.y = 0;
	Log_0.z = 0;

	x1 = LookPos.x;
	y1 = LookPos.y;
	z1 = LookPos.z;

	x0 = LookCenter.x;
	y0 = LookCenter.y;
	z0 = LookCenter.z;

	A = x0 - x1;
	B = y0 - y1;
	C = z0 - z1;

	LogC = ProjToSq(LookCenter, A, B, C);
	Log0 = ProjToSq(Log_0, A, B, C);
	LogU = ProjToSq(DirectionUp, A, B, C);

	u1 = LogU.x - Log0.x;
	u2 = LogU.y - Log0.y;
	u3 = LogU.z - Log0.z;

	v1 = u2*C - u3*B;
	v2 = u3*A - u1*C;
	v3 = u1*B - u2*A;

	NormU = sqrt(u1*u1 + u2*u2 + u3*u3);
	NormV = sqrt(v1*v1 + v2*v2 + v3*v3);

	u1 /= NormU;
	u2 /= NormU;
	u3 /= NormU;

	v1 /= -NormV;
	v2 /= -NormV;
	v3 /= -NormV;

	DirectionUp.x = u1;
	DirectionUp.y = u2;
	DirectionUp.z = u3;

	DirectionRight.x = v1;
	DirectionRight.y = v2;
	DirectionRight.z = v3;

	Scr_Cli = rectScreen.Width() / static_cast<double>(clientRect.Width());

	Log0_.x = LogC.x + u1*rectScreen.Heigth() / 2 - v1*rectScreen.Width() / 2;
	Log0_.y = LogC.y + u2*rectScreen.Heigth() / 2 - v2*rectScreen.Width() / 2;
	Log0_.z = LogC.z + u3*rectScreen.Heigth() / 2 - v3*rectScreen.Width() / 2;

	Log.x = v1*FizPoint.x*Scr_Cli - u1*FizPoint.y*Scr_Cli + Log0_.x;
	Log.y = v2*FizPoint.x*Scr_Cli - u2*FizPoint.y*Scr_Cli + Log0_.y;
	Log.z = v3*FizPoint.x*Scr_Cli - u3*FizPoint.y*Scr_Cli + Log0_.z;

	return Log;
}

double CGL3DView::ToUgol(CPoint3D p1, CPoint3D p2)
{
	return acos((p1.x*p2.x + p1.y*p2.y + p1.z*p2.z)
		/ (sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z)
			*sqrt(p2.x*p2.x + p2.y*p2.y + p2.z*p2.z)));

}

void CGL3DView::OnMouseMove(UNINT nFlags, CPoint point)
{
	if (m_nCurType == 4 && m_bLButtonPressed)
	{
		double Ug;
		double param = 0.05;

		CPoint3D NV2 = FizToLogical3D(point);
		CPoint3D NV3 = NV2;
		CPoint3D a, b;
		double mnoj = 100;
		CPoint3D Ve;

		a = mnoj*(NV3 - m_pPrevPosSp);

		b = LookPos - LookCenter;

		Ve.x = a.y*b.z - a.z*b.y;
		Ve.y = a.z*b.x - a.x*b.z;
		Ve.z = a.x*b.y - a.y*b.x;

		if (Ve.Trans(a, Ve, M_PI, LookCenter).x / b.x <= 0)
		{
			Ve *= -1;
		}

		Ug = sqrt((NV3.x - m_pPrevPosSp.x)*(NV3.x - m_pPrevPosSp.x)
			+ (NV3.y - m_pPrevPosSp.y)*(NV3.y - m_pPrevPosSp.y)
			+ (NV3.z - m_pPrevPosSp.z)*(NV3.z - m_pPrevPosSp.z))*param*RScene / (DistCam + RScene);


		if ((ScalPr(Ve, DirectionUp) > 0 && ScalPr(a, DirectionRight) < 0) ||
			(ScalPr(Ve, DirectionUp) < 0 && ScalPr(a, DirectionRight) > 0)) Ug *= -1;



		LookPos = LookPos.Trans(LookPos, Ve, Ug, CenterScene);
		LookCenter = LookCenter.Trans(LookCenter, Ve, Ug, CenterScene);
		VectorCam = (LookPos - LookCenter) / DistCam;

		//Ug*=-1;
		//DirectionUp=DirectionUp.Trans(DirectionUp, Ve, Ug, CenterScene);

		m_pPrevPos = point;
		m_pPrevPosSp = NV3;

		SetProjection();
		Invalidate(0);
	}

	//TODO: Борхес ковырял это
	if (m_nCurType == 2 && m_bLButtonPressed)
	{
		POINT cp;

		GetCursorPos(&cp);

		//сомнительно: разобраться
		double param = 0.00005, a;
		int dx = cp.x - 600, dy = cp.y - 400;
		CPoint3D LookCenter2;

		LookCenter2 = LookCenter + param*DistCam*dx*DirectionRight - param*DistCam*dy*DirectionUp;

		a = DistCam / (LookCenter2 - LookPos).Norm();

		LookCenter = LookPos + a*(LookCenter2 - LookPos);

		VectorCam = (LookPos - LookCenter) / DistCam;
		//////////////////
		SetProjection();
		Invalidate(0);
	}

	if (m_nCurType == 3)
	{
		POINT cp;

		GetCursorPos(&cp);

		double param = 0.00005, a;
		int dx = cp.x - 600, dy = cp.y - 400;
		CPoint3D LookCenter2;

		LookCenter2 = LookCenter + param*DistCam*dx*DirectionRight - param*DistCam*dy*DirectionUp;

		a = DistCam / (LookCenter2 - LookPos).Norm();

		LookCenter = LookPos + a*(LookCenter2 - LookPos);

		VectorCam = (LookPos - LookCenter) / DistCam;

		SetProjection();
		Invalidate(0);

	}

	CGLBasicView::OnMouseMove(nFlags, point);
}


void CGL3DView::OnArrow() {
	ChangeCursor(0, true);

	Invalidate(0);
	SetProjection();
	UpdateWindow();
}

void CGL3DView::OnUpdateArrow(CCmdUI *pCmdUI) {

	Invalidate(0);
	SetProjection();
	UpdateWindow();
}
void CGL3DView::OnZoom() {
	if (GetCursorType() == 1) ChangeCursor(0, true);
	else ChangeCursor(1, true);

	Invalidate(0);
	SetProjection();
	UpdateWindow();
}

void CGL3DView::OnUpdateZoom(CCmdUI *pCmdUI) {
	if (GetCursorType() == 1) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);

	Invalidate(0);
	SetProjection();
	UpdateWindow();
}

void CGL3DView::OnHand() {
	if (GetCursorType() == 2) ChangeCursor(0, true);
	else ChangeCursor(2, true);

	Invalidate(0);
	SetProjection();
	UpdateWindow();
}
void CGL3DView::OnUpdateHand(CCmdUI *pCmdUI) {
	if (GetCursorType() == 2) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);

	Invalidate(0);
	SetProjection();
	UpdateWindow();
}

void CGL3DView::OnAxis() {
	ShowAxis(!IsAxisVisible());
}
void CGL3DView::OnUpdateAxis(CCmdUI *pCmdUI) {
	pCmdUI->SetCheck(IsAxisVisible());
}

void CGL3DView::OnGrid() {
	ShowGrid(!IsGridVisible());
}
void CGL3DView::OnUpdateGrid(CCmdUI *pCmdUI) {
	pCmdUI->SetCheck(IsGridVisible());
}

void CGL3DView::OnColor() {

}
void CGL3DView::OnUpdateColor(CCmdUI *pCmdUI) {

}

void CGL3DView::OnLight() {

}
void CGL3DView::OnUpdateLight(CCmdUI *pCmdUI) {

}

void CGL3DView::OnPerspective() {
	ShowPerspective(!IsPerspective());
}
void CGL3DView::OnUpdatePerspective(CCmdUI *pCmdUI) {
	pCmdUI->SetCheck(IsPerspective());
}


void CGL3DView::OnFly() {
	SetTimer(50, 10, nullptr);

	ChangeCursor(3, true);
	::SetCursorPos(600, 400);

	Invalidate(0);
	SetProjection();
	UpdateWindow();
}
void CGL3DView::OnUpdateFly(CCmdUI *pCmdUI) {
	if (GetCursorType() == 3) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);

	Invalidate(0);
	SetProjection();
	UpdateWindow();
}

void CGL3DView::OnRotateFree() {
	if (GetCursorType() == 4) ChangeCursor(0, true);
	else ChangeCursor(4, true);

	Invalidate(0);
	SetProjection();
	UpdateWindow();

}
void CGL3DView::OnUpdateRotateFree(CCmdUI *pCmdUI) {
	if (GetCursorType() == 4) pCmdUI->SetCheck(1);
	else pCmdUI->SetCheck(0);

	Invalidate(0);
	SetProjection();
	UpdateWindow();

}


void CGL3DView::OnKeyDown(UNINT nChar, UNINT nRepCnt, UNINT nFlags)
{
	if (m_nCurType == 3)
	{
		switch (nChar)
		{
		case 87://W
			IsWPressed = true;
			break;

		case 65://A
			IsAPressed = true;
			break;

		case 83://S
			IsSPressed = true;
			break;

		case 68://D
			IsDPressed = true;
			break;

		case VK_SHIFT://Shift
			IsShiftPressed = true;
			break;

		case VK_LEFT://<-
			IsLeftPressed = true;
			break;

		case VK_RIGHT://->
			IsRightPressed = true;
			break;

		case VK_ESCAPE://Esc
			ChangeCursor(0, true);
			break;

		}

	}

	CGLBasicView::OnKeyDown(nChar, nRepCnt, nFlags);
}

void CGL3DView::OnKeyUp(UNINT nChar, UNINT nRepCnt, UNINT nFlags)
{
	switch (nChar)
	{
	case 87://W
		IsWPressed = false;
		break;

	case 65://A
		IsAPressed = false;
		break;

	case 83://S
		IsSPressed = false;
		break;

	case 68://D
		IsDPressed = false;
		break;

	case VK_SHIFT://Shift
		IsShiftPressed = false;
		break;

	case VK_LEFT://<-
		IsLeftPressed = false;
		break;

	case VK_RIGHT://->
		IsRightPressed = false;
		break;
	}

	CGLBasicView::OnKeyUp(nChar, nRepCnt, nFlags);
}


void CGL3DView::OnTimer(size_t nIDEvent)
{
	double param = 0.005;

	CPoint3D Ve;
	double Ug;

	if (GetCursorType() == 3)
	{
		::SetCursorPos(600, 400);

		if (IsWPressed)
		{
			CPoint3D asd = DirectionUp*param*DistCam;

			LookCenter += asd;

			Invalidate(0);
			SetProjection();
		}

		if (IsSPressed)
		{
			CPoint3D asd = DirectionUp*param*DistCam;

			LookCenter -= asd;

			Invalidate(0);
			SetProjection();
		}

		if (IsAPressed)
		{
			CPoint3D asd = DirectionRight*param*DistCam;

			LookCenter -= asd;

			Invalidate(0);
			SetProjection();
		}

		if (IsDPressed)
		{
			CPoint3D asd = DirectionRight*param*DistCam;

			LookCenter += asd;

			Invalidate(0);
			SetProjection();
		}

		if (IsLeftPressed)
		{
			Ug = -0.05;

			DirectionUp = DirectionUp.Trans(DirectionUp, VectorCam, Ug, CenterScene);
		}

		if (IsRightPressed)
		{
			Ug = 0.05;

			DirectionUp = DirectionUp.Trans(DirectionUp, VectorCam, Ug, CenterScene);
		}

		if (m_nCurType == 3)
		{

			if (m_bLButtonPressed)
			{
				param = 0.05;

				SetZoom(GetZoom() / (1 + param), false);
				SetProjection();
				Invalidate(0);

			}

			if (m_bRButtonPressed)
			{
				param = 0.05;

				SetZoom(GetZoom()*(1 + param), false);
				SetProjection();
				Invalidate(0);

			}
		}
	}

	else KillTimer(50);

	CGLBasicView::OnTimer(nIDEvent);
}


void CGL3DView::LimitArea()
{
	double D = (LookCenter - CenterScene).Norm();
	double a = RScene / D;
	CPoint3D asd;
	POINT cp;

	asd = (1 - a)*(CenterScene - LookCenter);

	GetCursorPos(&cp);

	LookCenter += asd;
	if (cp.x == 600 && cp.y == 400) LookPos += asd;

}
//TODO: подумать над реализацией 3Д Box без ссылки на 2Д 
void CGL3DView::SetDraftRect(double m_X, double mX, double m_Y, double mY, double m_Z, double mZ)
{
	/*Math::C2DRect rect;
	if (m_pDoc->GetBoundingBox(rect)) {	//чертёж 2д

		m_DraftRect = CRect3D(rect.left, rect.right, rect.bottom, rect.top, rect.left, rect.right);//CRect3D  Z=0 20
																								   //{m_X mX} = {m_Z mZ} справедливо в случае осесимметрической задачи , в других - параметр 

	}
	else {
	//*/
		m_DraftRect = CRect3D(m_X, mX, m_Y, mY, m_Z, mZ);
	//}
}


void CGL3DView::SetDraftRect(const CRect3D &rect, bool bCenter) {
	m_DraftRect = rect;
	if (bCenter) {
		//TODO: 3D SetCenter IN SetDraftRect

		//m_dCx = (rect.right + rect.left)*0.5;
		//m_dCy = (rect.top + rect.bottom)*0.5;
		//SetCentr(m_dCx, m_dCy);
	}
}
/*
void CGL3DView::OnShowAll(){

CRect3D m_DraftRect;//

if (m_pDoc->Get3DBoundingBox(m_DraftRect)){ //Какие параметры у BoundingBox?
SetDraftRect(0,60,0,10,0,8);
}else{
SetDraftRect(0, 30, 0, 30, 0, 30);	//заглушка, если ничего не загружено
}
SetZoom(1/1.1);

}
void CGL3DView::OnUpdateShowAll(CCmdUI *pCmdUI){
}*/
//--------------------------------------------------
/*LRESULT CGL3DView::OnInitialUpdate(WPARAM wParam, LPARAM lParam){
if(!m_pDoc)return 0;
SetZoom(1);
//CRectD rect;
CRect3D Rect;
if(m_pDoc->Get3DBoundingBox(Rect)) {
SendMessage(WMR_DRAFTRECT, (WPARAM)&Rect, (LPARAM)1);
}else{
//Если документ ещё не загружен, то ставим вид по-умолчанию
Rect = CRect3D(0, 30, 0, 10, 0, 4);
SendMessage(WMR_DRAFTRECT, (WPARAM)&Rect, (LPARAM)1);
};
return 0;
}*/

//------------------------------------------------------
void CGL3DView::On2DView() {

	Is2DViewShow = !Is2DViewShow;

}

void CGL3DView::OnUpdate2DView(CCmdUI *pCmdUI) {

	pCmdUI->SetCheck(Is2DViewShow);
	if (Is2DViewShow) {


	}
	else
	{


	}

}
//--------------------------------------
/*void CGL3DView::SetCentr(double x, double y,double z){

CRect2D scrRect;
GetScreenRect(scrRect);

if(m_DraftRect.minX <= scrRect.Width()){
m_dCx = min(max(x, m_DraftRect.left - scrRect.Width()*0.5), m_DraftRect.right + scrRect.Width()*0.5);
}else{
m_dCx = max(min(x, m_DraftRect.right + scrRect.Width()*0.5), m_DraftRect.left - scrRect.Width()*0.5);
}

if(m_DraftRect.minY <= scrRect.Heigth()){
m_dCy = min(max(y, m_DraftRect.bottom - scrRect.Heigth()*0.5), m_DraftRect.top + scrRect.Heigth()*0.5);
}else{
m_dCy = max(min(y, m_DraftRect.top + scrRect.Heigth()*0.5), m_DraftRect.bottom - scrRect.Heigth()*0.5);
}
m_dCz = 1;
}*/

void CGL3DView::OnShowAll() {

	CRect3D rect;

	if (m_pDoc->GetBoundingBox(rect)) {
		SetDraftRect(rect, true);
	}
	else {
		// X Y Z
		SetDraftRect(CRect3D(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0), true);	//заглушка, если ничего не загружено
	}
	SetZoom(1 / 1.1);
}

void CGL3DView::OnUpdateShowAll(CCmdUI *pCmdUI) {
}

void CGL3DView::OnCentr3D() {

	IsCentr3DShow = !IsCentr3DShow;

}

void CGL3DView::OnUpdateCentr3D(CCmdUI *pCmdUI) {

	pCmdUI->SetCheck(IsCentr3DShow);



	if (IsCentr3DShow) {
		// SetCentr(m_dCx, m_dCy,m_dCz);
		/*CenterScene.x=(m_DraftRect.minX);
		CenterScene.y=(m_DraftRect.minY);
		CenterScene.z=(m_DraftRect.minZ);

		LookCenter3=CenterScene; //помещаем фокус камеры в центр

		RScene=sqrt((m_DraftRect.maxX-m_DraftRect.minX)*(m_DraftRect.maxX-m_DraftRect.minX)
		+(m_DraftRect.maxY-m_DraftRect.minY)*(m_DraftRect.maxY-m_DraftRect.minY)
		+(m_DraftRect.maxZ-m_DraftRect.minZ)*(m_DraftRect.maxZ-m_DraftRect.minZ))/2;

		DistCam_default=RScene / sin(alpha/(2*rad));

		DistCam=DistCam_default;

		LookPos=LookCenter+DistCam*VectorCam;

		SetProjection();
		Invalidate(0);



		GetScreenRect(rectScreen);
		glOrtho(rectScreen.left,rectScreen.right,rectScreen.bottom,rectScreen.top,0,DistCam+10*RScene);


		CenterScene.x=0; //(m_DraftRect.minX+m_DraftRect.maxX)/2;
		CenterScene.y=0;// (m_DraftRect.minY+m_DraftRect.maxY)/2;
		CenterScene.z=0;// (m_DraftRect.minZ+m_DraftRect.maxZ)/2;

		LookCenter=CenterScene; //помещаем фокус камеры в центр

		RScene=sqrt((m_DraftRect.maxX-m_DraftRect.minX)*(m_DraftRect.maxX-m_DraftRect.minX)
		+(m_DraftRect.maxY-m_DraftRect.minY)*(m_DraftRect.maxY-m_DraftRect.minY)
		+(m_DraftRect.maxZ-m_DraftRect.minZ)*(m_DraftRect.maxZ-m_DraftRect.minZ))/2;

		DistCam_default=RScene / sin(alpha/(2*rad));

		DistCam=DistCam_default;

		LookPos=LookCenter+DistCam*VectorCam;

		double x0,y0,z0;
		double x1,y1,z1;
		double u1,u2,u3;
		double v1,v2,v3;
		double A,B,C;
		double NormU, NormV;
		C3DPoint Log_0, Log0, LogU;

		Log_0.x=0;
		Log_0.y=0;
		Log_0.z=0;

		x1=LookPos.x;
		y1=LookPos.y;
		z1=LookPos.z;

		x0=LookCenter.x;
		y0=LookCenter.y;
		z0=LookCenter.z;

		A=x0-x1;
		B=y0-y1;
		C=z0-z1;

		Log0=ProjToSq(Log_0,A,B,C);
		LogU=ProjToSq(DirectionUp,A,B,C);

		u1=LogU.x-Log0.x;
		u2=LogU.y-Log0.y;
		u3=LogU.z-Log0.z;

		v1=u2*C-u3*B;
		v2=u3*A-u1*C;
		v3=u1*B-u2*A;

		NormU=sqrt(u1*u1+u2*u2+u3*u3);
		NormV=sqrt(v1*v1+v2*v2+v3*v3);

		u1/=NormU;
		u2/=NormU;
		u3/=NormU;

		v1/=-NormV;
		v2/=-NormV;
		v3/=-NormV;

		DirectionUp.x=u1;
		DirectionUp.y=u2;
		DirectionUp.z=u3;

		DirectionRight.x=v1;
		DirectionRight.y=v2;
		DirectionRight.z=v3;

		SetProjection();
		Invalidate();

		}else{

		}

		*/
	}
}
//----------------------------------