// OpenGLControl.cpp: Implementierungsdatei
//
#include "stdafx.h"
#include "OpenGLControl.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

/////////////////////////////////////////////////////////////////////////////
// COpenGLControl

COpenGLControl::COpenGLControl()
{
	dc = nullptr;
	m_Scale = nullptr;
}

COpenGLControl::~COpenGLControl()
{
	delete dc;
}

void COpenGLControl::SetScale(C2DScale *scale)
{
	m_Scale = scale;
}


BEGIN_MESSAGE_MAP(COpenGLControl, CWnd)
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_CREATE()
	ON_WM_ERASEBKGND()
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// Behandlungsroutinen fьr Nachrichten COpenGLControl 


void COpenGLControl::InitGL()
{
	glShadeModel(GL_SMOOTH);
	glClearColor(1, 1, 1, 1);
	glClearDepth(1.0f);							
	glEnable(GL_DEPTH_TEST);					
	glDepthFunc(GL_LEQUAL);	
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
}

void COpenGLControl::DrawGLScene()
{
	glClear(GL_COLOR_BUFFER_BIT
		 |  GL_DEPTH_BUFFER_BIT);


	glLoadIdentity();

	//***************************
	// DRAWING CODE
	//***************************

	if (!m_Scale)
		return;
	
	size_t scaleCount = m_Scale->GetCount();

	for (size_t i = 0; i < scaleCount; i++)
	{
		C2DColor color = m_Scale->GetColor(m_Scale->GetValue(i));
		glColor3d(color.GetR(), color.GetG(), color.GetB());

		double step = 1.9 / scaleCount;	// TODO: Волшебные числа

		glBegin(GL_QUADS);
			glVertex2d(-0.95 + i * step, -0.3);
			glVertex2d(-0.95 + i * step + step, -0.3);
			glVertex2d(-0.95 + i * step + step, 0.3);
			glVertex2d(-0.95 + i * step, 0.3);
		glEnd();

		glColor3d(0,0, 0);
		glBegin(GL_LINES);
			glVertex2d(-0.95 + i * step, -0.4);
			glVertex2d(-0.95 + i * step, -0.6);
		glEnd();
	}

	glColor3d(0,0, 0);

	glBegin(GL_LINES);
		glVertex2d(0.95, -0.4);
		glVertex2d(0.95, -0.6);
	glEnd();

	glBegin(GL_LINES);
		glVertex2d(-0.95, -0.6);
		glVertex2d(0.95, -0.6);
	glEnd();

	SwapBuffers(dc->m_hDC);
}



void COpenGLControl::Create(CRect rect, CWnd *parent)
{
	CString className = AfxRegisterWndClass(
		CS_HREDRAW | CS_VREDRAW | CS_OWNDC,
		nullptr,
		(HBRUSH)GetStockObject(BLACK_BRUSH),
		nullptr);


	CreateEx(
		0,
		className,
		_T("OpenGL"),
		WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_CLIPCHILDREN,
		rect,
		parent,
		0);

	m_Rect = rect;
}

void COpenGLControl::OnPaint() 
{
	openGLDevice.makeCurrent();

	DrawGLScene();
}

void COpenGLControl::OnSize(UNINT nType, int cx, int cy) 
{
	CWnd::OnSize(nType, cx, cy);
	
	if (cy == 0)								
	{
		cy = 1;						
	}


	glViewport(0,0,cx,cy);	

	glMatrixMode(GL_PROJECTION);						
	glLoadIdentity();						

	
	glOrtho(-1.0f,1.0f,-1.0f,1.0f,1.0f,-1.0f);

	glMatrixMode(GL_MODELVIEW);						
	glLoadIdentity();
}


int COpenGLControl::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	dc = new CClientDC(this);

	openGLDevice.create(dc->m_hDC);
	InitGL();

	return 0;
}

int COpenGLControl::OnEraseBkgnd(CDC* pDC) 
{
	//COLORREF crBkgnd = pDC->GetBkColor();
	COLORREF crDialog = ::GetSysColor(COLOR_3DFACE);
	GLclampf gl_red = float(((crDialog << 8) >> 24) / 255.0);
	GLclampf gl_green = float(((crDialog << 16) >> 24) / 255.0);
	GLclampf gl_blue = float(((crDialog << 24) >> 24) / 255.0);
	glClearColor(gl_red, gl_green, gl_blue, 1);

	return 1;
}
