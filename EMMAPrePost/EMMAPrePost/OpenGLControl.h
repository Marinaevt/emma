#pragma once
#include "OpenGLDevice.h"
#include "DlgNewScale.h"
#include <gl\GL.h>
#include <gl\GLU.h>

/////////////////////////////////////////////////////////////////////////////
// Fenster COpenGLControl 

class COpenGLControl : public CWnd
{
// Konstruktion
public:
	COpenGLControl();

// Attribute
public:

// Operationen
public:

// Ьberschreibungen
	// Vom Klassen-Assistenten generierte virtuelle Funktionsьberschreibungen
	//{{AFX_VIRTUAL(COpenGLControl)
	public:
	//}}AFX_VIRTUAL

// Implementierung
public:
	void Create(CRect rect,CWnd* parent);
	virtual ~COpenGLControl();
	void SetScale(C2DScale *scale);




	// Generierte Nachrichtenzuordnungsfunktionen
protected:
	void InitGL();
	void DrawGLScene();
	OpenGLDevice openGLDevice;
	CClientDC* dc;
	C2DScale *m_Scale;
	CRect m_Rect;

	//{{AFX_MSG(COpenGLControl)
	afx_msg void OnPaint();
	afx_msg void OnSize(UNINT nType, int cx, int cy);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg int OnEraseBkgnd(CDC* pDC);
	//}}AFX_MSG

	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

