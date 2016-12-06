#pragma once
#include "EMMADoc.h"
#include "EMMAPrePost.h"

// CBasicView
#define	CBasicView_Parent CWnd
////////////////////////////////////////////////////////////////////////////////
class CBasicView : public CBasicView_Parent
{
public:
////////////////////////////////////////////////////////////////////////////////
	CBasicView();// Конструктор
	virtual ~CBasicView();// Деструктор
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

public:
	CEMMADoc *GetDocument(){return m_pDoc;}
	GLParam *GetParameter(){return m_Parameter;}
	virtual void SetDocument(CEMMADoc *pDoc);
	HTREEITEM GetTreeItem(){return m_pDoc ? m_pDoc->GetTreeItem() : nullptr;}
	virtual void DrawDoc(GLParam &parameter);
	virtual bool ShowToolBar(UNINT nID){return 0;}

protected:
	CEMMADoc *m_pDoc;
	GLParam *m_Parameter;									// Параметр отрисовки

protected:
	bool m_bLButtonPressed;// истина, если нажата левая кнопка мыши
	bool m_bRButtonPressed;// истина, если нажата правая кнопка мыши
	bool m_bMButtonPressed;// истина, если нажата средняя кнопка мыши
	bool m_bMouseIsInside;// истина, если курсор мыши находится внутри окна
	CPoint m_MousePos;// текущие координаты мыши


protected:
	afx_msg void OnLButtonDown(UNINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UNINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UNINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UNINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UNINT nFlags, CPoint point);
	afx_msg int OnMouseWheel(UNINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnMouseLeave();
	afx_msg void OnMButtonDown(UNINT nFlags, CPoint point);
	afx_msg void OnMButtonUp(UNINT nFlags, CPoint point);
	void TrackMouseEvent(UNLONG dwFlags);
	////////////////////////////////

	int GetCurrentMouseX(){return m_MousePos.x;};
	int GetCurrentMouseY(){return m_MousePos.y;};
	//__declspec(property(get=GetCurrentMouseX)) int CurrentMouseX;
	//__declspec(property(get=GetCurrentMouseY)) int CurrentMouseY;

	////////////////////////////////
	DECLARE_MESSAGE_MAP()
public:
	afx_msg int OnEraseBkgnd(CDC* pDC);
	afx_msg void OnDestroy();
	
};
////////////////////////////////////////////////////////////////////////////////


