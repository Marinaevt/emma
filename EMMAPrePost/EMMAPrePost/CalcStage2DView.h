#pragma once
#include "gl2dview.h"
#include "..\..\Common\2DPoint.h"
#include "RegisteredMessages.h"
#include "2DCalcStage.h"
#include "ChildFrm.h"

//Представление для CalcStage (вывод расчётов)
class CCalcStage2DView : public CGL2DView
{
	//friend class C2DCalcStagePane;

	DECLARE_DYNCREATE(CCalcStage2DView);
public:
	CCalcStage2DView(void);
	~CCalcStage2DView(void);
	virtual void SetDocument(CEMMADoc *pDoc);
	//virtual void SetCentr(double x, double y);	
	GLint* RectToViewport(const CRect &rect);

DECLARE_MESSAGE_MAP()
	afx_msg LRESULT OnDraftRect(WPARAM wParam, LPARAM lParam);
	afx_msg void OnSize(UNINT nType, int cx, int cy);
	afx_msg void OnMouseMove(UNINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UNINT nFlags, CPoint point);

protected:
	static const int RESIZE_LEFT = 1;
	static const int RESIZE_RIGHT = 2;
	static const int RESIZE_TOP = 4;
	static const int RESIZE_BOTTOM = 8;

	ptrdiff_t m_CurrentCursor;
	int m_ResizeState;
	bool m_IsScalePicked;
	bool m_IsScaleResizeMode;

	bool IsInScale(const CPoint &point);					// Попадание координат мыши в отображние шкалы
	void UpdateResizeState(const CPoint &point);			// Прверка на попадание мыши в края отображения шкалы
	void SetCursor(const CPoint &point);					// Устанавливает курсор
	void MoveScale(const CPoint &point);					// Двигает шкалу мышкой
	void ResizeScale(const CPoint &point);					// Меняет размеры шкалы мышкой
public:
	afx_msg int OnSetCursor(CWnd* pWnd, UNINT nHitTest, UNINT message);
	afx_msg void OnLButtonUp(UNINT nFlags, CPoint point);
	afx_msg void On3D();
	afx_msg LRESULT OnActivate(WPARAM wParam, LPARAM lParam);
};
