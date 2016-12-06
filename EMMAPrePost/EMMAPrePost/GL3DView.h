#pragma once
#include "Resource.h"
#include "GlBasicView.h"
#include "../../Common/3DPoint.h"
#include "../../Common/2DPoint.h"
#include "Results.h"

////////////////////////////////////////////////////////////////////
class CGL3DView : public CGLBasicView
{
public:
	CGL3DView(void);
	~CGL3DView(void);

public:
	void ShowAxis(bool bShow, bool bUpdate = true); 
	void ShowGrid(bool bShow, bool bUpdate = true);
	void ShowPerspective(bool bShow, bool bUpdate = true); 
	void ChangeCursor(int curType=0, bool bUpdate = true); // поменять тип курсора. 0 - стрелка, 1 - лупа, 2 - рука, 4 - поворот
	bool IsAxisVisible(){return m_bShowAxis;}
	bool IsGridVisible(){return m_bShowGrid;}
	bool IsPerspective(){return m_bShowPerspective;}

	double GetZoom(){return m_dZoom;}
	void SetZoom(double dZoom, bool bUpdate = true);
	CPoint3D FizToLogical3D(CPoint FizPoint);
	CPoint3D ProjToSq(CPoint3D point, double A, double B, double C);
	double ToUgol(CPoint3D p1, CPoint3D p2); // Переводит изменение логических координат (соответствующих положению мыши) в угол поворота сцены
	double ScalPr(CPoint3D a1, CPoint3D a2){return a1.x*a2.x+a1.y*a2.y+a1.z*a2.z;}; // скалярное произведение двух векторов
	void LimitArea(); // Ограничение обзора сцены

	CRect3D m_DraftRect; // Координаты чертежа
	double m_dZoom; //Zoom

	double DistCam; // Расстояние от камеры до её фокуса
	double DistCam_default; // Расстояние от камеры до её фокуса по умолчанию
	double RScene; // Радиус описанной вокруг сцены сферы

	bool Is2DViewShow;  //
	bool IsCentr3DShow; //
	bool IsLightOn;     //
	
//свежее
	void SetDraftRect(double m_X, double mX, double m_Y, double mY, double m_Z, double mZ);
	void SetDraftRect(const CRect3D &rect, bool bCenter = false); // Устанавливает границы чертежа

	CPoint3D VectorCam; // Единичный направляющий вектор от фокуса до камеры
	CPoint3D LookCenter; // Координаты фокуса камеры
	CPoint3D LookPos; // Координаты положения камеры  (LookPos=LookCenter+VectorCam*DistCam)
	
	double alpha;  // Угол обзора камеры в перспективной проекции вдоль наименьшей оси (в градусах) 
	//(y - если высота области отображения больше ширины, x - если ширина больше высоты)

	CPoint3D CenterScene; // Координаты центра сцены
	CPoint3D DirectionUp; // Координаты вектора направления вверх
	CPoint3D DirectionRight; // Координаты вектора направления вправо

	CPoint m_pPrevPos; // Предыдущее физическое положение мыши 
	CPoint3D m_pPrevPosSp; // Предыдущее логическое положение мыши
	bool m_bShowAxis; // true если рисуем оси
	bool m_bShowGrid; // true если рисуем сетку
	bool m_bShowPerspective; // true если отображаем с перспективой

	bool IsWPressed; // 
	bool IsAPressed; // 
	bool IsSPressed; // 
	bool IsDPressed; //
	bool IsShiftPressed; // 
	bool IsLeftPressed; // 
	bool IsRightPressed; // 
	virtual bool ShowToolBar(UNINT nID){return nID == IDR_3DVIEW_TB;}

	bool m_ShowLight3Dobj;

protected:
	virtual void DrawGrid(void);
	virtual void DrawAxis(void);
	virtual void SetProjection(void);
	virtual void DrawScene(void);// Отрисовка всей сцены
	bool GetScreenRect(CRect2D &rect);// Возвращает координаты видимой области

public:
	DECLARE_MESSAGE_MAP()
	afx_msg int OnMouseWheel(UNINT nFlags, short zDelta, CPoint pt);
	afx_msg int OnSetCursor(CWnd* pWnd, UNINT nHitTest, UNINT message);
	afx_msg void OnLButtonDown(UNINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UNINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UNINT nFlags, CPoint point);


	afx_msg void OnArrow();
	afx_msg void OnUpdateArrow(CCmdUI *pCmdUI);
	afx_msg void OnZoom();
	afx_msg void OnUpdateZoom(CCmdUI *pCmdUI);
	afx_msg void OnHand();
	afx_msg void OnUpdateHand(CCmdUI *pCmdUI);
	afx_msg void OnAxis();
	afx_msg void OnUpdateAxis(CCmdUI *pCmdUI);
	afx_msg void OnGrid();
	afx_msg void OnUpdateGrid(CCmdUI *pCmdUI);
	afx_msg void OnColor();
	afx_msg void OnUpdateColor(CCmdUI *pCmdUI);

	afx_msg void OnLight();
	afx_msg void OnUpdateLight(CCmdUI *pCmdUI);
	afx_msg void OnPerspective();
	afx_msg void OnUpdatePerspective(CCmdUI *pCmdUI);

	afx_msg void OnFly();
	afx_msg void OnUpdateFly(CCmdUI *pCmdUI);
	afx_msg void OnRotateFree();
	afx_msg void OnUpdateRotateFree(CCmdUI *pCmdUI);

	afx_msg void OnKeyDown(UNINT nChar, UNINT nRepCnt, UNINT nFlags);
	afx_msg void OnKeyUp(UNINT nChar, UNINT nRepCnt, UNINT nFlags);
	
	afx_msg void OnTimer(size_t nIDEvent);

	afx_msg LRESULT OnInitialUpdate(WPARAM wParam, LPARAM lParam);
	afx_msg void OnShowAll();
	afx_msg void OnUpdateShowAll(CCmdUI *pCmdUI);

	afx_msg void On2DView();  // для отображения 2D представления
	afx_msg void OnUpdate2DView(CCmdUI *pCmdUI);

	afx_msg void OnCentr3D();  // для отображения центра 3D представления
	afx_msg void OnUpdateCentr3D(CCmdUI *pCmdUI);
};
