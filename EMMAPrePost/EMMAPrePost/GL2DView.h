#pragma once
#include "Resource.h"
#include "GlBasicView.h"
#include "../../Common/2DPoint.h"
#include "../../Common/2DOutline.h"
#include "I2DSketch.h"

////////////////////////////////////////////////////////////////////
typedef struct STROKE_PARAM{
	int nType;	 // Тип штриха (1, 5, 10, 50, 100)
	int nWidth;  // Ширина штриха в пикселях
	int nLehgth; // Длина штриха в пикселях
	double ColorR[4]; // Цвет штриха (рулетка)
	double ColorG[4]; // Цвет штриха (сетка)
} STROKE_PARAM;

////////////////////////////////////////////////////////////////////
typedef struct AXIS_PARAM{
	double	dDelta; // Минимальный шаг
	double	dPCX;   // Величина пикселя
	int		nScaleWidth;   // Ширина полосы осей
	double	ColorRuler[4]; // Цвет рулетки
	double	ColorRulerB[4];// Цвет границы рулетки
	STROKE_PARAM Stroke100;
	STROKE_PARAM Stroke50;
	STROKE_PARAM Stroke10;
	STROKE_PARAM Stroke5;
	STROKE_PARAM Stroke1;
} AXIS_PARAM;

/**!
	2D представление для всех документов
*/
class CGL2DView : public CGLBasicView
{
public:
	CGL2DView(void);
	~CGL2DView(void);

public:
	void ShowAxis(bool bShow, bool bUpdate = true);			// Показать оси
	void ShowGrid(bool bShow, bool bUpdate = true);			// Показать сетку
	void ChangeCursor(int curType=0, bool bUpdate = true);	// Смена типа курсора
	int IsAxisVisible(){return m_bShowAxis;}
	int IsGridVisible(){return m_bShowGrid;}

	double GetZoom(){return m_dZoom;}
	void SetZoom(double dZoom, bool bUpdate = true);		// Установка масштаба
	virtual void SetCentr(double x, double y);				// Установка центра области чертежа
	Math::C2DPoint FizToLogical(CPoint FizPoint);			// Преобразование координат мыши
	Math::C2DPoint GetCentr(){return Math::C2DPoint(m_dCx, m_dCy);}
	void SetDraftRect(const CRect2D &rect, bool bCenter = false); // Устанавливает границы чертежа
	virtual bool ShowToolBar(UNINT nID){return nID == IDR_2DVIEW_TB;}
	bool Is3DViewShow; //
	bool View;

protected:
	enum{X_axis, Y_axis, XY_axis};
	enum{ST_axis, ST_grid};
	CRect2D m_DraftRect; // Координаты чертежа
	double m_dCx; // Координата центра X
	double m_dCy; // Координата центра Y
	double m_dZoom; // Увеличение
	CPoint m_pPrevPos; // Предыдущее положение мыши 
	bool	m_bShowAxis; // true если рисуем оси
	bool	m_bShowGrid; // true если рисуем сетку

	AXIS_PARAM m_AxisParam; // Параметры отображения осей	

protected:
	void GetStrokeParam(size_t n, STROKE_PARAM *pStrokeParam);
	void RecalcAxisParam(void); // Пересчет значений m_AxisParam.dDelta и m_AxisParam.dPCX
	virtual void DrawGrid(void);
	virtual void DrawAxis(void);
	virtual void DrawStroke(const STROKE_PARAM &StrokeParam, double dPos, int nDirection, int nStrokeType);
	virtual double DrawText(const STROKE_PARAM &StrokeParam, double dPos, int nDirection);
	void DrawRuling(size_t nDirection, size_t nStrokeType);
	virtual void SetProjection(void);	//Установки просмотра и проекции
	virtual void DrawScene(void); // Отрисовка всей сцены
	bool GetScreenRect(CRect2D &rect); // Возвращает координаты видимой области

public:
	DECLARE_MESSAGE_MAP()
	afx_msg int OnMouseWheel(UNINT nFlags, short zDelta, CPoint pt);
	afx_msg int OnSetCursor(CWnd* pWnd, UNINT nHitTest, UNINT message);	// Установка вида курсора
	afx_msg void OnLButtonDown(UNINT nFlags, CPoint point);				// Действие по левому клику на области чертежа
	afx_msg void OnRButtonDown(UNINT nFlags, CPoint point);				// Действие по правому клику на области чертежа
	afx_msg void OnMouseMove(UNINT nFlags, CPoint point);
	afx_msg LRESULT OnInitialUpdate(WPARAM wParam, LPARAM lParam);

	afx_msg void OnArrow();
	afx_msg void OnUpdateArrow(CCmdUI *pCmdUI);
	afx_msg void OnZoom();
	afx_msg void OnUpdateZoom(CCmdUI *pCmdUI);
	afx_msg void OnShowAll();
	afx_msg void OnUpdateShowAll(CCmdUI *pCmdUI);
	afx_msg void OnHand();
	afx_msg void OnUpdateHand(CCmdUI *pCmdUI);
	//Разделитель
	afx_msg void OnAxis();
	afx_msg void OnUpdateAxis(CCmdUI *pCmdUI);
	afx_msg void OnGrid();
	afx_msg void OnUpdateGrid(CCmdUI *pCmdUI);
	afx_msg void OnColor();
	afx_msg void OnUpdateColor(CCmdUI *pCmdUI);

	afx_msg void On3DView();  // для отображения 3D представления
	afx_msg void OnUpdate3DView(CCmdUI *pCmdUI);
};
