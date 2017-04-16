#pragma once
#include "Resource.h"
#include "2DSketch.h"
#include "2DTrajectory.h"
#include "../../Common/2DOutline.h"


//! Класс Инструментов
class C2DTool :	public C2DSketch
{
	GETTYPE (C2DTOOL)
public:
	
	C2DTool(CString strName = _T(""));
	virtual ~C2DTool(void);

	virtual bool IsPaneVisible(size_t nID) /*inl*/ {return nID == IDR_2DTOOL_PANE;}

	virtual void ChangeNodes(Math::C2DPoint *new_point);	// Перемещение чертежа Инструмента
	C2DTrajectory* GetTrajectory();	// Возвращает ссылку на Траекторию
	bool CopyOutline(C2DTool *tool);
	
	virtual void AddStickPoint(Math::C2DPoint ppoint, GLParam &parameter);
	void ResetCountourCache();
	
	CMFCPropertyGridProperty* AddProp(CStringW name, double eval, size_t pos);
	void FillPropList(CMFCPropertyGridCtrl *pGrid);	// Заполнение параметров (правая панель)
	//void UpdateProp(double *pEval, UNLONG *pId);		// Обновление параметров (правая панель)
	void UpdatePropList(CMFCPropertyGridCtrl *pGrid);	//обновляет данные документа (Doc) по новым значениям таблицы
	
	virtual void DrawGL(GLParam &parameter);
	bool GetBoundingBox(CRect2D &rect);

public:
	IO::C2DPoint m_StickPoint;	//Точка прилипания Траектории
	IO::CDouble m_FrictionCoeff; //Коэффициент трения (>0.5 - заделка, <0.5 - симметрия)
	
};

IODEF (C2DTool)