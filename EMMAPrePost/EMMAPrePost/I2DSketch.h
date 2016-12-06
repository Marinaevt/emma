#pragma once
#include "../../Common/2DMotion.h"

//! Интерфейс, который должен реализовать классы документов, представления которых унаследованы от ChildView
class I2DSketch
{
public:
	I2DSketch(void);
	~I2DSketch(void);
	virtual C2DOutline *GetOutline() = 0;

	virtual void SetOutline(C2DOutline *) = 0;
	virtual void SetContours() = 0;
	virtual void SetClearOutline() = 0;
	void ClearActiveObjects();	//очищаем все активные элементы

	//Делаем активной узел\кривую\дугу
	virtual void SetActiveObject(Math::C2DPoint *ppoint);
	void I2DSketch::SetTempLine(Math::C2DPoint *ppoint);
	void I2DSketch::SetTempArc(Math::C2DPoint *ppoint);
	virtual void SetTempFacet(Math::C2DPoint *ppoint, bool set);

	void CopyInM_OutlineLine(C2DFacet_2lines *out);
	void CopyInM_OutlineLineArc(C2DFacet_line_arc *out);
	void CopyInM_OutlineArcs(C2DFacet_2arcs *out);

	void CopyInOutlineLine(int c);
	void CopyInOutlineArc(int c);
	void CopyInOutline(int nC1,int nC2, int nP, short t1, short t2);
	//End Yura
	
	virtual CMFCPropertyGridProperty * AddProp(CStringW, double, size_t) = 0;

	virtual void DrawGL(GLParam &parameter);
	virtual void OnClick(int *cur);
	//virtual void OnClick();


	virtual void AddItem(Math::C2DPoint *ppoint, int cur);

	virtual void AddNode(Math::C2DPoint *ppoint);
	int AddFirstNode(Math::C2DPoint *ppoint);
	virtual void AddCurve(Math::C2DPoint *ppoint);
	virtual void AddCircleArc(Math::C2DPoint *ppoint);

	virtual void ChangeTime(GLParam &parameter);
	void SetSpeedOfTime(DBL dNew);
	bool GetBoundingBox(CRect2D &rect);

protected:
	//Возвращает номер узла и записывает в nCurve1 и nCurve2 номера кривых, на которых можно построить фаску
	//если курсор находится в точке pPoint
	//В случае, если фаску построить невозможно, возвращает -1;
	int FindFacetNode(Math::C2DPoint *pPoint, int &nCurve1, int &nCurve2);

protected:
	
	std::vector<int> m_CurvesActive;	//массив индексов активных кривых
	std::vector<int> m_NodesActive;		//массив индексов активных узлов

	int m_NumCurve;
	int m_NumNode;
	bool m_fAdd; //
	
	C2DOutline m_Outline;

public:
	bool m_fCtrl; //для Child2DView
	double m_dCurTime;
	double m_dTimeSpeed;
	bool m_fPlay; //для Child2DView и C2DTrajectory
};
