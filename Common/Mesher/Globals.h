#pragma once
#include "../2DPoint.h"
#include "../IOInterface.h"
#include "../GLParam.h"

/////////////////////////////////////////////////////////////////////////////
typedef struct FRONTINF {
	int vertex;//точка границы, к которой прилипла вершина. 
			   //(-1 если таковая отсутствует)

	int nLoop;//номер контура к которому принадлежит вершина,
			  //имеет смысл ТОЛЬКО если vertex == -1
			  //(-1 если таковая отсутствует)
	int rib;//номер ребра контура nLoop, которому принадлежит вершина,
			//имеет смысл ТОЛЬКО если vertex == -1 и nLoop >= 0
			//(-1 если таковая отсутствует)
}FRONTINF;


/////////////////////////////////////////////////////////////////////////////
class CLineLoop {
public:
	Math::C2DPoint **m_Points;
	int m_nSize;
	bool	m_bLoop;
	bool	m_bFront;

public:
	CLineLoop();
	CLineLoop(const CLineLoop &LL);
	~CLineLoop();

	int		Add(Math::C2DPoint *DP);
	void	Invert();

	const CLineLoop& operator = (const CLineLoop &LL);
};

/////////////////////////////////////////////////////////////////////////////
class CFrontier
{
protected:
	Math::C2DPoint *m_FreePoints;
	double	*m_FPWeights;
	int m_nFreePointsNum;

public:
	Math::C2DPoint *m_Points;
	double	*m_Weights;
	int	*m_bRib;
	int m_nPointsNum;

	CLineLoop *m_Loops;
	int m_nLoopsNum;

	double	m_TgAng;
	double	m_dCeiling;
	int		m_nType;
	double	m_dVal;

public:
	CFrontier();
	virtual ~CFrontier();

	bool SetDirection(Math::C2DPoint **Points, FRONTINF **Inf, int nSize);
	int GetFrLoop(const FRONTINF &Inf);
	double GetValue(const Math::C2DPoint &DP, const FRONTINF &FI);
	//int Divide(Math::C2DPoint ***Points, FRONTINF ***Inf, bool **Dang);	//нигде не используется

	int GetLoopsNum() { return m_nLoopsNum; }
	int GetPointsNum(int nLoop) { return m_Loops[nLoop].m_nSize; }

	int AddPoint(const Math::C2DPoint &DBLP, double dWeight = -1);
	int AddLoop(const CLineLoop &LL);
	int AddFreePoint(const Math::C2DPoint &DBLP, double dWeight);
	void CheckWeights();

	Math::C2DPoint	GetPoint(int nLoop, int nRib);
	int		GetPointPos(int nLoop, int nRib);
	bool	IsFront(int nLoop);
	bool	IsLoop(int nLoop);

	void	Draw(CDC *pDC, double Mult);
	int		HaveInside(double x, double y, int nLoop);// возвращает -1 если точка снаружи, 0 - если на контуре, 1 - если внутри
	bool	HaveInside(double x, double y);
	double	MinDistanceEx(const Math::C2DPoint &DP, Math::C2DPoint *pDPRes, int *pnLoop, int *pnRib);
	bool	Cross(const Math::C2DPoint &DP1, const Math::C2DPoint &DP2, const FRONTINF &inf1, const FRONTINF &inf2, Math::C2DPoint *pDPRes = nullptr, FRONTINF *infRes = nullptr);
	double	GetWeight(int index);
	void	FindLeftEl(int *i_left, int *i_right, FRONTINF *pFInf);
	bool	SameLine(const FRONTINF &FI1, const FRONTINF &FI2, FRONTINF *pRes = nullptr);
	void	CheckDirection();// Проверка и исправление направления обхода контуров
	int		FindCoveringLoop();// находим внешний контур
	void	CheckLoopDirection(int nLoop, bool bCW);// Проверка и исправление направления обхода контура l в соответствии с bCW

};


//! Класс конечного элемента
class CFinEl {
public:
	FRONTINF *m_inf[3];
	int     m_status[3];
	int     m_mark[3];
	Math::C2DPoint	*m_Points[3];

public:
	CFinEl();
	~CFinEl();
	bool HaveRib(Math::C2DPoint *pDP1, Math::C2DPoint *pDP2);
	bool FrontCrossControl(CFrontier *pFront, int *b_pDPk);
	bool HaveInside(double x, double y);
	void DrawGl(GLParam &parameter);
	Math::C2DPoint GetPi() { return Math::C2DPoint(*m_Points[0]); }
	Math::C2DPoint GetPj() { return Math::C2DPoint(*m_Points[1]); }
	Math::C2DPoint GetPk() { return Math::C2DPoint(*m_Points[2]); }
};
