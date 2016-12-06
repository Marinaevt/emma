#pragma once
#include "Globals.h"
#include "../2DOutline.h"

/////////////////////////////////////////////////////////////////////////////

class CRibList{
public:
	bool bFront;
	Math::C2DPoint *pDP1;
	Math::C2DPoint *pDP2;
	FRONTINF *pInf1;
	FRONTINF *pInf2;
	CRibList *pNext;
	double Distance;
};

/////////////////////////////////////////////////////////////////////////////
class CMesh;

/////////////////////////////////////////////////////////////////////////////
class CRibLine{
public:
	CRibList *m_pList;
	CRibList *m_pFinal;	
//	Math::C2DPoint m_Center;
	int m_next;
public:
	CRibLine();
	~CRibLine();

public:
	void ChangePoint(Math::C2DPoint *pOld, Math::C2DPoint *pNew, FRONTINF *pNewInf);
	bool Load(C2DOutline &Outline, CMesh *pMesh, DBL len);
	void Offset(long long Poffset, long long Foffset, long long l_min, long long l_max);

	bool AddRib(Math::C2DPoint *pDP1, Math::C2DPoint *pDP2, FRONTINF *Inf1, FRONTINF *Inf2, bool bDel = true);
	bool GetNext();
	Math::C2DPoint *GetR1();
	Math::C2DPoint *GetR2();
	FRONTINF *GetInf1();
	FRONTINF *GetInf2();
	bool IfVertArea( CFinEl *pFinEl, int *b_pDPk, double dRibLen);
	bool IfTooClose( CFinEl *pFinEl, int *b_pDPk);
	void IfCross( CFinEl *pFinEl, int *b_pDPk);
	void Clear();
protected:
	void Sort();
	bool DivAng(Math::C2DPoint **Points, FRONTINF **Inf, int nSize, CFrontier &Front);
};
/////////////////////////////////////////////////////////////////////////////
