#pragma once
#include "../2DPoint.h"


/////////////////////////////////////////////////////////////////////////////
class CLineStrip
{
protected:
	CTypedPtrArray<CPtrArray,Math::C2DPoint*> m_PointArray;
	bool m_bLoop;
public:
	CLineStrip(bool bLoop = false);
	virtual ~CLineStrip();

	ptrdiff_t AddNext(Math::C2DPoint *pRealPoint);
//	void Draw(CDC *pDC, double Mult);
	bool HaveInside(double x, double y);
	double MinDistance(Math::C2DPoint *pDP, Math::C2DPoint *pDPRes);
	double MinDistanceEx(Math::C2DPoint *pDP, Math::C2DPoint *pDPRes, int *i, double *Lambda);
	CTypedPtrArray<CPtrArray,Math::C2DPoint*> *GetArray();
	bool Belong(Math::C2DPoint *pDP);
//	bool Cross( Math::C2DPoint *pDP1, Math::C2DPoint *pDP2,
//				FRONTINF inf1, FRONTINF inf2 );

	bool Delete();
};
