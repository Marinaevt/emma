#include "stdafx.h"
#include "LineStrip.h"


/////////////////////////////////////////////////////////////////////////////
//CLineStrip 
/////////////////////////////////////////////////////////////////////////////
CLineStrip::CLineStrip(bool bLoop) {
	m_bLoop = bLoop;
}
CLineStrip::~CLineStrip() {
}

/////////////////////////////////////////////////////////////////////////////
bool CLineStrip::Delete() {
	m_PointArray.FreeExtra();
	return true;
}

/////////////////////////////////////////////////////////////////////////////
CTypedPtrArray<CPtrArray, Math::C2DPoint*> *CLineStrip::GetArray() {
	return &m_PointArray;
}

/////////////////////////////////////////////////////////////////////////////
ptrdiff_t CLineStrip::AddNext(Math::C2DPoint *RealPoint) {
	return m_PointArray.Add(RealPoint);
}

//////////////////////////////////////////////////////////////////////////////
/*void CLineStrip::Draw(CDC *pDC, double Mult){
int i = m_PointArray.GetSize() - 1;
Math::C2DPoint *pRealPoint;

if(i >= 0){
if(m_bLoop){
pRealPoint = m_PointArray.GetAt(0);
}
else{
pRealPoint = m_PointArray.GetAt(i--);
if(i<0) i=0;
}
m_PointArray.GetAt(i)->Draw(pDC, Mult, *pRealPoint);
pRealPoint = m_PointArray.GetAt(i);
while(i--){
m_PointArray.GetAt(i)->Draw(pDC, Mult, *pRealPoint);
pRealPoint = m_PointArray.GetAt(i);
}
}
}*/

/////////////////////////////////////////////////////////////////////////////
bool CLineStrip::HaveInside(double x, double y) {
	ptrdiff_t Size = m_PointArray.GetSize(), i;
	Math::C2DPoint *pRealPoint1, *pRealPoint2;
	bool ret = false,
		F_ward = false;

	if (Size <= 0) return false;

	pRealPoint2 = m_PointArray[Size - 1];
	if (fabs(pRealPoint2->m_x - x) < EPS) {
		//double x0;
		i = Size - 1;
		while (i--) {
			//x0 = m_PointArray[i]->m_x;
			if (m_PointArray[i]->m_x > x) {
				F_ward = true;
				break;
			}
			else {
				if (m_PointArray[i]->m_x < x) {
					F_ward = false;
					break;
				}
			}

		}
	}
	for (i = 0;i<Size;i++) {
		pRealPoint1 = pRealPoint2;
		pRealPoint2 = m_PointArray.GetAt(i);

		if ((fabs(pRealPoint1->m_x - x) < EPS && fabs(pRealPoint1->m_y - y) < EPS) ||
			(fabs(pRealPoint2->m_x - x) < EPS && fabs(pRealPoint2->m_y - y) < EPS) ||
			(fabs(pRealPoint1->m_x - x) < EPS && fabs(pRealPoint2->m_x - x) < EPS &&
				((pRealPoint1->m_y > y && pRealPoint2->m_y < y) ||
					(pRealPoint1->m_y < y && pRealPoint2->m_y > y)))) {
			return true;
		}

		if (pRealPoint1->m_x < x && pRealPoint2->m_x > x) {
			if ((x - pRealPoint1->m_x)*(pRealPoint2->m_y - pRealPoint1->m_y)>
				(y - pRealPoint1->m_y)*(pRealPoint2->m_x - pRealPoint1->m_x)) {
				ret = !ret;
			}
			else {
				if (fabs(
					(x - pRealPoint1->m_x)*(pRealPoint2->m_y - pRealPoint1->m_y) -
					(y - pRealPoint1->m_y)*(pRealPoint2->m_x - pRealPoint1->m_x)
					) < EPS) {
					return true;
				}
			}
		}

		if (pRealPoint1->m_x > x && pRealPoint2->m_x < x) {
			if ((x - pRealPoint1->m_x)*(pRealPoint2->m_y - pRealPoint1->m_y)<
				(y - pRealPoint1->m_y)*(pRealPoint2->m_x - pRealPoint1->m_x)) {
				ret = !ret;
			}
			else {
				if (fabs(
					(x - pRealPoint1->m_x)*(pRealPoint2->m_y - pRealPoint1->m_y) -
					(y - pRealPoint1->m_y)*(pRealPoint2->m_x - pRealPoint1->m_x)
					) < EPS) {
					return true;
				}
			}
		}

		if (fabs(pRealPoint2->m_x - x) < EPS) {
			if (pRealPoint1->m_x < x) {
				F_ward = false;
			}
			else {
				if (pRealPoint1->m_x > x) {
					F_ward = true;
				}
			}
		}
		if ((fabs(pRealPoint1->m_x - x) < EPS) &&
			((pRealPoint2->m_x < x &&  F_ward) ||
				(pRealPoint2->m_x > x && !F_ward)) &&
			(pRealPoint1->m_y > y)) {
			ret = !ret;
		}
	}
	return ret;
}

/////////////////////////////////////////////////////////////////////////////
double CLineStrip::MinDistance(Math::C2DPoint *pDP, Math::C2DPoint *pDPRes) {
	ptrdiff_t Size = m_PointArray.GetSize(), i;
	Math::C2DPoint *pRealPoint1, *pRealPoint2;
	double not_res, res = 0.0,
		x = 0.0, y = 0.0;
	bool bNO_RES = true;

	i = Size;
	if (m_bLoop) {
		pRealPoint2 = m_PointArray.GetAt(0);
	}
	else {
		pRealPoint2 = m_PointArray.GetAt(i--);
		if (i<0) i = 0;
	}
	while (i--) {
		pRealPoint1 = pRealPoint2;
		pRealPoint2 = m_PointArray.GetAt(i);
		not_res = pDP->GetDistance(*pRealPoint1, *pRealPoint2, pDPRes);
		if (bNO_RES || not_res < res) {
			res = not_res;
			x = pDPRes->m_x;
			y = pDPRes->m_y;
			bNO_RES = false;
		}
	}
	pDPRes->m_x = x;
	pDPRes->m_y = y;
	return res;
}

/////////////////////////////////////////////////////////////////////////////
double CLineStrip::MinDistanceEx(Math::C2DPoint *pDP, Math::C2DPoint *pDPRes, int *index, double *Lambda) {
	ptrdiff_t Size = m_PointArray.GetSize(), i;
	Math::C2DPoint *pRealPoint1, *pRealPoint2;
	double not_res, res = 0.0,
		x = 0.0, y = 0.0;
	bool bNO_RES = true;

	i = Size;
	if (m_bLoop) {
		pRealPoint2 = m_PointArray.GetAt(0);
	}
	else {
		pRealPoint2 = m_PointArray.GetAt(i--);
		if (i<0) i = 0;
	}
	while (i--) {
		pRealPoint1 = pRealPoint2;
		pRealPoint2 = m_PointArray.GetAt(i);
		not_res = pDP->GetDistance(*pRealPoint1, *pRealPoint2, pDPRes);
		if (bNO_RES || not_res < res) {
			res = not_res;
			x = pDPRes->m_x;
			y = pDPRes->m_y;
			bNO_RES = false;
			*index = i;
		}
	}
	pDPRes->m_x = x;
	pDPRes->m_y = y;
	*Lambda = pDPRes->Splitting(*m_PointArray[*index], *m_PointArray[(*index + 1) % Size]);

	return res;
}

//////////////////////////////////////////////////////////////////////////////
bool CLineStrip::Belong(Math::C2DPoint *pDP) {
	ptrdiff_t i = m_PointArray.GetSize();

	if (i <= 0) return false;

	while (i--) {
		if (fabs(pDP->m_x - m_PointArray.GetAt(i)->m_x) < EPS &&
			fabs(pDP->m_y - m_PointArray.GetAt(i)->m_y) < EPS) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////
/*bool CLineStrip::Cross(Math::C2DPoint *pDP1, Math::C2DPoint *pDP2, FRONTINF inf1, FRONTINF inf2){
int i, Size  = m_PointArray.GetSize();
Math::C2DPoint *pDPfr1, *pDPfr2, DPRes;

pDPfr2 = m_PointArray[0];
i = Size;
while(i--){
pDPfr1 = pDPfr2;
pDPfr2 = m_PointArray[i];
if(	DPRes.Cross(pDP1, pDP2, pDPfr1, pDPfr2)
&&
inf1.rib != i && inf1.vertex != i && inf1.vertex != (i+1)%Size &&
inf2.rib != i && inf2.vertex != i && inf2.vertex != (i+1)%Size
){
return true;
}
}
return false;
}
*/
/////////////////////////////////////////////////////////////////////////////
