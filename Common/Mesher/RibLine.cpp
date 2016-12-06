#include "stdafx.h"
#include "RibLine.h"
#include "Mesh.h"

/////////////////////////////////////////////////////////////////////////////
//CRibLine
CRibLine::CRibLine() {
	m_pFinal = nullptr;
	m_pList = nullptr;
	m_next = 0;
}
CRibLine::~CRibLine() {
	while (GetNext());
}
/////////////////////////////////////////////////////////////////////////////
bool CRibLine::AddRib(Math::C2DPoint *pDP1, Math::C2DPoint *pDP2, FRONTINF *pInf1, FRONTINF *pInf2, bool bDel)
{
	CRibList *pNewRib;
	CRibList *pTemp1 = m_pList, *pPrior = nullptr;
	CRibList *pDel = pTemp1, *pDelPr = pPrior;

	while (pDel != nullptr) {
		if (bDel &&
			(pDel->pDP1 == pDP2 && pDel->pDP2 == pDP1) ||
			(pDel->pDP1 == pDP1 && pDel->pDP2 == pDP2)) {

			if (pDel == m_pFinal) {
				m_pFinal = pDelPr;
			}
			if (pDelPr != nullptr) {
				pDelPr->pNext = pDel->pNext;
			}
			else {
				m_pList = pDel->pNext;
			}
			delete pDel;
			return false;
		}
		pDelPr = pDel;
		pDel = pDel->pNext;
	}

	pNewRib = new CRibList;
	pNewRib->pDP1 = pDP1;
	pNewRib->pDP2 = pDP2;
	pNewRib->pInf1 = pInf1;
	pNewRib->pInf2 = pInf2;

	static int i = 0;
	i++;
	bool bFr = true;
	if (m_pList) bFr = (m_pList->pInf1->vertex != -1 || m_pList->pInf1->rib != -1) || (m_pList->pInf2->vertex != -1 || m_pList->pInf2->rib != -1);
	if (i % 2 || bFr) {
		pNewRib->pNext = nullptr;
		if (m_pFinal == nullptr) {
			m_pFinal = m_pList = pNewRib;
		}
		else {
			m_pFinal->pNext = pNewRib;
			m_pFinal = pNewRib;
		}
	}
	else {
		pNewRib->pNext = m_pList;
		m_pList = pNewRib;
		if (m_pFinal == nullptr) {
			m_pFinal = m_pList = pNewRib;
		}
	}

	return true;
}
/////////////////////////////////////////////////////////////////////////////
bool CRibLine::GetNext() {
	if (m_pList == nullptr) {
		return false;
	}

	CRibList *pTemp = m_pList;
	m_pList = m_pList->pNext;
	delete(pTemp);
	
	if (m_pList != nullptr) {
		return true;
	}

	return false;
}
/////////////////////////////////////////////////////////////////////////////
Math::C2DPoint *CRibLine::GetR1() {
	if (m_pList == nullptr) return nullptr;
	return m_pList->pDP1;

}
/////////////////////////////////////////////////////////////////////////////
Math::C2DPoint *CRibLine::GetR2() {
	if (m_pList == nullptr) return nullptr;
	return m_pList->pDP2;
}
/////////////////////////////////////////////////////////////////////////////
FRONTINF *CRibLine::GetInf1() {
	if (m_pList == nullptr) return nullptr;
	return m_pList->pInf1;
}
/////////////////////////////////////////////////////////////////////////////
FRONTINF *CRibLine::GetInf2() {
	if (m_pList == nullptr) return nullptr;
	return m_pList->pInf2;
}

/////////////////////////////////////////////////////////////////////////////
bool CRibLine::IfVertArea(CFinEl *pFinEl, int *b_pDPk, double dRibLen)
{
	const double kff = 0.25;
	double  distance, temp1, distance1;
	CRibList *pList, *pRL;
	FRONTINF *pInf;
	Math::C2DPoint *pDPNearest = nullptr;
	Math::C2DPoint DP1, DP2, DPRes;

	distance = dRibLen;
	distance1 = dRibLen * 3;

	pList = pRL = m_pList;
	double d0;
	while (pList != nullptr) {
		DP1 = *pList->pDP1;
		DP2 = *pList->pDP2;

		d0 = pFinEl->GetPk().Len(DP1) + kff*DP1.GetDistance(pFinEl->GetPi(), pFinEl->GetPj(), &DPRes);
		temp1 = DP1.GetDistance(pFinEl->GetPi(), pFinEl->GetPj(), &DPRes);
		if (
			pList->pDP1 != pFinEl->m_Points[0] &&
			pList->pDP1 != pFinEl->m_Points[1] &&
			pList->pDP1 != pFinEl->m_Points[2] &&
			!DP1.Right(pFinEl->GetPi(), pFinEl->GetPj()) &&
			d0 < distance && temp1 < distance1)
		{
			distance1 = temp1;
			pDPNearest = pList->pDP1;
			pInf = pList->pInf1;
			pRL = pList;
		}

		d0 = pFinEl->GetPk().Len(DP2) + kff*DP2.GetDistance(pFinEl->GetPi(), pFinEl->GetPj(), &DPRes);
		temp1 = DP2.GetDistance(pFinEl->GetPi(), pFinEl->GetPj(), &DPRes);
		if (
			pList->pDP2 != pFinEl->m_Points[0] &&
			pList->pDP2 != pFinEl->m_Points[1] &&
			pList->pDP2 != pFinEl->m_Points[2] &&
			!DP2.Right(pFinEl->GetPi(), pFinEl->GetPj()) &&
			d0 < distance && temp1 < distance1)
		{
			distance1 = temp1;
			pDPNearest = pList->pDP2;
			pInf = pList->pInf2;
			pRL = pList;
		}
		pList = pList->pNext;
	}

	if (pDPNearest != nullptr) {
		*b_pDPk = 0;
		pFinEl->m_Points[2] = pDPNearest;
		pFinEl->m_inf[2] = pInf;
		return true;
	}
	return false;
}


/////////////////////////////////////////////////////////////////////////////
bool CRibLine::IfTooClose(CFinEl *pFinEl, int *b_pDPk) {
	if (!*b_pDPk)return true;
	CRibList *pList = m_pList;

	Math::C2DPoint DP1, DP2;
	Math::C2DPoint	DPRes, DPTemp;
	Math::C2DPoint	DPi = pFinEl->GetPi();
	Math::C2DPoint	DPj = pFinEl->GetPj();
	Math::C2DPoint	DPk = pFinEl->GetPk();
	double distance = DPi.Len(DPj)*0.25, temp;

	Math::C2DPoint *pDP1, *pDP2;
	FRONTINF *pInf1, *pInf2;

	Math::C2DPoint *pDP = nullptr;
	FRONTINF *pInf = nullptr;
	while (pList != nullptr) {

		if (pList->pDP1 != pFinEl->m_Points[0] || pList->pDP2 != pFinEl->m_Points[1]) {
			DP1 = *pList->pDP1;
			DP2 = *pList->pDP2;

			temp = DPk.GetDistance(DP1, DP2, &DPRes);
			if (temp < distance) {
				pDP1 = pList->pDP1;  pDP2 = pList->pDP2;
				pInf1 = pList->pInf1;   pInf2 = pList->pInf2;
				if (pFinEl->m_Points[0] == pDP2 || pFinEl->m_Points[1] == pDP2 || DP2.Right(DPi, DPj)) {
					if (!DP1.Right(DPi, DPj)) {
						pDP = pDP1;
						pInf = pInf1;
					}
				}
				else if (pFinEl->m_Points[0] == pDP1 || pFinEl->m_Points[1] == pDP1 || DP1.Right(DPi, DPj)) {
					if (!DP2.Right(DPi, DPj)) {
						pDP = pDP2;
						pInf = pInf2;
					}
				}
				else if (DPRes.Len(DP1) < DPRes.Len(DP2)) {
					pDP = pDP1;
					pInf = pInf1;
				}
				else {
					pDP = pDP2;
					pInf = pInf2;
				}
				if (pDP == pFinEl->m_Points[0] || pDP == pFinEl->m_Points[1]) {
					CDlgShowError cError(ID_ERROR_RIBLINE_PDP_NULL); //_T("pDP = nullptr;"));
					pDP = nullptr;
				}
			}
		}
		pList = pList->pNext;
	}

	if (pDP != nullptr) {
		*b_pDPk = 0;
		//distance = temp; не используется
		pFinEl->m_Points[2] = pDP;
		pFinEl->m_inf[2] = pInf;
	}

	return (pDP != nullptr);
}
/////////////////////////////////////////////////////////////////////////////
void CRibLine::Clear()
{
	//if(m_pList) 
	delete[] m_pList;
	//if(m_pFinal) 
	delete[] m_pFinal;
}
/////////////////////////////////////////////////////////////////////////////
void CRibLine::IfCross(CFinEl *pFinEl, int *b_pDPk) {
	Math::C2DPoint DPi = pFinEl->GetPi(),
		DPj = pFinEl->GetPj(),
		DPk = pFinEl->GetPk(),
		DP1, DP2, DPTemp;
	Math::C2DPoint *pDPi = pFinEl->m_Points[0],
		*pDPj = pFinEl->m_Points[1],
		*pDPk = pFinEl->m_Points[2],
		*pDP1, *pDP2;

	bool bNextStep = true;
	CRibList *pList = m_pList;

	//Проверка пересечения с фронтом построения
	while (pList != nullptr) {
		pDP1 = pList->pDP1;
		pDP2 = pList->pDP2;
		DP1 = *pDP1;
		DP2 = *pDP2;
		DPk = *pDPk;
		if (pDP1 != pDPk && pDP2 != pDPk) {
			if (pDP1 != pDPi && pDP2 != pDPi &&	DPTemp.Cross(DPi, DPk, DP1, DP2)) {
				*b_pDPk = 0;
				if ((pDPj == pDP2 || DP2.Right(DPi, DPj)) && !DP1.Right(DPi, DPj)) {
					pFinEl->m_Points[2] = pDPk = pDP1;
					pFinEl->m_inf[2] = pList->pInf1;
				}
				else {
					if ((pDPj == pDP1 || DP1.Right(DPi, DPj)) && !DP2.Right(DPi, DPj)) {
						pFinEl->m_Points[2] = pDPk = pDP2;
						pFinEl->m_inf[2] = pList->pInf2;
					}
					else {
						if (DP1.GetDistance(DPi, DPj, &DPTemp) < DP2.GetDistance(DPi, DPj, &DPTemp)) {
							pFinEl->m_Points[2] = pDPk = pDP1;
							pFinEl->m_inf[2] = pList->pInf1;
						}
						else {
							pFinEl->m_Points[2] = pDPk = pDP2;
							pFinEl->m_inf[2] = pList->pInf2;
						}
					}
				}
				bNextStep = false;
			}

			if (bNextStep && pDP1 != pDPj && pDP2 != pDPj && pDP1 != pDPk && pDP2 != pDPk && DPTemp.Cross(DPj, DPk, DP1, DP2)) {
				*b_pDPk = 0;
				if ((pDPi == pDP2 || DP2.Right(DPi, DPj)) && !DP1.Right(DPi, DPj)) {
					pFinEl->m_Points[2] = pDPk = pDP1;
					pFinEl->m_inf[2] = pList->pInf1;
				}
				else {
					if ((pDPi == pDP1 || DP1.Right(DPi, DPj)) && !DP2.Right(DPi, DPj)) {
						pFinEl->m_Points[2] = pDPk = pDP2;
						pFinEl->m_inf[2] = pList->pInf2;
					}
					else {
						if (DP1.GetDistance(DPi, DPj, &DPTemp) < DP2.GetDistance(DPi, DPj, &DPTemp)) {
							pFinEl->m_Points[2] = pDPk = pDP1;
							pFinEl->m_inf[2] = pList->pInf1;
						}
						else {
							pFinEl->m_Points[2] = pDPk = pDP2;
							pFinEl->m_inf[2] = pList->pInf2;
						}
					}
				}
				bNextStep = false;
			}
		}

		if (bNextStep) {//Если случилось пересечение, необходимо искать заново
			pList = pList->pNext;
		}
		else {
			pList = m_pList;
			bNextStep = true;
		}
	}


	//Проверка случая, когда участок фронта построения полностью оказывается внутри нового элемента
	pList = m_pList;
	while (pList != nullptr) {
		pDP1 = pList->pDP1; pDP2 = pList->pDP2;
		if (pDPi != pDP1 && pDPj != pDP1 && pDPk != pDP1 &&	pFinEl->HaveInside(pDP1->m_x, pDP1->m_y)) {
			*b_pDPk = 0;
			pFinEl->m_inf[2] = pList->pInf1;
			pFinEl->m_Points[2] = pDPk = pDP1;
		}
		if (pDPi != pDP2 && pDPj != pDP2 && pDPk != pDP2 &&	pFinEl->HaveInside(pDP2->m_x, pDP2->m_y)) {
			*b_pDPk = 0;
			pFinEl->m_inf[2] = pList->pInf2;
			pFinEl->m_Points[2] = pDPk = pDP2;
		}
		pList = pList->pNext;
	}
}

/////////////////////////////////////////////////////////////////////////////
void CRibLine::Offset(long long Poffset, long long Foffset, long long l_min, long long l_max) {
	CRibList *pList = m_pList;

	// используем (long long) или (pointer) для x64
	//Type           ILP64   LP64   LLP64
	//long             64     64      32
	//long long        64     64      64
	//pointer          64     64      64

	while (pList != nullptr) {
		if ((long long)pList->pDP1 >= l_min && (long long)pList->pDP1 < l_max) {
			pList->pDP1 = (Math::C2DPoint *)((long long)pList->pDP1 + Poffset);
			pList->pInf1 = (FRONTINF *)((long long)pList->pInf1 + Foffset);
		}

		if ((long long)pList->pDP2 >= l_min && (long long)pList->pDP2 < l_max) {
			pList->pDP2 = (Math::C2DPoint *)((long long)pList->pDP2 + Poffset);
			pList->pInf2 = (FRONTINF *)((long long)pList->pInf2 + Foffset);
		}
		pList = pList->pNext;
	}
}

//! Загрузка 
bool CRibLine::Load(C2DOutline &Outline, CMesh *pMesh, DBL len)
{
	if (Outline.GetContourCount() == 0) return false;

	for (size_t l = 0; l < Outline.GetContourCount(); l++)
	{
		for (size_t p = 0; p < Outline.GetCurveCount(); p++)
		{
			C2DCurve *pCur = Outline.GetCurve(p);
			Math::C2DPoint DP0 = pCur->GetStartNode()->GetPoint();
			Math::C2DPoint DP1 = pCur->GetEndNode()->GetPoint();

			Math::C2DPoint DP = DP0;
			Math::C2DPoint DPa;

			DBL l1 = pMesh->GetRibLen(DP);	//получаем длину ребра сетки
			DBL l0 = 0.0;
			FRONTINF FI = { -1, l, int(p) - 1 };	//p-1

			pMesh->AddPoint(pCur->GetStartNode()->GetPoint(), FI);

			while (l1 < pCur->CalcLength())
			{
				if (l1 + pMesh->GetRibLen(DP) > pCur->CalcLength()) {
					l1 = 0.5*(l0 + pCur->CalcLength());
				}
				pCur->GetPoint(l1, DP, DPa);
				FRONTINF infDP = { -1, l, int(p) - 1 };	//p-1
				pMesh->AddPoint(DP, infDP);

				l0 = l1;
				l1 += pMesh->GetRibLen(DP); //
			}
		}
	}

	bool bIfAddSuccess = true;
	for (size_t i = 1; i<pMesh->m_nPointsNum; i++)
	{
		bIfAddSuccess = bIfAddSuccess && AddRib(&pMesh->m_Points[i], &pMesh->m_Points[i - 1], &pMesh->m_FrInf[i], &pMesh->m_FrInf[i - 1], false);
	}
	bIfAddSuccess = bIfAddSuccess && AddRib(&pMesh->m_Points[0], &pMesh->m_Points[pMesh->m_nPointsNum - 1], &pMesh->m_FrInf[0], &pMesh->m_FrInf[pMesh->m_nPointsNum - 1], false);

	///////////////////////////////////////////////////////
	/*
	int nPN0 = pMesh->m_nPointsNum;
	for(int v = 0; v < Outline.GetNodeCount(); v++)
	{
	FRONTINF FI = {v, -1, -1};
	pMesh->AddPoint(Outline.GetNode(v)->GetPoint(), FI);
	}

	for(int l = 0; l < Outline.GetContourCount(); l++)
	{
	for(int p = 0; p <= Outline.GetCurveCount()-1; p++)
	{
	Math::C2DPoint DP0 = Outline.GetCurve(p)->GetStartNode()->GetPoint();
	Math::C2DPoint DP1 = Outline.GetCurve(p)->GetEndNode()->GetPoint();
	Math::C2DPoint DP = DP0;
	Math::C2DPoint DPa;
	FRONTINF infDP0 = {DP0, -1, -1};
	int nPNl = pMesh->m_nPointsNum;
	int nPos = nPN0 + Outline.GetPointPos(l, p);
	int len = 1 //Front.GetValue(DP, infDP);
	if(len < Outline.GetCurve(p)->CalcLength())
	{
	DP.Rad(DP0, DP1, len);
	FRONTINF infDP = {-1, l, p-1};
	while(len < Outline.GetCurve(p)->CalcLength())
	{
	DP.Rad(DP0, DP1, len);
	int nPos1 = pMesh->AddPoint(DP, infDP);
	AddRib(&pMesh->m_Points[nPos], &pMesh->m_Points[nPos1], &pMesh->m_FrInf[nPos], &pMesh->m_FrInf[nPos1], false);
	//if(!Outline.GetContour(l)->IsClosed())
	AddRib(&pMesh->m_Points[nPos1], &pMesh->m_Points[nPos], &pMesh->m_FrInf[nPos1], &pMesh->m_FrInf[nPos], false);
	nPos = nPos1;
	DP0 = DP;
	len += 10; //Front.GetValue(DP, infDP);
	if(len*0.5 >= Outline.GetCurve(p)->CalcLength())break;
	DP.Rad(DP0, DP1, len*0.5);
	}
	Math::C2DPoint DP_0 = Outline.GetCurve(p)->GetStartNode()->GetPoint();
	double dKff = Outline.GetCurve(p)->CalcLength()/(Outline.GetCurve(p)->CalcLength()+len);
	for(int i = nPNl; i < pMesh->m_nPointsNum; i++)
	{
	DP = pMesh->m_Points[i];
	DP = DP_0 + (DP - DP_0)*dKff;
	pMesh->m_Points[i] = DP;
	}
	int nPos1 = nPN0 + Outline.GetPointPos(l,p);
	AddRib(&pMesh->m_Points[nPos], &pMesh->m_Points[nPos1], &pMesh->m_FrInf[nPos], &pMesh->m_FrInf[nPos1], false);
	AddRib(&pMesh->m_Points[nPos1], &pMesh->m_Points[nPos], &pMesh->m_FrInf[nPos1], &pMesh->m_FrInf[nPos], false);
	}
	}
	}
	*/
	///////////////////////////////


	return bIfAddSuccess;
}
/////////////////////////////////////////////////////////////////////////////
//bool CRibLine::Load(CFrontier &Front, CMesh *pMesh){
//	Front.CheckWeights();
//	int nPN0 = pMesh->m_nPointsNum;
//	for(int v=0; v<Front.m_nPointsNum; v++){
//		FRONTINF FI = {v, -1, -1};
//		pMesh->AddPoint(Front.m_Points[v], FI);
//	}
//	for(int l = 0; l<Front.GetLoopsNum(); l++)
//	{
//		for(int p = 1; p<Front.GetPointsNum(l) + (Front.IsLoop(l) ? 1:0); p++)
//		{
//			int nPNl = pMesh->m_nPointsNum;
//			int nPos = nPN0 + Front.GetPointPos(l, p-1);
//			Math::C2DPoint DP0 = Front.GetPoint(l, p-1);
//			Math::C2DPoint DP1 = Front.GetPoint(l, p);
//			Math::C2DPoint DP;
//			FRONTINF infDP0 = {Front.GetPointPos(l, p-1), -1, -1};
//			double len = Front.GetValue(DP0, infDP0)*0.5;
//			if(len < DP0.Len(DP1)){
//				DP.Rad(DP0, DP1, len);
//				FRONTINF infDP = {-1, l, p-1};
//				len = Front.GetValue(DP, infDP);
//				while(len < DP0.Len(DP1)){
//					DP.Rad(DP0, DP1, len);
//					int nPos1 = pMesh->AddPoint(DP, infDP);
//					//AddRib(&pMesh->m_Points[nPos], &pMesh->m_Points[nPos1], &pMesh->m_FrInf[nPos], &pMesh->m_FrInf[nPos1], 1, !Front.IsFront(l)); 
//					AddRib(&pMesh->m_Points[nPos], &pMesh->m_Points[nPos1], &pMesh->m_FrInf[nPos], &pMesh->m_FrInf[nPos1], !Front.IsFront(l), false); 
//					if(!Front.IsFront(l))AddRib(&pMesh->m_Points[nPos1], &pMesh->m_Points[nPos], &pMesh->m_FrInf[nPos1], &pMesh->m_FrInf[nPos], !Front.IsFront(l), false); 
//					nPos = nPos1;
//
//					DP0 = DP;
//					len = Front.GetValue(DP, infDP);
//					if(len*0.5 >= DP0.Len(DP1))break;
//					DP.Rad(DP0, DP1, len*0.5);
//					len = Front.GetValue(DP, infDP);
//				}
//				Math::C2DPoint DP_0 = Front.GetPoint(l, p-1);
//				double dKff = DP_0.Len(DP1)/(DP_0.Len(DP0)+len);
//				for(int i=nPNl; i<pMesh->m_nPointsNum; i++){
//					DP = pMesh->m_Points[i];
//					DP = DP_0 + (DP - DP_0)*dKff; 
//					pMesh->m_Points[i] = DP;
//				}
//				int nPos1 = nPN0 + Front.GetPointPos(l, p);
//				AddRib(&pMesh->m_Points[nPos], &pMesh->m_Points[nPos1], &pMesh->m_FrInf[nPos], &pMesh->m_FrInf[nPos1], !Front.IsFront(l), false); 
//				if(!Front.IsFront(l))AddRib(&pMesh->m_Points[nPos1], &pMesh->m_Points[nPos], &pMesh->m_FrInf[nPos1], &pMesh->m_FrInf[nPos], !Front.IsFront(l), false); 
//			}
//		}
//	}
//	return true;
//}

/////////////////////////////////////////////////////////////////////////////

void CRibLine::Sort() {
	/*
	CRibList *pRL1, *pRL2;
	CRibList *pTemp1 = m_pList;

	while(pTemp1!=nullptr && pTemp1->pNext!=nullptr && !pTemp1->pNext->bFront){
	if( min((long)pTemp1->pDP1, (long)pTemp1->pDP2) <= min((long)pTemp1->pNext->pDP1, (long)pTemp1->pNext->pDP2)){
	pTemp1 = pTemp1->pNext;
	}else{
	pRL2 = pRL1 = pTemp1->pNext;
	while(pRL2->pNext!=nullptr && min((long)pRL2->pDP1, (long)pRL2->pDP2) == min((long)pRL2->pNext->pDP1, (long)pRL2->pNext->pDP2)){
	pRL2 = pRL2->pNext;
	}
	pTemp1->pNext = pRL2->pNext;
	pTemp1 = m_pList;
	while(pTemp1->pNext!=nullptr && !pTemp1->pNext->bFront && min((long)pTemp1->pNext->pDP1, (long)pTemp1->pNext->pDP2) < min((long)pRL2->pDP1, (long)pRL2->pDP2)){
	pTemp1 = pTemp1->pNext;
	}
	pRL2->pNext = pTemp1->pNext;
	pTemp1->pNext = pRL1;
	}
	}//*/
}

/////////////////////////////////////////////////////////////////////////////
void CRibLine::ChangePoint(Math::C2DPoint *pOld, Math::C2DPoint *pNew, FRONTINF *pNewInf) {
	CRibList *pTemp = m_pList;
	FRONTINF *pOldInf = nullptr;

	while (pTemp != nullptr) {
		if (pTemp->pDP1 == pOld) {
			pTemp->pDP1 = pNew;
			pOldInf = pTemp->pInf1;
			pTemp->pInf1 = pNewInf;
		}
		if (pTemp->pDP2 == pOld) {
			pTemp->pDP2 = pNew;
			pOldInf = pTemp->pInf2;
			pTemp->pInf2 = pNewInf;
		}
		pTemp = pTemp->pNext;
	}
	if (pOldInf) {
		delete pOld;
		delete pOldInf;
	}
	Sort();
}

/////////////////////////////////////////////////////////////////////////////
bool CRibLine::DivAng(Math::C2DPoint **Points, FRONTINF **Inf, int nSize, CFrontier &Front) {
	int n = 0;
	for (int p = 1; p<Front.m_nPointsNum; p++) {
		if (Front.m_Points[p].m_x < Front.m_Points[n].m_x)n = p;
	}
	//FRONTINF FI = {n,-1,-1}; //не используется
	//int nFr = Front.GetFrLoop(FI); //не используется

	for (int p = 0; p< Front.m_nPointsNum; p++) {
		FRONTINF FI = { p, -1, -1 };
		int nLoop = Front.GetFrLoop(FI);
		if (nLoop >= 0) {
			Math::C2DPoint *pDP = nullptr, *pDP1 = nullptr, *pDP2 = nullptr;
			FRONTINF *pFI1, *pFI2, *pFI;
			for (int i = 0; i<nSize; i++) {
				if (Inf[i * 2 + 1]->vertex == p) {
					pDP1 = Points[i * 2];
					pFI1 = Inf[i * 2];
					if (pDP && pDP != Points[i * 2 + 1])return false;
					pDP = Points[i * 2 + 1];
					pFI = Inf[i * 2 + 1];
				}
				if (Inf[i * 2]->vertex == p) {
					pDP2 = Points[i * 2 + 1];
					pFI2 = Inf[i * 2 + 1];
					if (pDP && pDP != Points[i * 2])return false;
					pDP = Points[i * 2];
					pFI = Inf[i * 2];
				}
				if (pDP1 != nullptr && pDP2 != nullptr)break;
			}
			if (pDP1 != nullptr && pDP2 != nullptr) {
				Math::C2DPoint DP1(*pDP1), DP2(*pDP2), DP(*pDP), DPR;
				if (!DP.Right(DP1, DP2)) continue;
				DP1.Rad(DP, DP1, 1);
				DP2.Rad(DP, DP2, 1);

				DPR = (DP1 + DP2)*0.5;
				DPR.Rad(DP, DPR, Front.GetValue(DP, FI));

				FRONTINF *pINF = new FRONTINF;
				pINF->nLoop = pINF->rib = pINF->vertex = -1;
				Math::C2DPoint *pDPk = new Math::C2DPoint;
				// 1234
				*pDPk = DPR;

				CFinEl FinEl;
				FinEl.m_Points[0] = pDP1;
				FinEl.m_Points[1] = pDP;
				FinEl.m_Points[2] = pDPk;
				FinEl.m_inf[0] = pFI1;
				FinEl.m_inf[1] = pFI;
				FinEl.m_inf[2] = pINF;

				int b_pDPk = 1;
				//	IfVertArea(&FinEl, &b_pDPk, &Front, Front.GetValue(DP, FI));
				//	IfTooClose(&FinEl, &b_pDPk, &Front);
				IfCross(&FinEl, &b_pDPk);

				FinEl.m_Points[0] = pDP;
				FinEl.m_Points[1] = pDP2;
				FinEl.m_inf[0] = pFI;
				FinEl.m_inf[1] = pFI2;

				//	IfVertArea(&FinEl, &b_pDPk, &Front, Front.GetValue(DP, FI));
				//	IfTooClose(&FinEl, &b_pDPk, &Front);
				IfCross(&FinEl, &b_pDPk);

				//TODO: просто удаляем?
				if (!b_pDPk) {
					delete pDPk;
					delete pINF;
				}
				AddRib(pDP, FinEl.m_Points[2], pFI, FinEl.m_inf[2], 1);
			}
		}
	}
	return true;
}

