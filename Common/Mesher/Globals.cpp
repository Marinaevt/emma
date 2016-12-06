#include "stdafx.h"
#include "Globals.h" 

/////////////////////////////////////////////////////////////////////////////
//CLineLoop 
/////////////////////////////////////////////////////////////////////////////
CLineLoop::CLineLoop() {
	m_nSize = 0;
	m_Points = nullptr;
	m_bLoop = true;
	m_bFront = true;
}

const CLineLoop& CLineLoop::operator=(const CLineLoop &LL) {
	m_nSize = LL.m_nSize;
	if (m_Points) {
		delete[] m_Points;
		m_Points = nullptr;
	}
	if (m_nSize>0) {
		m_Points = new Math::C2DPoint *[LL.m_nSize];
		memcpy(m_Points, LL.m_Points, m_nSize);
	}
	//Копируем флаги
	m_bLoop = LL.m_bLoop;
	m_bFront = LL.m_bFront;
	return *this;
}

CLineLoop::CLineLoop(const CLineLoop &LL) {
	m_nSize = LL.m_nSize;
	if (m_Points) {
		delete[] m_Points;
		m_Points = nullptr;
	}
	if (m_nSize>0) {
		m_Points = new Math::C2DPoint *[LL.m_nSize];
		memcpy(m_Points, LL.m_Points, m_nSize);
	}
	//Копируем флаги
	m_bLoop = LL.m_bLoop;
	m_bFront = LL.m_bFront;
}

CLineLoop::~CLineLoop() {
	delete[] m_Points;
}

int CLineLoop::Add(Math::C2DPoint *pDP) {
	Math::C2DPoint **ppOld = m_Points;
	++m_nSize;
	m_Points = new Math::C2DPoint *[m_nSize];
	if (ppOld) {
		m_Points = (Math::C2DPoint **)memcpy((void *)m_Points, (void *)ppOld, sizeof(Math::C2DPoint *)*(m_nSize - 1));
		delete[] ppOld;
	}
	m_Points[m_nSize - 1] = pDP;
	return m_nSize;
}

void CLineLoop::Invert() {
	for (int i = 0; i<m_nSize / 2; i++) {
		Math::C2DPoint *tmp = m_Points[m_nSize - 1 - i];
		m_Points[m_nSize - 1 - i] = m_Points[i];
		m_Points[i] = tmp;
	}
}

/////////////////////////////////////////////////////////////////////////////
//CFrotier 
/////////////////////////////////////////////////////////////////////////////
CFrontier::CFrontier() {
	m_nLoopsNum = 0;
	m_nPointsNum = 0;
	m_nFreePointsNum = 0;
	m_Points = nullptr;
	m_Weights = nullptr;
	m_bRib = nullptr;
	m_FreePoints = nullptr;
	m_FPWeights = nullptr;
	m_Loops = nullptr;

	m_dCeiling = 1;
	m_nType = 0;
	m_TgAng = 0.0;
	m_dVal = 0.0;
}

CFrontier::~CFrontier() {
	//if(m_Points!=nullptr)
	delete[] m_Points;
	delete[] m_Weights;
	delete[] m_bRib;
	delete[] m_FreePoints;
	delete[] m_FPWeights;
	delete[] m_Loops;
}

/////////////////////////////////////////////////////////////////////////////
int CFrontier::AddPoint(const Math::C2DPoint &DP, double dWeight) {
	Math::C2DPoint *pOld = m_Points;
	double *pOldW = m_Weights;
	int *pOldR = m_bRib;
	++m_nPointsNum;
	m_Points = new Math::C2DPoint[m_nPointsNum];
	m_Weights = new double[m_nPointsNum];
	m_bRib = new int[m_nPointsNum];

	// используем (long long) или (pointer) для x64
	//Type           ILP64   LP64   LLP64
	//long             64     64      32
	//long long        64     64      64
	//pointer          64     64      64

	long long offset = (long long)m_Points - (long long)pOld;
	if (pOld) {
		m_Points = (Math::C2DPoint *)memcpy(m_Points, pOld, sizeof(Math::C2DPoint)*(m_nPointsNum - 1));
		m_Weights = (double *)memcpy(m_Weights, pOldW, sizeof(double)*(m_nPointsNum - 1));
		m_bRib = (int *)memcpy(m_bRib, pOldR, sizeof(int)*(m_nPointsNum - 1));
		delete[] pOld;
		delete[] pOldW;
		delete[] pOldR;
	}
	m_Points[m_nPointsNum - 1] = DP;
	m_Weights[m_nPointsNum - 1] = dWeight;
	m_bRib[m_nPointsNum - 1] = (dWeight>0) ? 1 : 0;

	for (int l = 0; l<m_nLoopsNum; l++) {
		for (int i = 0; i<m_Loops[l].m_nSize; i++) {
			m_Loops[l].m_Points[i] = (Math::C2DPoint *)((long long)m_Loops[l].m_Points[i] + offset);
		}
	}

	return m_nPointsNum;
}
/////////////////////////////////////////////////////////////////////////////
int CFrontier::AddFreePoint(const Math::C2DPoint &DP, double dWeight) {
	Math::C2DPoint *pOldPoints = m_FreePoints;
	double *pOldWeights = m_FPWeights;

	++m_nFreePointsNum;
	m_FreePoints = new Math::C2DPoint[m_nFreePointsNum];
	m_FPWeights = new double[m_nFreePointsNum];

	if (pOldPoints) {
		m_FreePoints = (Math::C2DPoint *)memcpy(m_FreePoints, pOldPoints, sizeof(Math::C2DPoint)*(m_nFreePointsNum - 1));
		m_FPWeights = (double *)memcpy(m_FPWeights, pOldWeights, sizeof(double)*(m_nFreePointsNum - 1));
		delete[] pOldPoints;
		delete[] pOldWeights;
	}
	m_FreePoints[m_nFreePointsNum - 1] = DP;
	m_FPWeights[m_nFreePointsNum - 1] = dWeight;

	return m_nFreePointsNum;
}

/////////////////////////////////////////////////////////////////////////////
int CFrontier::AddLoop(const CLineLoop &LL) {
	CLineLoop *pOld = m_Loops;
	++m_nLoopsNum;
	m_Loops = new CLineLoop[m_nLoopsNum];
	if (pOld) {
		m_Loops = (CLineLoop *)memcpy(m_Loops, pOld, sizeof(CLineLoop)*(m_nLoopsNum - 1));
		for (int i = 0; i<m_nLoopsNum - 1; i++) {
			pOld[i].m_Points = nullptr;
			pOld[i].m_nSize = 0;
		}
		delete[] pOld;
	}
	memcpy(&m_Loops[m_nLoopsNum - 1], &LL, sizeof(CLineLoop));
	m_Loops[m_nLoopsNum - 1].m_Points = new Math::C2DPoint *[LL.m_nSize];
	for (int i = 0; i<LL.m_nSize; i++) {
		m_Loops[m_nLoopsNum - 1].m_Points[i] = LL.m_Points[i];
	}
	return m_nLoopsNum;
}

//////////////////////////////////////////////////////////////////////////////
Math::C2DPoint CFrontier::GetPoint(int nLoop, int nRib) {

	if (m_Loops[nLoop].m_bLoop) {
		nRib = (nRib + m_Loops[nLoop].m_nSize) % m_Loops[nLoop].m_nSize;
	}
	else {
		nRib = max(0, min(nRib, m_Loops[nLoop].m_nSize - 1));
	}
	return	*m_Loops[nLoop].m_Points[nRib];
}

//////////////////////////////////////////////////////////////////////////////
int CFrontier::GetPointPos(int nLoop, int nRib) {
	if (nLoop<0 || nLoop >= m_nLoopsNum) {
		return 0;
	}
	if (m_Loops[nLoop].m_bLoop) {
		nRib = (nRib + m_Loops[nLoop].m_nSize) % m_Loops[nLoop].m_nSize;
	}
	else {
		nRib = max(0, min(nRib, m_Loops[nLoop].m_nSize - 1));
	}

	return	((long long)m_Loops[nLoop].m_Points[nRib] - (long long)m_Points) / sizeof(Math::C2DPoint);
}
//////////////////////////////////////////////////////////////////////////////
bool CFrontier::IsFront(int nLoop) {
	if (nLoop<0 || nLoop >= m_nLoopsNum) {
		return false;
	}
	return m_Loops[nLoop].m_bFront;
}

//////////////////////////////////////////////////////////////////////////////
bool CFrontier::IsLoop(int nLoop) {
	if (nLoop<0 || nLoop >= m_nLoopsNum) {
		return false;
	}
	return m_Loops[nLoop].m_bLoop;
}

//////////////////////////////////////////////////////////////////////////////
void CFrontier::Draw(CDC *pDC, double Mult)
{

	CPen pen(PS_SOLID, 1, RGB(0, 150, 0));
	CPen *pPen = pDC->SelectObject(&pen);

	int l = m_nLoopsNum;

	while (l--)
	{
		int r = m_Loops[l].m_nSize;
		if (!m_Loops[l].m_bLoop) r--;
		if (r>0)
		{
			Math::C2DPoint RealPoint = GetPoint(l, r);
			while (r--)
			{
				///////////////////////////////////////
				//выпилить
				//GetPoint(l,r).Draw(pDC, Mult, RealPoint);
				///////////////////////////////////////
				RealPoint = GetPoint(l, r);
			}
		}
	}
	pDC->SelectObject(pPen);
}//*/

 /*
 Проверяет расположение точки относительно контура
 Возвращает -1 если точка снаружи, 0 - если на контуре, 1 - если внутри
 */
int CFrontier::HaveInside(double x, double y, int nLoop) {

	if (!m_Loops[nLoop].m_bLoop) return -1;

	int nTop = 0, nBottom = 0;
	Math::C2DPoint DP1 = GetPoint(nLoop, 0), DP(x, y);
	bool bLeft = DP1.m_x < x ? true : false;
	int r = m_Loops[nLoop].m_nSize;

	while (r--)
	{
		if (DP == DP1)return 0;
		Math::C2DPoint DP2 = DP1;
		DP1 = GetPoint(nLoop, r);
		if ((DP1.m_x < x && !bLeft) || (DP1.m_x >= x && bLeft))
		{
			bLeft = !bLeft;
			//double y0 = (DP1.m_y-DP2.m_y)/(DP1.m_x-DP2.m_x)*(x-DP2.m_x) + DP2.m_y;
			//double dVal = (DP1.m_y-DP2.m_y)*(x-DP2.m_x)*(DP1.m_x-DP2.m_x) > (DP1.m_x-DP2.m_x)*(DP1.m_x-DP2.m_x)*(y - DP2.m_y);
			double dVal = (DP1.m_y - DP2.m_y)*(x - DP2.m_x) - (DP1.m_x - DP2.m_x)*(y - DP2.m_y);
			if (dVal > 0) {
				nTop++;
			}
			else {
				if (dVal < 0) {
					nBottom++;
				}
				else {
					return 0;
				}
			}
		}
		else if (fabs(DP1.m_x - x) < EPS && fabs(DP2.m_x - x) < EPS && (y - DP1.m_y)*(y - DP2.m_y) <= 0) {
			return 0;
		}
	}
	return 2 * (nTop % 2) - 1;
}


/*
Проверяет расположение точки относительно ВСЕХ контуров

*/
bool CFrontier::HaveInside(double x, double y)
{
	int N = 0;
	for (int l = 0; l<m_nLoopsNum; l++)
	{
		if (m_Loops[l].m_bFront && m_Loops[l].m_bLoop)
		{
			int nInside = HaveInside(x, y, l);
			if (nInside == 1) {
				N++;
			}
			else {
				if (nInside == 0) {
					return true;
				}
			}
		}
	}
	if (N % 2 != 0) {
		return true;
	}
	return false;
}



/////////////////////////////////////////////////////////////////////////////
/*double CFrontier::MinDistance(Math::C2DPoint *pDP, Math::C2DPoint *pDPRes){
Math::C2DPoint *pRealPoint1, *pRealPoint2;
double not_res, res, x, y;
bool NO_RES = true;

pRealPoint2 = m_PointArray.GetAt(Size-1);
for(int i=0;i<Size;i++){
pRealPoint1 = pRealPoint2;
pRealPoint2 = m_PointArray.GetAt(i);
not_res = pDP->GetDistance(pRealPoint1, pRealPoint2, pDPRes);
if(NO_RES || not_res < res){
res = not_res;
x=pDPRes->m_x;
y=pDPRes->m_y;
NO_RES = false;
}
}
pDPRes->m_x = x;
pDPRes->m_y = y;
return res;
}//*/

/////////////////////////////////////////////////////////////////////////////
double CFrontier::MinDistanceEx(const Math::C2DPoint &DP, Math::C2DPoint *pDPRes, int *pnLoop, int *pnRib) {

	Math::C2DPoint DP1, DP2;
	double temp, res = 0, x = 0, y = 0;
	bool bFirst = true;
	int nLoop = *pnLoop, nRib = *pnRib;

	for (int l = 0; l<m_nLoopsNum; l++) {
		//fprintf(f,"  l = %d \r\n", l);
		int r = m_Loops[l].m_nSize;
		if (!m_Loops[l].m_bLoop)r--;
		DP1 = GetPoint(l, r);
		while (r--) {
			DP2 = DP1;
			DP1 = GetPoint(l, r);
			if (l != nLoop || r != nRib) {
				temp = DP.GetDistance(DP1, DP2, pDPRes);
				if (temp < res || bFirst) {
					bFirst = false;
					res = temp;
					x = pDPRes->m_x;
					y = pDPRes->m_y;
					*pnLoop = l;
					*pnRib = r;
				}
			}
		}
	}
	pDPRes->m_x = x;
	pDPRes->m_y = y;
	return res;
}

//////////////////////////////////////////////////////////////////////////////
bool CFrontier::Cross(const Math::C2DPoint &DP1, const Math::C2DPoint &DP2, const FRONTINF &inf1, const FRONTINF &inf2, Math::C2DPoint *pDPRes, FRONTINF *infRes) {
	Math::C2DPoint p1, p2, DPRes;
	double len, tmp;
	bool bCross = false;

	for (int l = 0; l<m_nLoopsNum; l++) {
		int r = m_Loops[l].m_nSize;
		if (!m_Loops[l].m_bLoop)r--;
		p1 = GetPoint(l, r);
		while (r--) {
			FRONTINF Inf = { -1, l, r };
			p2 = p1;
			p1 = GetPoint(l, r);
			if (!SameLine(inf1, Inf) && !SameLine(inf2, Inf) && DPRes.Cross(p1, p2, DP1, DP2)) {
				tmp = DP1.Len(DPRes);
				if (!bCross || tmp < len) {
					bCross = true;
					len = tmp;
					if (pDPRes) {
						*pDPRes = DPRes;
					}
					if (infRes) {
						*infRes = Inf;
					}
				}
			}
		}
	}

	return bCross;
}

/////////////////////////////////////////////////////////////////////////////
double CFrontier::GetWeight(int index) {
	double temp, distance = 0, dLambda;
	bool bFirst = true;

	Math::C2DPoint DP1, DP2, DPRes, Targ = m_Points[index];
	FRONTINF inf1 = { -1,-1,-1 }, inf2 = { -1,-1,-1 }, infTarg = { -1,-1,-1 }, infRes = { -1,-1,-1 };

	infTarg.vertex = index;
	infTarg.rib = inf1.rib = inf2.rib = infRes.vertex = -1;
	infTarg.nLoop = inf1.nLoop = inf2.nLoop = -1;

	for (int l = 0; l<m_nLoopsNum; l++) {

		int r = m_Loops[l].m_nSize;
		if (!m_Loops[l].m_bLoop) r--;
		DP1 = GetPoint(l, r);

		while (r--) {
			DP2 = DP1;
			DP1 = GetPoint(l, r);
			int i = GetPointPos(l, r);
			int i1 = GetPointPos(l, r + 1);
			if (i != index && i1 != index) {
				inf1.vertex = i;
				inf2.vertex = i1;
				infRes.nLoop = l;
				infRes.rib = r;
				infRes.vertex = -1;
				temp = Targ.GetDistance(DP1, DP2, &DPRes);
				if ((temp < distance || bFirst) && HaveInside((DPRes.m_x + Targ.m_x)*0.5, (DPRes.m_y + Targ.m_y)*0.5)) {
					dLambda = DPRes.Splitting(DP1, DP2);
					if (fabs(dLambda) < EPS) {
						infRes = inf1;
					}
					else {
						if (fabs(dLambda - 1) < EPS) {
							infRes = inf2;
						}
					}
					if (!Cross(DPRes, Targ, infRes, infTarg)) {
						distance = temp;
						bFirst = false;
					}
				}
			}
			else {
				temp = DP1.Len(DP2);
				if (temp < distance || bFirst) {
					distance = temp;
					bFirst = false;
				}
			}
		}
	}
	return distance;
}

/////////////////////////////////////////////////////////////////////////////
void CFrontier::FindLeftEl(int *i_left, int *i_right, FRONTINF *pFInf) {
	int i;
	Math::C2DPoint DPTemp, DPRes, DP;
	double minY = 0;
	bool bFirst = true;

	i = *i_left = m_nPointsNum - 1;

	while (i--) {
		if (m_Points[i].m_x < m_Points[*i_left].m_x) {
			*i_left = i;
		}
	}
	for (int l = 0; l<m_nLoopsNum; l++) {
		if (m_Loops[l].m_bFront && m_Loops[l].m_bLoop) {
			for (int r = 0; r<m_Loops[l].m_nSize; r++) {
				if (GetPointPos(l, r) == *i_left) {
					int r1 = (r + m_Loops[l].m_nSize - 1) % m_Loops[l].m_nSize;
					int r2 = (r + 1) % m_Loops[l].m_nSize;
					DPRes.Rad(GetPoint(l, r), GetPoint(l, r1), 1);
					if (DPRes.m_y<minY || bFirst) {
						minY = DPRes.m_y;
						*i_right = GetPointPos(l, r1);
						pFInf->nLoop = l;
						pFInf->rib = r1;
						pFInf->vertex = -1;
					}
					DPRes.Rad(GetPoint(l, r), GetPoint(l, r2), 1);
					if (DPRes.m_y<minY) {
						minY = DPRes.m_y;
						*i_right = GetPointPos(l, r2);
						pFInf->nLoop = l;
						pFInf->rib = r;
						pFInf->vertex = -1;
					}
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
bool CFrontier::SameLine(const FRONTINF &FI1, const FRONTINF &FI2, FRONTINF *pRes) {
	if ((FI1.vertex == -1 && FI1.nLoop == -1) || (FI2.vertex == -1 && FI2.nLoop == -1))return false;

	if (pRes) {
		if (FI1.nLoop != -1) *pRes = FI1;
		else* pRes = FI2;
	}

	if (FI1.nLoop != -1 && FI1.nLoop == FI2.nLoop && FI1.rib == FI2.rib) return true;

	if (FI1.nLoop != -1 && (GetPointPos(FI1.nLoop, FI1.rib) == FI2.vertex ||
		GetPointPos(FI1.nLoop, FI1.rib + 1) == FI2.vertex)) return true;

	if (FI2.nLoop != -1 && (GetPointPos(FI2.nLoop, FI2.rib) == FI1.vertex ||
		GetPointPos(FI2.nLoop, FI2.rib + 1) == FI1.vertex)) return true;

	if (FI2.vertex != -1 && FI1.vertex != -1) {
		for (int l = 0; l<m_nLoopsNum; l++) {
			int r = m_Loops[l].m_nSize;
			while (r--) {
				if (FI1.vertex == GetPointPos(l, r)) {
					break;
				}
			}

			if (r != -1) {
				if (FI2.vertex == GetPointPos(l, r - 1)) {
					if (pRes) {
						pRes->vertex = -1;
						pRes->nLoop = l;
						pRes->rib = (r != 0) ? r - 1 : m_Loops[l].m_nSize - 1;
					}
					return true;
				}
				if (FI2.vertex == GetPointPos(l, r + 1)) {
					if (pRes) {
						pRes->vertex = -1;
						pRes->nLoop = l;
						pRes->rib = r;
					}
					return true;
				}
			}
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////
// Проверка и исправление направления обхода контуров
void CFrontier::CheckDirection() {
	int nCoveringLoop = FindCoveringLoop();

	if (nCoveringLoop < 0) {
		CDlgShowError cError(ID_ERROR_GLOBALS_NCOVERINGLOOP); //_T("nCoveringLoop < 0"));
		return;
	}

	CheckLoopDirection(nCoveringLoop, 0);// внешний контур должен быть против часовой стрелки
	for (int nLoop = 0; nLoop<m_nLoopsNum; nLoop++) {
		if (nLoop != nCoveringLoop) CheckLoopDirection(nLoop, 1);
	}
}

int CFrontier::FindCoveringLoop() {// находим внешний контур

	for (int l = 0; l<m_nLoopsNum; l++) {
		if (m_Loops[l].m_bFront) {//проверяем только граничные контуры
			
			// граничный контур должен быть замкнутым
			if (m_Loops[l].m_bLoop == false) {
				CDlgShowError cError(ID_ERROR_GLOBALS_MBLOOP_FALSE); //_T("m_bLoop is false"));
				return 0;
			}

			bool bCovering = true;
			for (int i = 0; i<m_nPointsNum; i++) {
				if (HaveInside(m_Points[i].m_x, m_Points[i].m_y, l) == -1) {
					bCovering = false;
					break;
				}
			}
			if (bCovering) return l;
		}
	}
	return -1;
}

void CFrontier::CheckLoopDirection(int nLoop, bool bCW) {// Проверка и исправление направления обхода контура l в соответствии с bCW
	if (!m_Loops[nLoop].m_bFront) return;//проверяем только граничные контуры
	
	// граничный контур должен быть замкнутым
	if (m_Loops[nLoop].m_bLoop == false) {
		CDlgShowError cError(ID_ERROR_GLOBALS_MBLOOP_FALSE);
		return;
	}

	int nLeft = 0;
	for (int i = 1; i<m_Loops[nLoop].m_nSize; i++) {
		if (GetPoint(nLoop, i).m_x < GetPoint(nLoop, nLeft).m_x)nLeft = i;
	}
	Math::C2DPoint LeftPoint = GetPoint(nLoop, nLeft);
	Math::C2DPoint PrevPoint = GetPoint(nLoop, nLeft - 1);
	Math::C2DPoint NextPoint = GetPoint(nLoop, nLeft + 1);
	LeftPoint.m_x -= 1;//Смещаем левую точку еще левее, на случай, если все три лежат на одной вертикальной прямой
	bool bActualCW = NextPoint.Right(PrevPoint, LeftPoint);
	if (bCW != bActualCW) {
		m_Loops[nLoop].Invert();
	}
}

/////////////////////////////////////////////////////////////////////////////
void CFrontier::CheckWeights() {
	for (int i = 0; i<m_nPointsNum; i++) {
		double minW = GetWeight(i)*1.5;
		if (m_Weights[i] <= 0 || m_Weights[i] > minW) m_Weights[i] = minW;
	}
}


/////////////////////////////////////////////////////////////////////////////
/* // нигде не используется 
int CFrontier::Divide(Math::C2DPoint ***Points, FRONTINF ***Inf, bool **Dang) {
	int nSize = 0;
	CheckWeights();
	for (int l = 0; l<m_nLoopsNum; l++) {
		for (int p = 1; p<m_Loops[l].m_nSize + (m_Loops[l].m_bLoop ? 1 : 0); p++) {	//оператор ?: был без скобок
			Math::C2DPoint DP0 = GetPoint(l, p - 1);
			Math::C2DPoint DP1 = GetPoint(l, p);
			Math::C2DPoint DP;
			FRONTINF infDP0 = { p - 1, -1, -1 };
			double len = GetValue(DP0, infDP0)*0.5;
			if (len < DP0.Len(DP1)) {
				DP.Rad(DP0, DP1, len);
				FRONTINF infDP = { -1, l, p - 1 };
				len = GetValue(DP, infDP);
				while (len < DP0.Len(DP1)) {
					nSize++;
					DP.Rad(DP0, DP1, len);
					DP0 = DP;
					len = GetValue(DP, infDP)*0.5;
					if (len >= DP0.Len(DP1))break;
					DP.Rad(DP0, DP1, len);
					len = GetValue(DP, infDP);
				}
			}
			nSize++;
		}
	}

	*Points = new Math::C2DPoint *[nSize * 2];
	*Inf = new FRONTINF *[nSize * 2];
	*Dang = new bool[nSize];
	int nR = 0;

	for (int l = 0; l<m_nLoopsNum; l++) {
		for (int p = 1; p<m_Loops[l].m_nSize + (m_Loops[l].m_bLoop ? 1 : 0); p++) {	//оператор ?: был без скобок
			int nR0 = nR;
			Math::C2DPoint DP0 = GetPoint(l, p - 1);
			Math::C2DPoint DP1 = GetPoint(l, p);
			Math::C2DPoint DP;
			FRONTINF infDP0 = { GetPointPos(l, p - 1), -1, -1 };
			(*Points)[nR * 2] = &m_Points[GetPointPos(l, p - 1)];
			FRONTINF *pFI = new FRONTINF;
			*pFI = infDP0;
			(*Inf)[nR * 2] = pFI;
			double len = GetValue(DP0, infDP0)*0.5;
			if (len < DP0.Len(DP1)) {
				DP.Rad(DP0, DP1, len);
				FRONTINF infDP = { -1, l, p - 1 };
				len = GetValue(DP, infDP);
				while (len < DP0.Len(DP1)) {

					DP.Rad(DP0, DP1, len);
					Math::C2DPoint *pDP = new Math::C2DPoint;
					pFI = new FRONTINF;
					*pFI = infDP;
					//наверное можно проще 1234
					*pDP = DP;
					(*Points)[nR * 2 + 1] = (*Points)[nR * 2 + 2] = pDP;
					(*Inf)[nR * 2 + 1] = (*Inf)[nR * 2 + 2] = pFI;
					(*Dang)[nR] = !m_Loops[l].m_bFront;
					nR++;

					DP0 = DP;
					len = GetValue(DP, infDP)*0.5;
					if (len >= DP0.Len(DP1))break;
					DP.Rad(DP0, DP1, len);
					len = GetValue(DP, infDP);
				}
			}
			(*Points)[nR * 2 + 1] = &m_Points[GetPointPos(l, p)];
			pFI = new FRONTINF;
			FRONTINF FInf = { GetPointPos(l, p), -1, -1 };
			*pFI = FInf;
			(*Inf)[nR * 2 + 1] = pFI;
			(*Dang)[nR] = !m_Loops[l].m_bFront;
			nR++;

			if (nR - nR0 > 1) {
				for (int i = 0; i<4; i++) {	// счётчик i нигде не используется!
					for (int n = nR0; n<nR - 1; n++) {
						Math::C2DPoint DPA(*(*Points)[2 * n]);
						Math::C2DPoint DPB(*(*Points)[2 * n + 1]);
						Math::C2DPoint DPC(*(*Points)[2 * n + 3]);
						
						FRONTINF FI = { -1, l, p - 1 };
						double AB, AC;
						AB = GetValue((DPA + DPB)*0.5, FI);
						AC = AB + GetValue((DPB + DPC)*0.5, FI);
						AB /= AC;
						DPB.Rad(DPA, DPC, AB*DPA.Len(DPC));
						*(*Points)[2 * n + 1] = DPB;
					}
				}
			}

		}
	}

	for (int i = 0; i<m_nPointsNum; i++) {
		Math::C2DPoint *pDP = new Math::C2DPoint;
		FRONTINF *pFI = new FRONTINF;
		bool bSet = true;
		for (int nP = 0; nP<nSize * 2; nP++) {
			if ((*Points)[nP] == &m_Points[i]) {
				if (bSet) {
					*pDP = m_Points[i];
					*pFI = *(*Inf)[nP];
					bSet = false;
				}
				delete (*Inf)[nP];
				(*Points)[nP] = pDP;
				(*Inf)[nP] = pFI;
			}
		}
	}

	if (!SetDirection(*Points, *Inf, nSize)) {
		if (nSize>0) {
			for (int i = 0; i<nSize * 2; i++) {
				if ((*Points)[i] != nullptr) {
					Math::C2DPoint *pDP = (*Points)[i];
					for (int j = i + 1;j<nSize * 2;j++) {
						if ((*Points)[j] == pDP) (*Points)[j] = nullptr;
					}
					delete (*Points)[i];
					delete (*Inf)[i];
				}
			}
			delete[](*Points);
			delete[](*Inf);
			delete[](*Dang);
		}
		nSize = 0;
	}

	return nSize;
}
// нигде не используется */ 

/////////////////////////////////////////////////////////////////////////////
double CFrontier::GetValue(const Math::C2DPoint &DP, const FRONTINF &FI) {
	if (m_nType == 1)return m_dCeiling;

	double dVal = m_dCeiling, dTemp, dDistance;
	double dW1;
	Math::C2DPoint DP1, DPRes;

	int i = m_nFreePointsNum;
	while (i--) {
		DP1 = m_FreePoints[i];
		dW1 = m_FPWeights[i];
		dDistance = DP.Len(DP1);
		dTemp = dW1 + dDistance*m_TgAng;
		if (dTemp < dVal)	dVal = dTemp;
	}

	for (int l = 0; l<m_nLoopsNum; l++) {
		for (int p = 1; p<m_Loops[l].m_nSize + 1; p++) {
			if (p<(m_Loops[l].m_nSize + (m_Loops[l].m_bLoop ? 1 : 0)) && m_bRib[GetPointPos(l, p - 1)] && m_bRib[GetPointPos(l, p)]) {
				Math::C2DPoint DP_0 = GetPoint(l, p - 1);
				Math::C2DPoint DP_1 = GetPoint(l, p);

				//расстояние не должно быть близко к 0.
				double dLenDP0DP1 = DP_0.Len(DP_1);
				double a = 0.0;
				if (dLenDP0DP1 > EPS) {
					a = DP_0.Len(DP1) / dLenDP0DP1;
				}

				dW1 = m_Weights[GetPointPos(l, p - 1)] * (1 - a) + m_Weights[GetPointPos(l, p)] * a;
			}
			else {
				DP1 = GetPoint(l, p - 1);
				dW1 = m_Weights[GetPointPos(l, p - 1)];
			}
			dDistance = DP.Len(DP1);
			dTemp = dW1 + dDistance*m_TgAng;
			if (dTemp < dVal)	dVal = dTemp;

			if (dVal <= 0) {
				CDlgShowError cError(ID_ERROR_GLOBALS_DVAL); //_T("dVal <= 0"));
				return 0.0;
			}

		}
	}
	return dVal;
}


/////////////////////////////////////////////////////////////////////////////
int CFrontier::GetFrLoop(const FRONTINF &Inf) {
	if (Inf.nLoop >= 0) {
		if (m_Loops[Inf.nLoop].m_bFront)return Inf.nLoop;
		return -1;
	}
	if (Inf.vertex >= 0) {
		for (int l = 0; l<m_nLoopsNum; l++) {
			if (m_Loops[l].m_bFront) {
				for (int p = 0; p<m_Loops[l].m_nSize; p++) {
					if (Inf.vertex == GetPointPos(l, p))return l;
				}
			}
		}
	}
	return -1;
}

/////////////////////////////////////////////////////////////////////////////
bool CFrontier::SetDirection(Math::C2DPoint **Points, FRONTINF **Inf, int nSize) {
	int n = 0;

	for (int i = 1; i<m_nPointsNum; i++) {
		if (m_Points[i].m_x < m_Points[n].m_x)n = i;
	}

	FRONTINF FI = { n,-1,-1 };
	int nFr = GetFrLoop(FI);

	if (nFr<0)return false;


	for (int l = 0; l<m_nLoopsNum; l++) {
		if (m_Loops[l].m_bFront) {
			bool bFirst = true;
			int minx = n;
			if (l != nFr) {
				minx = GetPointPos(l, 0);
				for (int p = 1; p<m_Loops[l].m_nSize; p++) {
					if (m_Loops[l].m_Points[p]->m_x < m_Points[minx].m_x)minx = GetPointPos(l, p);
				}
			}
			bool bPPT;
			Math::C2DPoint DP(m_Points[minx]), DP1;
			for (int i = 0; i<nSize; i++) {
				if (l == GetFrLoop(*Inf[i * 2]) && l == GetFrLoop(*Inf[i * 2 + 1])) {
					if (Inf[i * 2]->vertex == minx) {
						Math::C2DPoint DP2(*Points[i * 2 + 1]);
						if (bFirst || DP2.Right(DP, DP1) || (fabs(DP2.m_x - DP1.m_x) < EPS && DP2.m_y < DP1.m_y)) {
							DP1 = DP2;
							bFirst = false;
							bPPT = true;
						}
					}
					if (Inf[i * 2 + 1]->vertex == minx) {
						Math::C2DPoint DP2(*Points[i * 2]);
						if (bFirst || DP2.Right(DP, DP1) || (fabs(DP2.m_x - DP1.m_x) < EPS && DP2.m_y < DP1.m_y)) {
							DP1 = DP2;
							bFirst = false;
							bPPT = false;
						}
					}
				}
			}
			if (!bFirst && ((l == nFr && !bPPT) || (l != nFr && bPPT))) {
				for (int i = 0; i<nSize; i++) {
					if (l == GetFrLoop(*Inf[i * 2]) && l == GetFrLoop(*Inf[i * 2 + 1])) {
						Math::C2DPoint *pDP = Points[i * 2];
						Points[i * 2] = Points[i * 2 + 1];
						Points[i * 2 + 1] = pDP;
						FRONTINF *pInf = Inf[i * 2];
						Inf[i * 2] = Inf[i * 2 + 1];
						Inf[i * 2 + 1] = pInf;
					}
				}
			}
		}
	}
	return true;
}


/////////////////////////////////////////////////////////////////////////////
//CFinEl
CFinEl::CFinEl() {
	for (int i = 0; i<3;i++) {
		m_inf[i] = nullptr;
		m_Points[i] = nullptr;
		m_status[i] = 0;
		m_mark[i] = 0;
	}
}
CFinEl::~CFinEl()
{
}



/////////////////////////////////////////////////////////////////////////////
bool CFinEl::HaveRib(Math::C2DPoint *pDP1, Math::C2DPoint *pDP2) {
	for (int i = 0; i<3; i++) {
		if (pDP1 == m_Points[i] && pDP2 == m_Points[(i + 1) % 3]) {
			return true;
		}
	}
	return false;

}

/////////////////////////////////////////////////////////////////////////////
bool CFinEl::FrontCrossControl(CFrontier *pFront, int *b_pDPk) {

	Math::C2DPoint DPi, DPj, DPk, DPRes, DPTemp;
	Math::C2DPoint DP;
	int i, nearest;
	double distance, temp;
	bool bFirst = true;
	bool bRet = true;

	DPi = GetPi();
	DPj = GetPj();
	DPk = GetPk();


	nearest = -1;
	i = pFront->m_nPointsNum;
	while (i--) {
		DP = Math::C2DPoint(pFront->m_Points[i]);
		if (i != m_inf[0]->vertex && i != m_inf[1]->vertex && i != m_inf[2]->vertex && HaveInside(DP.m_x, DP.m_y)) {
			temp = DP.GetDistance(DPi, DPj, &DPRes);
			if (bFirst || temp < distance) {
				distance = temp;
				nearest = i;
				bFirst = false;
			}
		}
	}
	if (nearest != -1) {
		if (*b_pDPk == 0) {
			m_Points[2] = new Math::C2DPoint;
			m_inf[2] = new FRONTINF;
			*b_pDPk = 1;
		}
		*m_Points[2] = pFront->m_Points[nearest];
		m_inf[2]->vertex = nearest;
		m_inf[2]->nLoop = -1;
		m_inf[2]->rib = -1;
		DPk = GetPk();
		bRet = false;
	}


	Math::C2DPoint DPCentr(0.5*(DPi.m_x + DPj.m_x), 0.5*(DPi.m_y + DPj.m_y));
	FRONTINF OwnRib;
	if (!pFront->SameLine(*m_inf[0], *m_inf[1], &OwnRib)) {
		OwnRib.nLoop = -1;
		OwnRib.rib = -1;
		OwnRib.vertex = -1;
	}

	FRONTINF InfCross;
	if (pFront->Cross(DPCentr, DPk, OwnRib, *m_inf[2], &DPRes, &InfCross)) {
		if (*b_pDPk == 0) {
			m_Points[2] = new Math::C2DPoint;
			m_inf[2] = new FRONTINF;
			*b_pDPk = 1;
		}
		Math::C2DPoint DPx1, DPx2;
		bool bDPx1 = false, bDPx2 = false;
		int l = InfCross.nLoop, r = InfCross.rib;
		if (!pFront->SameLine(InfCross, *m_inf[0])) {
			if (DPx1.Cross(DPi, DPk, pFront->GetPoint(l, r), pFront->GetPoint(l, r + 1))) bDPx1 = true;
		}
		if (!pFront->SameLine(InfCross, *m_inf[1])) {
			if (DPx2.Cross(DPj, DPk, pFront->GetPoint(l, r), pFront->GetPoint(l, r + 1))) bDPx2 = true;
		}

		if (bDPx1 && bDPx2) {
			DPRes = (DPx1 + DPx2)*0.5;
		}
		else {
			if (bDPx1) DPRes = DPx1;
			if (bDPx2) DPRes = DPx2;
		}
		//1234
		*m_Points[2] = DPRes;
		*m_inf[2] = InfCross;
		DPk = GetPk();
		bRet = false;
	}
	return bRet;

}

/////////////////////////////////////////////////////////////////////////////
bool CFinEl::HaveInside(double x, double y) {
	Math::C2DPoint DP1, DP2;
	int nTop = 0, nBottom = 0;
	bool bLeft;

	DP1 = *m_Points[0];
	bLeft = DP1.m_x < x ? true : false;
	int i = 3;
	double y0;
	while (i--) {
		DP2 = DP1;
		DP1 = *m_Points[i];
		if ((DP1.m_x < x && !bLeft) || (DP1.m_x >= x && bLeft)) {
			bLeft = !bLeft;

			//Здесь всё нормально, знаменатель != 0
			y0 = (DP1.m_y - DP2.m_y) / (DP1.m_x - DP2.m_x)*(x - DP2.m_x) + DP2.m_y;
			if (y0>y) nTop++;
			else if (y0<y) nBottom++;
			else return true;
		}
	}

	if (nTop % 2 != 0) {
		return true;
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////
//Отрисовка конечных элементов (красный)
void CFinEl::DrawGl(GLParam &parameter)
{
	Math::C2DPoint DP0, DP1, DP2;
	int j = 0;

	glColor3d(1, 0, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); //GL_FILL для заполнения, GL_LINE для линий
	glBegin(GL_TRIANGLES);
	glVertex2dv(*m_Points[0]);
	glVertex2dv(*m_Points[1]);
	glVertex2dv(*m_Points[2]);
	glEnd();

	//Вершины границы
	glPointSize(6);
	glBegin(GL_POINTS);
	for (int i = 0; i<3; i++) {
		if (m_inf[i]->nLoop != -1 || m_inf[i]->vertex != -1)	glVertex2dv(*m_Points[i]);
	}
	glEnd();

	glPointSize(1);
	for (int i = 0; i<3; i++)
	{
		DP0 = *m_Points[i];
		DP1 = *m_Points[(i + 1) % 3];
		/*if(m_inf[i]->vertex != -1){
		//pDC->FillSolidRect(DP.m_x*Mult-100, DP.m_y*Mult-100, 200, 200, RGB(0, 150, 0));
		}else if(m_inf[i]->rib != -1){
		//pDC->FillSolidRect(DP.m_x*Mult-100, DP.m_y*Mult-100, 200, 200, RGB(0, 255, 0));
		}*/

		////////////////////////////////////////////
		/*
		glBegin(GL_LINES);
		glVertex2dv(DP);
		glVertex2dv(DP1);
		glEnd();
		*/
		////////////////////////////////////////////
		DP0.m_x = 0.5*(m_Points[i]->m_x + m_Points[(i + 1) % 3]->m_x);
		DP0.m_y = 0.5*(m_Points[i]->m_y + m_Points[(i + 1) % 3]->m_y);
		DP1.Rad(DP0, *m_Points[(i + 2) % 3], 0.2*DP0.Len(*m_Points[(i + 2) % 3]));

		if (m_status[i] && m_mark[i])
		{
			j++;
			////////////////////////////////////////////
			glColor3d(0, 1, 0);
			glBegin(GL_LINES);
			glVertex2dv(DP0);
			glVertex2dv(DP1);
			glEnd();
			////////////////////////////////////////////;
		}

		//if errror{{
		/*if(m_status[i] && ((i==2 && j==0) || (j==3)))
		{
		pDC->Rectangle( (int)(DP.m_x*Mult - 1000), (int)(DP.m_y*Mult - 1000),
		(int)(DP.m_x*Mult + 1000), (int)(DP.m_y*Mult + 1000));
		}*///if errror}}
	}
}//*/
 /////////////////////////////////////////////////////////////////////////////
