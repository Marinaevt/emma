#include "stdafx.h"
#include "Mesh.h"
#include "../Geometry/Mesher.h"
#include "../Geometry/GeometryUtils.h"


IOIMPL(CMesh, С2DMESH)

CMesh::CMesh() {

	m_Elements = nullptr;
	m_Points = nullptr;
	m_FrInf = nullptr;

	m_nElNum = 0;
	m_nElNumMax = 0;
	m_nPointsNum = 0;
	m_nPointsNumMax = 0;
	m_dRibLen = 1.0;
	m_dMeshRad = 1.0;
	m_IsSetPoint = false;
	m_nMeshType = ID_MESH_TYPE_MSW_COMMON;	//тоже, что и ID_MESH_TYPE_FRONT !!!

	RegisterMember(&m_MeshPoint);
	RegisterMember(&m_dRibLen);
	RegisterMember(&m_dMeshRad);
	RegisterMember(&m_IsSetPoint);
	RegisterMember(&m_nMeshType);
}

CMesh::~CMesh() {

	delete[] m_Elements;
	m_Elements = nullptr;

	delete[] m_Points;
	m_Points = nullptr;

	delete[] m_FrInf;
	m_FrInf = nullptr;
}

const Math::C2DPoint& CMesh::GetNode(size_t nNode) const {

	if (nNode >= m_nPointsNum) {
		CDlgShowError cError(ID_ERROR_MESH_NNODE_WRONG); //_T("nNode is wrong"));
		return Math::C2DPoint::Zero;
	}
	return m_Points[nNode];
}

const CFinEl& CMesh::GetElement(size_t nEl) const {

	if (nEl >= m_nElNum) {
		CDlgShowError cError(ID_ERROR_MESH_NEL_WRONG); //_T("nEl is wrong"));
		return m_Elements[0];	//TODO: вернуть пустой объект
	}

	return m_Elements[nEl];
}

/////////////////////////////////////////////////////////////////////////////
int	CMesh::AddPoint(const Math::C2DPoint &DP, const FRONTINF &FI) {
	m_nPointsNum++;
	//Проверяем на количество узлов в КЭ сетке
	if (m_nPointsNum > m_nPointsNumMax) {
		m_nPointsNumMax = m_nPointsNum + 10;	//TODO: Magic Number

		//Сохраняем старые массивы
		Math::C2DPoint *pOld = m_Points;
		FRONTINF *pFFOld = m_FrInf;
		
		//Новые массивы (бОльшего размера)
		m_Points = new Math::C2DPoint[m_nPointsNumMax];
		m_FrInf = new FRONTINF[m_nPointsNumMax];

		//Используем long long или pointer для x64 компиляторов (вдруг наступит светлое будущее)
		//Type           ILP64   LP64   LLP64
		//long             64     64      32
		//long long        64     64      64
		//pointer          64     64      64

		//Смещения массива точек и фронта
		long long offset = (long long)m_Points - (long long)pOld;
		long long FFoffset = (long long)m_FrInf - (long long)pFFOld;

		//Копируем старые массивы в новые
		if (pOld) {
			m_Points = (Math::C2DPoint *)memcpy((void *)m_Points, (void *)pOld, sizeof(Math::C2DPoint)*(m_nPointsNum - 1));
			m_FrInf = (FRONTINF *)memcpy((void *)m_FrInf, (void *)pFFOld, sizeof(FRONTINF)*(m_nPointsNum - 1));
			delete[] pOld;
			delete[] pFFOld;
		}

		//Смещаем узлы и фронты у Элементов
		for (size_t i = 0; i < m_nElNum; i++) {
			for (int j = 0; j < 3; j++) {
				m_Elements[i].m_Points[j] = (Math::C2DPoint *)((long long)m_Elements[i].m_Points[j] + offset);
				m_Elements[i].m_inf[j] = (FRONTINF *)((long long)m_Elements[i].m_inf[j] + FFoffset);
			}
		}
		//Смещаем узлы и фронты у RibLine
		m_RibLine.Offset(offset, FFoffset, (long long)m_Points - offset, (long long)(&m_Points[m_nPointsNum - 1]) - offset);
	}
	m_Points[m_nPointsNum - 1] = DP;
	m_FrInf[m_nPointsNum - 1] = FI;
	return m_nPointsNum - 1;

}

/////////////////////////////////////////////////////////////////////////////
int CMesh::AddElement(const CFinEl &FinEl) {
	m_nElNum++;
	if (m_nElNum > m_nElNumMax) {
		m_nElNumMax = m_nElNum + 10;

		//Сохраняем старые массивы
		CFinEl *pOld = m_Elements;
		//Новые массивы (бОльшего размера)
		m_Elements = new CFinEl[m_nElNumMax];
		
		//Копируем старые массивы в новые
		if (pOld) {
			m_Elements = (CFinEl *)memcpy(m_Elements, pOld, sizeof(CFinEl)*(m_nElNum - 1));
			delete[] pOld;
		}
	}
	m_Elements[m_nElNum - 1] = FinEl;
	return m_nElNum - 1;
}
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
//void CMesh::ConstructMap (){
//	int Size, i, j;
//
//
//	Size = m_FinElArray.GetSize();
//
//	for(i=0; i<Size; i++){
//		j = i;
//		while( !SetMarks(j) && j>=0){
//			m_FinElArray[j]->m_mark[0] = !m_FinElArray[j]->m_mark[0];
//			j = FindNeighbor(j, 0);
//		}
//	}
//}

//! Удаление сетки (без точки сгущения)
void CMesh::DeleteMesh()
{
	//if(m_Elements){
	delete[] m_Elements;
	m_Elements = nullptr;
	//}
	//if(m_Points){
	delete[] m_Points;
	m_Points = nullptr;
	//}
	//if(m_FrInf){
	delete[] m_FrInf;
	m_FrInf = nullptr;
	//}
	m_RibLine.Clear();
	m_nPointsNum = m_nPointsNumMax = m_nElNumMax = m_nElNum = 0;

	//m_dRibLen = 1.0;
}

//!Удаление точки сгущения
void CMesh::DeleteMeshPoint() {
	//m_MeshPoint = nullptr;
	m_dMeshRad = 0.0;
	m_IsSetPoint = false;
}

//! Возвращаем желаемую длину ребра элемента в данной точке с учетом сгущений
double CMesh::GetRibLen(const Math::C2DPoint &point) {

	double dRibLen = m_dRibLen;
	if (m_IsSetPoint) {

		if (m_MeshPoint().Len(point) < m_dMeshRad && m_MeshPoint().Len(point) > 0.5*m_dMeshRad) {
			dRibLen = 0.5*m_dRibLen;
		}
		if (m_MeshPoint().Len(point) < 0.5*m_dMeshRad) {
			dRibLen = 0.25*m_dRibLen;
		}

		//TODO: специфичный случай
		/*
		dRibLen = min(
			max(0.45*m_dRibLen, min(m_dRibLen, m_dRibLen*0.45 + 0.4*(point.y-15.) )),
			max(m_dRibLen*0.1, m_dRibLen*0.2*m_MeshPoint.Len(point))
			);
		if(point.x<6 && point.y>15) dRibLen *= 2.;
		//*/
	}
	return dRibLen;
}

bool CMesh::GenerateMesh(double dRibLen, C2DOutline* Outline) {
	
	bool ret = false;
	//DeleteMesh();

	switch (m_nMeshType) {
	case ID_MESH_TYPE_FRONT_COMMON:
		ret = GenerateMeshFront(dRibLen, 1.0, Outline);
		break;
	case ID_MESH_TYPE_MSW_COMMON :
		ret = GenerateMeshMSW(dRibLen, Outline);
		break;
	default:
		ret = GenerateMeshMSW(dRibLen, Outline);
		break;
	}

	return ret;
}

//! Создание сетки (тип, длина разбиения, угол, чертёж)
bool CMesh::GenerateMeshFront(double dRibLen, double Ang, C2DOutline* Outline)
{

	// MeshType - тип сетки
	// dRibLen - длина ребра сетки
	// Ang - угол конуса сгущения

	if (dRibLen <= 0.0) { return false; }
	m_dRibLen = dRibLen;

	CFinEl *pFinEl;

	if (!m_RibLine.Load(*Outline, this, dRibLen)) {
		return false;
	}

	int nStep = 0;
	//int nMaxStep = (Outline->GetSumLengthOfAllCurves() / m_dRibLen) *3;	//???
	//	CRibList* current = m_RibLine.m_pList;

	while (m_RibLine.m_pList)
	{
		if (m_RibLine.m_pList == nullptr)
			break;

		pFinEl = NewElement();
		m_RibLine.GetNext();
		m_RibLine.AddRib(pFinEl->m_Points[0], pFinEl->m_Points[2], pFinEl->m_inf[0], pFinEl->m_inf[2]);
		m_RibLine.AddRib(pFinEl->m_Points[2], pFinEl->m_Points[1], pFinEl->m_inf[2], pFinEl->m_inf[1]);

		//GetApp()->UpdateAllViews(); //TODO: Отладочная отрисовка сетки (поэлементно, долго)
		nStep++;
		
	}

	Optimize(5);

	return true;
}
bool CMesh::GenerateMeshMSW(double dRibLen, C2DOutline* Outline)
{

	// dRibLen - длина ребра сетки

	if (dRibLen <= 0.0) { return false; }
	m_dRibLen = dRibLen;

	// заполняем пул непересекающихся отрезков
	FemLibrary::NonIntersectedCurvesPool pool;
	for (size_t p = 0; p < Outline->GetCurveCount(); p++)
	{
		// кривые разбиваются на несколько отрезков таким же образом, как и в CRibLine::Load
		C2DCurve *pCur = Outline->GetCurve(p);
		Math::C2DPoint DP = pCur->GetStartNode()->GetPoint(), DPa;

		auto lastPoint = pCur->GetStartNode()->GetPoint();
		if(pCur->GetType()!= C2DLINESEG)
		{
			auto totalLen = pCur->CalcLength();
			for (DBL l1 = GetRibLen(DP); l1 < totalLen; l1 += GetRibLen(DP))
			{
				if (l1 > totalLen - GetRibLen(DP))
					l1 = 0.5*(l1 + totalLen - GetRibLen(DP));
				pCur->GetPoint(l1, DP, DPa);
				pool.AddLine(FemLibrary::LineSeg(lastPoint.x, lastPoint.y, DP.x, DP.y));

				lastPoint = DP;
			}
		}
		DP = pCur->GetEndNode()->GetPoint();
		pool.AddLine(FemLibrary::LineSeg(lastPoint.x, lastPoint.y, DP.x, DP.y));
	}

	// ищем контуры в массиве непересекеающихся отрезков
	auto contours = FemLibrary::GeometryUtils::FindContours(pool.Lines);
	
	// триангулируем каждый контур поотдельности
	FemLibrary::Mesh mesh;
	FemLibrary::MeshingParameters meshing_parameters;
	meshing_parameters.RibLength = m_dRibLen;
	for (size_t contourInd = 0; contourInd < contours.size(); contourInd++)
	{
		auto result = mesh.generateMesh(contours[contourInd], meshing_parameters);
		if (!result.success)
			return false;

		auto np = result.points.size();
		std::vector<int> pointIndexes(np, -1);
		for (int i = 0; i < np; ++i)
			pointIndexes[i] = AddPoint(CPoint2D(result.points[i].X, result.points[i].Y), FRONTINF{ -1, static_cast<int>(contourInd), -1 });

		auto nt = result.triangleIndexes.size() / 3;
		for (int i = 0; i < nt; ++i){
			CFinEl FinEl;
			int i0 = pointIndexes[result.triangleIndexes[i * 3    ]];
			int i1 = pointIndexes[result.triangleIndexes[i * 3 + 1]];
			int i2 = pointIndexes[result.triangleIndexes[i * 3 + 2]];
			FinEl.m_Points[0] = &m_Points[i0];
			FinEl.m_Points[1] = &m_Points[i1];
			FinEl.m_Points[2] = &m_Points[i2];
			FinEl.m_inf[0] = &m_FrInf[i0];
			FinEl.m_inf[1] = &m_FrInf[i1];
			FinEl.m_inf[2] = &m_FrInf[i2];
			AddElement(FinEl);
		}
	}

	return true;
}
/////////////////////////////////////////////////////////////////////////////
void CMesh::DrawGl(GLParam &parameter) const
{

	//точки
	if (m_IsSetPoint()) {
		glPointSize(5);
		glColor3d(1, 0, 1);
		glBegin(GL_POINTS);
		glVertex2dv(m_MeshPoint());
		glEnd();
	}

	/*
	for(int i = 0; i < this->m_nPointsNum; i++)
	{
		glColor3d((1-0.01*i),0,(0.01*i));
		glBegin(GL_POINTS);
			glVertex2dv(this->m_Points[i]);
		glEnd();
	} //*/

	//сетка
	glColor3d(1, 0, 0);
	for (size_t i = 0; i < m_nElNum; i++)
	{
		m_Elements[i].DrawGl(parameter);
	}
	glColor3d(0, 1, 0);//*/

//фронт
	//glBegin(GL_LINES);
	//CRibList* current = m_RibLine.m_pList;
	//do
	//{
	//	if(current)
	//	{		
	//		glVertex2dv(*current->pDP1);
	//		glVertex2dv(*current->pDP2);
	//		current = current->pNext;
	//	}

	//}
	//while(current != nullptr);
	//glEnd();

	glColor3d(0, 0, 1);
	CRibList *pList = m_RibLine.m_pList;

	glBegin(GL_LINES);
	while (pList) {
		glVertex2d(pList->pDP1->x, pList->pDP1->y);
		glVertex2d(pList->pDP2->x, pList->pDP2->y);
		pList = pList->pNext;
	}
	glEnd();

	glPointSize(5);
	glBegin(GL_POINTS);
	pList = m_RibLine.m_pList;
	while (pList) {
		glVertex2d(pList->pDP1->x, pList->pDP1->y);
		glVertex2d(pList->pDP2->x, pList->pDP2->y);
		pList = pList->pNext;
	}
	glEnd();
}
/////////////////////////////////////////////////////////////////////////////
CFinEl *CMesh::NewElement() {

	//const double dKff = 0.8; //не используется
	Math::C2DPoint DPi, DPj, DPk, DPRes, DPTemp;
	Math::C2DPoint *pDPi, *pDPj;
	int b_pDPk = 1;

	DPi = *m_RibLine.GetR1();
	DPj = *m_RibLine.GetR2();
	//double rib_len = DPi.Len(DPj);	//не используется

	// В DPk помещаем новую точку, образующую треугольник с DPi и DPj
	DPk.Third(DPi, DPj);

	FRONTINF FK = { -1,-1,-1 };
	AddPoint(DPk, FK);		//Помещаем новую точку во фронт

	int nDPk = m_nPointsNum - 1;
	CFinEl FinEl;
	FinEl.m_Points[0] = pDPi = m_RibLine.GetR1();
	FinEl.m_Points[1] = pDPj = m_RibLine.GetR2();
	FinEl.m_Points[2] = &m_Points[nDPk];
	FinEl.m_inf[0] = m_RibLine.GetInf1();
	FinEl.m_inf[1] = m_RibLine.GetInf2();
	FinEl.m_inf[2] = &m_FrInf[nDPk];

	DPTemp = (DPi + DPj)*0.5;	//Средняя точка на отрезке [DPi, DPj]
	DPRes = DPk;

	// Устанавливаем длину ребра, с учётом сгущения вокруг точек
	double dRibLen = GetRibLen(DPk);

	//! Устанавливаем новую точку КЭ сетки с учётом сгущения
	double kff = sin(M_PI / 3); 
	DPk.Rad(DPTemp, DPk, dRibLen*kff);

	*FinEl.m_Points[2] = DPk;

	m_RibLine.IfCross(&FinEl, &b_pDPk);
	m_RibLine.IfVertArea(&FinEl, &b_pDPk, dRibLen*0.6);
	m_RibLine.IfCross(&FinEl, &b_pDPk);
	//m_RibLine.IfTooClose(&FinEl, &b_pDPk);

	// используем (long long) или (pointer) для x64
	//Type           ILP64   LP64   LLP64
	//long             64     64      32
	//long long        64     64      64
	//pointer          64     64      64

	if (b_pDPk == 0) {
		if ((long long)FinEl.m_Points[2] >= (long long)(m_Points) && (long long)FinEl.m_Points[2] < (long long)(&m_Points[m_nPointsNum])) {
			m_nPointsNum--;
		}
		else {
			m_Points[nDPk] = *FinEl.m_Points[2];
			m_FrInf[nDPk] = *FinEl.m_inf[2];
			m_RibLine.ChangePoint(FinEl.m_Points[2], &m_Points[nDPk], &m_FrInf[nDPk]);
			FinEl.m_Points[2] = &m_Points[nDPk];
			FinEl.m_inf[2] = &m_FrInf[nDPk];
		}
	}
	int nElemNum = AddElement(FinEl);
	return &m_Elements[nElemNum];
}


/////////////////////////////////////////////////////////////////////////////
void CMesh::ConstructMap()
{
	for (size_t i = 0; i < m_nElNum; i++) //проходим по всем элементам
	{
		size_t j = i;
		while (!SetMarks(j))
		{
			m_Elements[j].m_mark[0] = !m_Elements[j].m_mark[0];
			if (FindNeighbor(j, 0) < 0) {
				break;
			}
			j = FindNeighbor(j, 0);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
int CMesh::FindNeighbor(size_t num, int rib, int *NeigRib)
{
	Math::C2DPoint *pDP1, *pDP2;

	pDP1 = m_Elements[num].m_Points[rib];
	pDP2 = m_Elements[num].m_Points[(rib + 1) % 3];


	for (size_t n = 0; n < m_nElNum; n++)
	{
		for (int i = 0; i < 3; i++)
		{
			if (m_Elements[n].m_Points[i] == pDP2 && m_Elements[n].m_Points[(i + 1) % 3] == pDP1)
			{
				if (NeigRib != nullptr) *NeigRib = i;
				return n;
			}
		}
	}
	return -1;
}
/////////////////////////////////////////////////////////////////////////////
bool CMesh::SetMarks(int index)
{
	int Neighbor, N_Rib, Rib;
	int b_Prior = 1, 
		First = 1; //b_Old; //не используется

	for (Rib = 0; Rib < 3; Rib++)
	{
		//b_Old = b_Prior; //не используется
		Neighbor = FindNeighbor(index, Rib, &N_Rib);
		if (Neighbor == -1)
		{
			First = 0;
			b_Prior = m_Elements[index].m_mark[Rib] = !b_Prior;
		}
		else
		{
			if (m_Elements[Neighbor].m_status[N_Rib])
			{
				b_Prior = m_Elements[index].m_mark[Rib] = !(m_Elements[Neighbor].m_mark[N_Rib]);
			}
			else
			{
				First = 0;
				b_Prior = m_Elements[index].m_mark[Rib] = !b_Prior;
			}
		}
		m_Elements[index].m_status[Rib] = 1;
	}
	if ((m_Elements[index].m_mark[0] && m_Elements[index].m_mark[1] && m_Elements[index].m_mark[2]) ||
		(!m_Elements[index].m_mark[0] && !m_Elements[index].m_mark[1] && !m_Elements[index].m_mark[2]))
	{
		if (!First)
		{
			m_Elements[index].m_mark[0] = !m_Elements[index].m_mark[0];
			return true;
		}
		return false;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////
void CMesh::Renumber(Math::C2DPoint ***Stars, int *StarSize)
{
	int *NewNums = new int[m_nPointsNum];
	int *NewNums_1 = new int[m_nPointsNum];
	size_t k = 0, k0 = 0;
	for (size_t s = 0; s < m_nPointsNum; s++)
	{
		for (int j = 0; j < StarSize[s]; j++)
		{
			Math::C2DPoint DPs(m_Points[s]);		//Узлы 
			Math::C2DPoint DPj(*Stars[s][j]);		//Звёзды
			Math::C2DPoint DP;
			DP.Rad(DPs, DPj, 0.1*DPs.Len(DPj));		//Формируем новую точку
		}
	}

	NewNums[0] = 0;
	NewNums_1[0] = 0;	//TODO: Нигде не используем
	while (k0 < m_nPointsNum)
	{
		int ii = 0;
		int NewMinS = 2;
		while (ii < StarSize[NewNums[k0]])
		{
			int MinS = NewMinS;
			for (int i = 0; i < StarSize[NewNums[k0]]; i++)
			{
				int N = (int)((long long)Stars[NewNums[k0]][i] - (long long)m_Points) / sizeof(Math::C2DPoint);
				if (StarSize[N] == MinS)
				{
					bool bNew = true;
					for (size_t l = 0; l <= k; l++)
					{
						if (N == NewNums[l])
						{
							bNew = false;
							break;
						}
					}
					if (bNew)
					{
						NewNums[++k] = N;
						NewNums_1[N] = k;

						Math::C2DPoint DPs(m_Points[NewNums[k0]]);
						Math::C2DPoint DPj(*Stars[NewNums[k0]][i]);
						Math::C2DPoint DP;
						DP.Rad(DPs, DPj, DPs.Len(DPj));
					}
					ii++;
				}
				else if (StarSize[N] > MinS && (NewMinS == MinS || StarSize[N] < NewMinS))
				{
					NewMinS = StarSize[N];
				}
			}
		}
		k0++;
		if (k0 > k && k0 < m_nPointsNum)
		{
			delete[] NewNums;
			delete[] NewNums_1;
			return;
		}
	}

	Math::C2DPoint *NewPoints = new Math::C2DPoint[m_nPointsNum];
	FRONTINF  *NewInf = new FRONTINF[m_nPointsNum];

	for (size_t i = 0; i < m_nPointsNum; i++)
	{
		NewPoints[i] = m_Points[NewNums[i]];
		NewInf[i] = m_FrInf[NewNums[i]];
	}

	for (size_t e = 0; e < m_nElNum; e++)
	{
		for (int j = 0; j < 3; j++)
		{
			int N = (int)((long long)m_Elements[e].m_Points[j] - (long long)m_Points) / sizeof(Math::C2DPoint);
			int Ne = (int)((long long)m_Elements[e].m_inf[j] - (long long)m_FrInf) / sizeof(FRONTINF);

			if (Ne != N) {
				CDlgShowError cError(ID_ERROR_MESH_NE_NOT_N); //_T("Ne != N"));
				return;
			}

			m_Elements[e].m_Points[j] = &NewPoints[NewNums_1[N]];
			m_Elements[e].m_inf[j] = &NewInf[NewNums_1[N]];
		}
	}
	Math::C2DPoint *pOldPoints = m_Points;
	FRONTINF *pOldInf = m_FrInf;
	m_Points = NewPoints;
	m_FrInf = NewInf;

	delete[] pOldPoints;
	delete[] pOldInf;
	delete[] NewNums;
	delete[] NewNums_1;
}

/////////////////////////////////////////////////////////////////////////////
void CMesh::Optimize(int n) {

	//FRONTINF FI = {-1, -1, -1}; //не используется
	Math::C2DPoint ***Stars = new Math::C2DPoint**[m_nPointsNum];
	int *StarSize = new int[m_nPointsNum];
	bool *Front = new bool[m_nPointsNum];

	for (int k = 0; k < n; k++) {
		for (size_t i = 0; i < m_nPointsNum; i++) {
			if (k == 0) {
				StarSize[i] = GetStar(i, &Stars[i]);
				Front[i] = !(m_FrInf[i].nLoop == -1 && m_FrInf[i].vertex == -1);
			}

			if (!Front[i]) {
				bool bFlag = true;
				for (int j = 0; j < StarSize[i]; j++) {
					Math::C2DPoint DP0(*Stars[i][j]), DP1(*Stars[i][(j + 1) % StarSize[i]]), DP2(*Stars[i][(j + 2) % StarSize[i]]);
					if (!DP1.Right(DP0, DP2)) {
						bFlag = false;
						break;
					}
				}

				if (bFlag) {
					Math::C2DPoint DP(0, 0);
					for (int j = 0; j < StarSize[i]; j++) {
						DP.m_x += Stars[i][j]->m_x;
						DP.m_y += Stars[i][j]->m_y;
					}
					DP.m_x /= static_cast<double>(StarSize[i]);
					DP.m_y /= static_cast<double>(StarSize[i]);
					m_Points[i] = DP;
				}
			}
		}
	}
	Renumber(Stars, StarSize);

	for (size_t i = 0; i < m_nPointsNum; i++) {
		delete[] Stars[i];
	}
	delete[] Stars;
	delete[] Front;
	delete[] StarSize;
}
/////////////////////////////////////////////////////////////////////////////

void CMesh::SetOutline()
{
}

//!
int CMesh::GetStar(int nPoint, Math::C2DPoint ***Star) {
	int nENum = 0;
	for (size_t i = 0; i < m_nElNum; i++) {
		for (int j = 0; j < 3; j++) {
			if (m_Elements[i].m_Points[j] == &m_Points[nPoint])nENum++;
		}
	}
	Math::C2DPoint **StarTmp = new Math::C2DPoint *[nENum * 2];
	int p = 0;
	int nPNum = 0;
	for (size_t i = 0; i < m_nElNum; i++) {
		for (int j = 0; j < 3; j++) {
			if (m_Elements[i].m_Points[j] == &m_Points[nPoint]) {
				StarTmp[p++] = m_Elements[i].m_Points[(j + 1) % 3];
				StarTmp[p++] = m_Elements[i].m_Points[(j + 2) % 3];
				nPNum += 2;
				for (int pp = 0; pp < p - 2; pp++) {
					if (StarTmp[pp] == StarTmp[p - 1] || StarTmp[pp] == StarTmp[p - 2])nPNum--;
				}
			}
		}
	}
	(*Star) = new Math::C2DPoint *[nPNum];
	(*Star)[0] = StarTmp[0];
	(*Star)[1] = StarTmp[1];
	if (nPNum == nENum) {
		for (int i = 2; i < nENum; i++) {
			//bool bFind = false; //не используется
			for (int j = 1; j < nENum; j++) {
				if (StarTmp[2 * j] == (*Star)[i - 1]) {
					(*Star)[i] = StarTmp[2 * j + 1];
					//bFind = true; //не используется
					break;
				}
			}
		}
	}
	else {
		p = 2;
		for (int i = 2; i < nENum * 2; i++) {
			bool bNew = true;
			for (int j = 0; j < p; j++) {
				if ((*Star)[j] == StarTmp[i]) {
					bNew = false;
					break;
				}
			}
			if (bNew) (*Star)[p++] = StarTmp[i];
		}
	}
	delete[] StarTmp;
	return nPNum;
}


void CMesh::SetMeshType(UNINT id) {
	m_nMeshType() = id;
}


void CMesh::WriteToLog() const
{
	/*DLOG(CString(_T("CMesh")), log_info);
	DLOG(CString(_T("Point count: "))+AllToString(m_nPointsNum)+CString(_T(", Element count: "))+AllToString(m_nElNum), log_info);
	CString tmp;

	for (int i=0; i<m_nElNum; i++){
		tmp += CString("Element[")+AllToString(i)+CString("]");
		for (int j=0; j<3; j++){
			tmp += CString("(")+AllToString(GetElement(i).m_Points[j]->x)+CString(",")+AllToString(GetElement(i).m_Points[j]->y)+CString(") ");
		}
		DLOG(tmp, log_info);
		tmp="";
	}*/
}