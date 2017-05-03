#include "StdAfx.h"
#include "2DMesh.h"

// --- отладка битовых функций ---
void WriteBitFuncsToLog()
{
	DLOG(CString(_T("[dec]100725389 = [bin]0000 0110 0000 0000 1111 0010 1000 1101")), log_info);

	DLOG(CString(_T("--- BitCount(src, stop): количество установленных бит от младшего до бита с номером stop ---")), log_info);
	DLOG(CString(_T("BitCount(100725389) = ")) + AllToString(BitCount(100725389)), log_info);
	DLOG(CString(_T("BitCount(100725389, 8) = ")) + AllToString(BitCount(100725389, 8)), log_info);
	DLOG(CString(_T("BitCount(0) = ")) + AllToString(BitCount(0)), log_info);
	DLOG(CString(_T("BitCount(8, 3) = ")) + AllToString(BitCount(8, 3)), log_info);
	DLOG(CString(_T("BitCount(8, 4) = ")) + AllToString(BitCount(8, 3)), log_info);

	DLOG(CString(_T("--- Index2BitIndex(src, index): номер index-го установленного бита (-1, если его нет) ---")), log_info);
	DLOG(CString(_T("Index2BitIndex(100725389, 0) = ")) + AllToString(Index2BitIndex(100725389, 0)), log_info);
	DLOG(CString(_T("Index2BitIndex(100725389, 1) = ")) + AllToString(Index2BitIndex(100725389, 1)), log_info);
	DLOG(CString(_T("Index2BitIndex(100725389, 2) = ")) + AllToString(Index2BitIndex(100725389, 2)), log_info);
	DLOG(CString(_T("Index2BitIndex(100725389, 3) = ")) + AllToString(Index2BitIndex(100725389, 3)), log_info);
	DLOG(CString(_T("Index2BitIndex(100725389, 11) = ")) + AllToString(Index2BitIndex(100725389, 11)), log_info);

	DLOG(CString(_T("--- HasBit(int src, int bits): есть ли хотя бы один бит из bits в src")), log_info);
	DLOG(CString(_T("HasBit(100725389, 1) = ")) + AllToString(HasBit(100725389, 1)), log_info);
	DLOG(CString(_T("HasBit(100725389, 2) = ")) + AllToString(HasBit(100725389, 2)), log_info);
	DLOG(CString(_T("HasBit(100725389, 4) = ")) + AllToString(HasBit(100725389, 4)), log_info);
	DLOG(CString(_T("HasBit(100725389, 8) = ")) + AllToString(HasBit(100725389, 8)), log_info);


	DLOG(CString(_T("--- GetBitIndex(int src, int bits): порядковый номер первого установленного бита из bits в src ---")), log_info);
	DLOG(CString(_T("GetBitIndex(100725389, 1) = ")) + AllToString(GetBitIndex(100725389, 1)), log_info);
	DLOG(CString(_T("GetBitIndex(100725389, 2) = ")) + AllToString(GetBitIndex(100725389, 2)), log_info);
	DLOG(CString(_T("GetBitIndex(100725389, 4) = ")) + AllToString(GetBitIndex(100725389, 4)), log_info);
	DLOG(CString(_T("GetBitIndex(100725389, 8) = ")) + AllToString(GetBitIndex(100725389, 8)), log_info);

	DLOG(CString(_T("---")), log_info);
}

IOIMPL(C2DMesh3, CIO2DMESH3)

// возвращает количество установленных бит в src начиная с нулевого бита и заканчивая битом под номером stop
int BitCount(int src, UNINT stop)
{
	int ret = 0;

	for (UNINT i = 0; i < stop; i++)
	{
		ret += src & 0x1;
		src >>= 1;

		if (!src) // ранний выход, если дальше единиц точно нет
			return ret;
	}

	return ret;
}

// возвращает номер index-го установленного бита (-1, если его нет)
int Index2BitIndex(int src, int index)
{
	int ret = -1;

	for (int i = 0; i < 32; i++)
	{
		if ((src & 0x1))
		{
			if (!index)
			{
				ret = i;
				break;
			}
			src >>= 1;
			index -= 1;
		}
	}

	return ret;
}

// возвращает число с установленным битом index
int BitIndex2Number(int index)
{
	return (0x1 << index);
}

// возвращает true, если (src & bits)
bool HasBit(int src, int bits)
{
	if ((src & bits) != 0) {
		return true;
	}
	return false;
}

// возвращает порядковый номер первого установленного бита из bits в src
int GetBitIndex(int src, int bits)
{
	int tmp = src & bits;
	if (!tmp)
		return 0;

	int ret = 0;

	while (!(tmp & 0x1))
	{
		tmp >>= 1;
		ret++;
	}

	return ret;
}

///////////	C2DElement3

int C2DElement3::HaveEdge(ptrdiff_t node1, ptrdiff_t node2) const
{
	if (node1 == n0)
	{
		if (node2 == n1) return 0;
		if (node2 == n2) return 2;
	}
	else if (node1 == n1)
	{
		if (node2 == n0) return 0;
		if (node2 == n2) return 1;
	}
	else if (node1 == n2)
	{
		if (node2 == n0) return 2;
		if (node2 == n1) return 1;
	}
	return -1;
}

bool C2DElement3::IsHaveNode(ptrdiff_t nNode) {
	if (n0 == nNode) return true;
	if (n1 == nNode) return true;
	if (n2 == nNode) return true;
	return false;
}

///////////	C2DMeshInterface

void C2DMeshInterface::RegisterMembers() {
	RegisterMember(&m_nodes);
	RegisterMember(&m_nfields);
	RegisterMember(&m_efields);
	RegisterMember(&m_nodefields);
	RegisterMember(&m_bordernodes);
	RegisterMember(&m_borderelems);
	//RegisterMember(&m_bordernodes_marks);
}

bool C2DMeshInterface::InitNodes(size_t nodes, int nfields)
{
	if (!nfields)
		return false;

	bool ret;

	m_nodes().resize(nodes);
	ret = m_nodefields().resize(nodes * BitCount(nfields));

	if (!ret)
		return ret;

	m_nfields() = nfields;

	return true;
}

bool C2DMeshInterface::InitNodeFields(int nfields)
{
	if (!nfields)
		return false;

	bool ret = m_nodefields().resize(m_nodes().size() * BitCount(nfields));

	if (!ret)
		return false;


	m_nfields() = nfields;

	return true;
}

Math::C2DPoint C2DMeshInterface::GetBorderNode(size_t nBorderNode) const {

	//if (nBorderNode < 0 || nBorderNode >= m_bordernodes.GetSize()) {
	if (nBorderNode >= m_bordernodes.GetSize()) {
		CDlgShowError cError(ID_ERROR_2DMESH_NBORDERNODE); //_T("nBorderNode is wrong"));
		return Math::C2DPoint::Zero;
	}
	return m_nodes[m_bordernodes[nBorderNode]];
}

DBL& C2DMeshInterface::GetNField(size_t node, int field)
{
	//size_t flds = m_nfields();
	int index = GetBitIndex(m_nfields(), field);
	int count = BitCount(m_nfields());

	return m_nodefields().cell2d(node, index, count);
}

DBL C2DMeshInterface::GetNField(size_t node, int field) const
{
	//size_t flds = m_nfields();
	int index = GetBitIndex(m_nfields(), field);
	int count = BitCount(m_nfields());

	return m_nodefields().cell2d(node, index, count);
}

void C2DMeshInterface::AddBorderNode(int node)
{
	for (size_t i = 0; i < m_bordernodes().size(); i++)
		if (node == m_bordernodes()[i])
			return;

	m_bordernodes().push_back(node);
}

void C2DMeshInterface::AddBorderElem(int elem)
{
	for (size_t i = 0; i < m_borderelems().size(); i++)
		if (elem == m_borderelems()[i])
			return;

	m_borderelems().push_back(elem);
}

Math::C2DRect C2DMeshInterface::GetBoundingBox() const
{
	Math::C2DRect rect;

	const std::vector<Math::C2DPoint> &nodes = m_nodes();

	if (m_nodes().size())
	{
		rect.SetPoint(nodes[0]);

		for (size_t i = 1; i < nodes.size(); i++)
			rect.AddPoint(nodes[i]);
	}

	return rect;
}

bool C2DMeshInterface::GetBoundingBox(Math::C2DRect &rect) {
	if (m_nodes().size())
	{
		rect.SetPoint(m_nodes()[0]);

		for (size_t i = 1; i < m_nodes().size(); i++) {
			rect.AddPoint(m_nodes()[i]);
		}
		return true;
	}
	return false;
}

//Функция возвращает вектор, направленный по нормали к поверхности в граничной точке nBoundaryNode
Math::C2DPoint C2DMeshInterface::GetNormalS(ptrdiff_t nBoundaryNode, int nSMode) const {
	ptrdiff_t nBorderLen = GetBorderSize() - 1;

	if (nBorderLen < 1) {
		CDlgShowError cError(ID_ERROR_2DMESH_NBORDERLEN_0); //_T("nBorderLen <= 0"));
		return Math::C2DPoint::Zero;	//нулевой вектор
	}

	Math::C2DPoint P0 = GetBorderNode(nBoundaryNode);
	Math::C2DPoint P1 = GetBorderNode(nBoundaryNode < nBorderLen - 1 ? nBoundaryNode + 1 : 0);

	Math::C2DPoint n = Math::C2DPoint(P1.y - P0.y, P0.x - P1.x)*0.5; // получаем P1-P0 и поворачиваем на 90 градусов cw

	if (nSMode == 1) {
		n *= 0.25*(3 * P1.x + P0.x);
	}
	else if (nSMode == -1) {
		n *= 0.25*(P1.x + 3 * P0.x);
	}

	return n;
}

//Функция, возвращает площадь кольца вокруг граничного узла nBoundaryNode
DBL C2DMeshInterface::GetCircleSquare(ptrdiff_t nBoundaryNode) const {
	ptrdiff_t nBorderLen = GetBorderSize() - 1;
	
	if (nBorderLen < 1) {
		CDlgShowError cError(ID_ERROR_2DMESH_NBORDERLEN_0); //_T("nBorderLen <= 0"));
		return 0.0;	//нулевая площадь
	}

	Math::C2DPoint P0 = GetBorderNode(nBoundaryNode);
	Math::C2DPoint P1 = GetBorderNode(nBoundaryNode < nBorderLen - 1 ? nBoundaryNode + 1 : 0);	//вправо на 1 узел
	Math::C2DPoint P2 = GetBorderNode(nBoundaryNode > 0 ? nBoundaryNode - 1 : nBorderLen - 1);	//влево

	//Получаем точки на 1/2 расстояния
	Math::C2DPoint n1 = Math::C2DPoint(P0 + P1)*0.5; // получаем (P1+P0)/2
	Math::C2DPoint n2 = Math::C2DPoint(P0 + P2)*0.5; // получаем (P2+P0)/2

	//TODO: случай прямого угла
	//Осесимметричный случай
	DBL dS_upper = M_PI*n1.x*n1.x;
	DBL dS_lower = M_PI*n2.x*n2.x;

	return fabs(dS_upper - dS_lower);
}


void C2DMeshInterface::ReorderNodes(const std::vector<size_t>& neworder)
{

	std::vector<size_t> temp(neworder);

	for (size_t i = 0; i < neworder.size(); i++) {
		temp[neworder[i]] = i;
	}

	/*
	for (size_t i = 0; i < nodeN_; i++) {
	// вариант 1
	newNodes[temp[i]] = nodes_[i];

	// вариант 2
	//newNodes[i] = nodes_[newOrder[i]];
	}
	*/

	// узлы
	CIO2DPoints tmp1;
	tmp1().resize(m_nodes().size());
	for (size_t i = 0; i < m_nodes().size(); i++)
	{
		tmp1()[temp[i]] = m_nodes()[i];
	}

	for (size_t i = 0; i < m_nodes().size(); i++)
	{
		m_nodes()[i] = tmp1()[i];
	}

	/*
	// граничные узлы
	IO::CIntArray tmp2;
	tmp2().SetSize()
	for (size_t i = 0; i < m_bordernodes().GetCount(); i++)
	{
	tmp2()[i] = neworder[m_bordernodes()[i]];
	}
	m_bordernodes() = tmp2;

	// поля узлов
	IO::CCustomDoubleArray tmp3;
	*/
}

/////////	C2DMesh3

C2DMesh3::C2DMesh3(const CMesh &mesh) {
	RegisterMembers();
	//Копирование узлов
	m_nodes().resize(mesh.GetNodesNum());
	for (size_t i = 0; i<m_nodes().size(); i++) {
		m_nodes()[i] = mesh.GetNode(i);
	}
	//Копирование эл-ов
	m_elements().resize(mesh.GetElNum());
	for (size_t i = 0; i<m_elements().size(); i++) {
		for (int j = 0; j<3; j++) {//TODO: Переписать после исправления CFinEl
			m_elements()[i].n[j] = mesh.GetNodeInEl(i, j);
		}
	}
}

void C2DMesh3::RegisterMembers() {
	C2DMeshInterface::RegisterMembers();
	RegisterMember(&m_elements);
	RegisterMember(&m_elemfields);
}

/// выделяет память под элементы и поля элементов
bool C2DMesh3::InitElements(size_t elements, int efields)
{
	if (!efields)
		return false;

	bool ret;

	m_elements().resize(elements);
	ret = m_elemfields().resize(elements * BitCount(efields));


	if (!ret)
		return ret;

	m_efields() = efields;

	return true;
}

bool C2DMesh3::InitElemFields(int efields)
{
	if (!efields)
		return false;

	bool ret = m_elemfields().resize(m_elements().size() * BitCount(efields));
	if (!ret)
		return false;


	m_efields() = efields;

	return true;
}

DBL& C2DMesh3::GetEField(ptrdiff_t elem, int field)
{
	return m_elemfields().cell2d(elem, GetBitIndex(m_efields(), field), BitCount(m_efields()));
}

//! Заполнение матрицы жёсткости (rp -> std::vector)
void C2DMesh3::FillSTM(BANDMATRIX& m, std::vector<DBL>& rp, const CMaterial& mat, DBL dt, size_t startE, size_t endE)
{
	//LOGGER.Init(CString("..\\..\\Logs\\C2DMesh.cpp_FillSTM.txt"));

	C2DPoints& nodes = m_nodes();
	std::vector<C2DElement3>& elems = m_elements();

	//if (!Math::CheckBoundsInc(0, (int)elems.size(), startE, endE))
	if (startE > elems.size() || endE > elems.size())
		return;

	if (!HasBit(m_efields, eFields::mju))
		return;

	for (size_t i = startE; i < endE; i++) {
		C2DElement3& e = elems[i];

		Math::C2DPoint& n0 = nodes[e._0]; // 0 - i
		Math::C2DPoint& n1 = nodes[e._1]; // 1 - j
		Math::C2DPoint& n2 = nodes[e._2]; // 2 - m k

		DBL B[3], C[3];
		B[0] = n1.m_y - n2.m_y; //Y[j]-Y[k];
		B[1] = n2.m_y - n0.m_y; //Y[k]-Y[i];
		B[2] = n0.m_y - n1.m_y; //Y[i]-Y[j];
		C[0] = n2.m_x - n1.m_x; //X[k]-X[j];
		C[1] = n0.m_x - n2.m_x; //X[i]-X[k];
		C[2] = n1.m_x - n0.m_x; //X[j]-X[i];
		DBL sm = abs(n0.m_x * B[0] + n2.m_x * B[2] + n1.m_x * B[1]) * 0.5; //(X[i]*(B[0]) + X[k]*(B[2]) + X[j]*(B[1]))/2;

		DBL dMju = GetEField(i, eFields::mju); // вычисляется в зависимости от интенсивности напряжений и деформаций
		DBL dKk = mat.CalcK(dt, dMju);
		DBL dK = mat.GetK(dt, dMju);

		//DBL K1 = 4.0 / 3.0 * mju + K * dt;
		//DBL K2 = K * dt - 2.0 / 3.0 * mju;

		//Для осесимметричных задач{{
		DBL x = (n0.m_x + n1.m_x + n2.m_x) / 3.0;
		//Для осесимметричных задач}}//*/

		for (int s = 0; s < 3; s++) {
			ptrdiff_t ii = e.n[s];
			DBL b1 = B[s];
			DBL c1 = C[s];

			for (int t = 0; t < 3; t++) {
				ptrdiff_t jj = e.n[t];
				DBL b2 = B[t];
				DBL c2 = C[t];

				if (ii >= jj) {

					DBL g11 = ((2 * dMju + dKk)*b1*b2 + dMju*c1*c2) / 4 / sm;
					DBL g12 = (dKk*b2*c1 + dMju*c2*b1) / 4 / sm;
					DBL g21 = (dKk*b1*c2 + dMju*c1*b2) / 4 / sm;
					DBL g22 = ((2 * dMju + dKk)*c1*c2 + dMju*b1*b2) / 4 / sm;


					//Для осесимметричных задач{{
					g11 = x*g11 + (dKk*(b1 + b2) / 6 + (2 * dMju + dKk)*sm / 9 / x);	//исправил b1*b2 на b1+b2 (в "Сверхпластичности" ошибка)
					g12 = x*g12 + dKk*c1 / 6;
					g21 = x*g21 + dKk*c2 / 6;
					g22 = x*g22;
					//Для осесимметричных задач}}//*/

					m.cell(jj * 2, ii * 2) += g11;
					m.cell(jj * 2, ii * 2 + 1) += g12;
					if (ii>jj)
						m.cell(jj * 2 + 1, ii * 2) += g21;
					m.cell(jj * 2 + 1, ii * 2 + 1) += g22;
					if (isnan(m.cell(jj * 2, ii * 2))) {
						g11 = g11;
					}
				}


				/*//Для плоских задач{{
				rp[2 * jj]     -= GetEField(i, eFields::smo_d) * dK * b2 / 6;
				rp[2 * jj + 1] -= GetEField(i, eFields::smo_d) * dK * c2 / 6;
				//Для плоских задач}}
				/*/
				//Для осесимметричных задач{{
				rp[2 * jj] -= GetEField(i, eFields::smo_d) * dK * (b2 * x / 2 + sm / 3.0) / 3;
				rp[2 * jj + 1] -= x*GetEField(i, eFields::smo_d) * dK * c2 / 6;
				//Для осесимметричных задач}}//*/
				//rp[2 * jj] -= 2 * M_PI * x * sm * (GetEField(i, eFields::smo_d)*dK * (b2 / (2 * sm) + 1 / (3 * x)));
				//rp[2 * jj + 1] -= 2 * M_PI * x * sm * (GetEField(i, eFields::smo_d)*dK * c2 / (2 * sm));
			}
		}
	}
}


void C2DMesh3::SortNodesInElements(bool ccw)
{

	C2DPoints& nodes = m_nodes();
	std::vector<C2DElement3>& elems = m_elements();

	for (size_t i = 0; i < elems.size(); i++)
	{
		Math::C2DPoint &n0 = nodes[elems[i].n0],
			&n1 = nodes[elems[i].n1],
			&n2 = nodes[elems[i].n2];

		if (n0.Right(n1, n2) && ccw)
			Math::swap(elems[i].n0, elems[i].n1);
	}

}

namespace std {
	template <> struct hash<CEdge> {
		std::size_t operator()(const CEdge& p) const { return CEdge::GetHash(p); }
	};
}

//! Заполнение граничных узлов
void C2DMesh3::FillBorderNodes()
{
	size_t elem_size = m_elements().size();
	if (elem_size == 0)
		return;

	std::unordered_set<CEdge> edges_hash;

	//формируем список уникальных рёбер (граничные == уникальные), которые принадлежат только 1 эл-ту
	for (size_t i = 0; i < elem_size; i++)
	{
		C2DElement3& el = m_elements()[i];


		//в элементе 3 ребра
		for (int j = 0; j < 3; j++) {
			//получаем индексы вида (0,1), (1,2), (2,0)
			CEdge edge_from_element(el.n[j],  el.n[FemLibrary::Utils::NextInd(j, 3)]);
			CEdge edge_from_element_other_way(el.n[FemLibrary::Utils::NextInd(j, 3)], el.n[j]);	//то же ребро, но "в другую сторону"

			//ищем оба таких ребра
			std::unordered_set<CEdge>::const_iterator got = edges_hash.find(edge_from_element);
			//std::unordered_set<CEdge>::const_iterator got_other_way = edges_hash.find(edge_from_element_other_way);
			
			if (got == edges_hash.end()){ //&& got_other_way == edges_hash.end()) {	//ничего не нашли
				edges_hash.insert(edge_from_element);
			}
			else {

				edges_hash.erase(got);	//удаляем, т.к. повтор, значит ребро не внешнее \ не граничное
				/*
				if(got == edges_hash.end()){
					edges_hash.erase(got_other_way);	//удаляем, т.к. повтор, значит ребро не внешнее \ не граничное
				}
				if (got_other_way == edges_hash.end()) {
					edges_hash.erase(got);	//удаляем, т.к. повтор, значит ребро не внешнее \ не граничное
				}*/
			}

		}
	}

	//очищаем всё, что было раньше
	m_bordernodes().clear();

	std::unordered_set<CEdge>::const_iterator it = edges_hash.begin();
	// добавляем первые узлы 
	m_bordernodes().push_back(it->n0);
	m_bordernodes().push_back(it->n1);
	//текущий "последний узел" запоминаем, чтобы выстроить последовательно узлы
	ptrdiff_t cur_node = it->n1;
	edges_hash.erase(it);
	it = edges_hash.begin();

	while (edges_hash.size() > 1)	// 1 чтобы отсечь второе добавление первого узла, когда замкнём
	{
		if (it->n0 == cur_node)
		{
			m_bordernodes().push_back(it->n1);
			cur_node = it->n1;
			edges_hash.erase(it);
			it = edges_hash.begin();
			continue;
		}
		if (it->n1 == cur_node)
		{
			m_bordernodes().push_back(it->n0);
			cur_node = it->n0;
			edges_hash.erase(it);
			it = edges_hash.begin();
			continue;
		}
		++it;
	}

	int i = 1;
	//ORIGIAL CODE
	/*

	std::list<CEdge> edges;

	for (size_t i = 0; i < m_elements().size(); i++)
	{
	C2DElement3& el = m_elements()[i];

	CEdge ed0(el.n0, el.n1),
	ed1(el.n1, el.n2),
	ed2(el.n2, el.n0);	// рёбра элемента
	bool flg0 = false, flg1 = false, flg2 = false;

	//
	std::list<CEdge>::iterator it = edges.begin();	// указатель на список граничных рёбер (принадлежит только 1-му эл-ту)
	std::list<CEdge>::iterator it1 = it;

	while (it != edges.end())		//
	{
	if ((*it) == ed0)
	{
	edges.erase(it);
	flg0 = true;
	it = it1;
	}
	else if ((*it) == ed1)
	{
	edges.erase(it);
	flg1 = true;
	it = it1;
	}
	else if ((*it) == ed2)
	{
	edges.erase(it);
	flg2 = true;
	it = it1;
	}
	else				// если нет таких рёбер
	{
	it1 = it;		// добавляем в список
	++it;
	}
	}

	if (!flg0) edges.push_back(ed0);
	if (!flg1) edges.push_back(ed1);
	if (!flg2) edges.push_back(ed2);
	}

	std::list<CEdge>::iterator it = edges.begin();
	
	m_bordernodes().clear();
	*/
	
	//COMMENTED ORIGINAL CODE
	/*
	while (it != edges.end())
	{
	AddBorderNode(it->n0);
	AddBorderNode(it->n1);
	it++;
	}


	it = edges.begin();
	m_bordernodes().RemoveAll();
	*/
	
	//ORIGIAL CODE
	/*
	// добавляем узлы с учетом их порядка на границе (против часовой стрелки)
	// т.к. контур замкнут
	m_bordernodes().push_back(it->n0);
	m_bordernodes().push_back(it->n1);

	ptrdiff_t cur_node = it->n1;
	edges.erase(it);

	it = edges.begin();

	while (edges.size() > 1) // чтобы отсечь второе добавление нулевого узла
	{
		if (it->n0 == cur_node)
		{
			m_bordernodes().push_back(it->n1);
			cur_node = it->n1;
			edges.erase(it);
			it = edges.begin();
			continue;
		}
		++it;
	}
	*/
}


void C2DMesh3::CalcFormCoeffs()
{
	//LOGGER.Init(CString("..\\..\\Logs\\C2DMesh.cpp_CalcFromCoeff.txt"));
	//this->WriteToLog();

	for (size_t e = 0; e < m_elements().size(); e++)
	{
		C2DElement3& es = m_elements()[e];
		Math::C2DPoint& n0 = m_nodes()[es.n0]; // 0 - i
		Math::C2DPoint& n1 = m_nodes()[es.n1]; // 1 - j
		Math::C2DPoint& n2 = m_nodes()[es.n2]; // 2 - k

		DBL B[3];

		B[0] = n1.y - n2.y;//nodes[j].y - nodes[k].y;
		B[1] = n2.y - n0.y;//nodes[k].y - nodes[i].y;
		B[2] = n0.y - n1.y;//nodes[i].y - nodes[j].y;

						   /*
						   C[0] = nodes[k].x - nodes[j].x;
						   C[1] = nodes[i].x - nodes[k].x;
						   C[2] = nodes[j].x - nodes[i].x;
						   */

						   // координаты центра
						   //a = (n0.x + n1.x + n2.x) / 3.0;
						   //b = (n0.y + n1.y + n2.y) / 3.0;
		GetEField(e, eFields::cm_x) = (n0.x + n1.x + n2.x) / 3.0;//(nodes[i].x + nodes[j].x + nodes[k].x) / 3.0;
		GetEField(e, eFields::cm_y) = (n0.y + n1.y + n2.y) / 3.0;//(nodes[i].y + nodes[j].y + nodes[k].y) / 3.0;

																 // площадь элемента
		GetEField(e, eFields::sqr) = abs(n0.x*B[0] + n1.x*B[1] + n2.x*B[2])*0.5; //(nodes[i].x * B[0] + nodes[j].x * B[1] + nodes[k].x * B[2]) * 0.5;
	}

}

void C2DMesh3::DrawGL(GLParam &parameter)
{
	if (parameter.GetField()>0) {
		if (m_elements().size() == 0) return;

		glColor3f(0.0, 0.0, 0.0);

		/* C2DMesh enum
		cm_x =  (size_t)0x00001, // координата x центра
		cm_y =  (size_t)0x00002, // координата y центра

		def   = (size_t)0x00004, // деформация
		def_s = (size_t)0x00008, // скорость деформаций

		str   = (size_t)0x00010, // напряжение
		avg_d = (size_t)0x00020, // средние деформации

		eps =   (size_t)0x00040, // модуль Юнга

		mju =   (size_t)0x00100, // коэффициент вязкости
		sqr  =  (size_t)0x01000, // площадь элемента

		int_d  = (size_t)0x02000, // интенсивность деформаций
		int_ds = (size_t)0x04000, // интенсивность скорости деформаций
		int_s  = (size_t)0x08000, // интенсивность напряжений
		int_ss = (size_t)0x10000, //

		hyd_p  = (size_t)0x20000, // накопленное гидростатическое давление
		smo_d  = (size_t)0x40000, // сглаженные средние деформации
		//*/
		int field_num = parameter.GetField();//int_d;	//сменил имя на field_num чтобы не теряться при поиске по "field"
		parameter.m_PolygonBuilder.Draw(*this, field_num, *parameter.GetScale());

		if (parameter.m_ScaleVisualizer.IsScaleVisible()) {
			parameter.m_ScaleVisualizer.DrawScale(*parameter.GetScale());
		}

		glColor4f(0, 0, 0, 1.);
	}
	else if (parameter.GetField() == 0) {
		glColor3f(0.0, 0.0, 0.0);
		for (size_t i = 0; i<m_elements.GetSize(); i++) {
			glBegin(GL_LINE_LOOP);
			glVertex2d(m_nodes[m_elements[i].n0].m_x, m_nodes[m_elements[i].n0].m_y);
			glVertex2d(m_nodes[m_elements[i].n1].m_x, m_nodes[m_elements[i].n1].m_y);
			glVertex2d(m_nodes[m_elements[i].n2].m_x, m_nodes[m_elements[i].n2].m_y);
			glEnd();
		}
	}

	//рисуем узлы КЭ сетки
	glColor3f(0.5, 0.5, 0.5);
	glPointSize(3);
	glBegin(GL_POINTS);
	for (size_t i = 0; i<m_nodes().size(); i++) {
		glVertex2d(m_nodes[i].m_x, m_nodes[i].m_y);
	}
	glEnd();
	glPointSize(1);
	
	//else {
		//glColor3f(0.0, 0.0, 0.0);
		//glBegin(GL_LINE_LOOP);
		//for (size_t i = 0; i<m_bordernodes.GetSize(); i++) {
		//	glVertex2d(m_nodes[m_bordernodes[i]].m_x, m_nodes[m_bordernodes[i]].m_y);
		//}
		//glEnd();
	//} 
	
}

//------------------------------------------------------------------------------------------------------------
void C2DMesh3::DrawGL3D(GLParam &parameter)
{
	if (parameter.GetField()>0) {
		if (m_elements().size() == 0) return;


		//// C2DMesh enum
		////cm_x =  (size_t)0x00001, // координата x центра
		////cm_y =  (size_t)0x00002, // координата y центра

		////def   = (size_t)0x00004, // деформация
		////def_s = (size_t)0x00008, // скорость деформаций

		////str   = (size_t)0x00010, // напряжение
		////avg_d = (size_t)0x00020, // средние деформации

		////eps =   (size_t)0x00040, // модуль Юнга

		////mju =   (size_t)0x00100, // коэффициент вязкости
		////sqr  =  (size_t)0x01000, // площадь элемента

		////int_d  = (size_t)0x02000, // интенсивность деформаций
		////int_ds = (size_t)0x04000, // интенсивность скорости деформаций
		////int_s  = (size_t)0x08000, // интенсивность напряжений
		////int_ss = (size_t)0x10000, // 

		////hyd_p  = (size_t)0x20000, // накопленное гидростатическое давление
		////smo_d  = (size_t)0x40000, // сглаженные средние деформации
		////
		int field_num = parameter.GetField();//int_d;	//сменил имя на field_num чтобы не теряться при поиске по "field"
		parameter.m_PolygonBuilder.Draw(*this, field_num, *parameter.GetScale());

	}

}

void C2DMesh3::ReorderNodes(const std::vector<size_t>& neworder)
{
	C2DMeshInterface::ReorderNodes(neworder);

	std::vector<size_t> temp(neworder), temp1(neworder);

	for (size_t i = 0; i < neworder.size(); i++) {
		temp[neworder[i]] = i;
	}

	for (size_t i = 0; i < neworder.size(); i++) {
		temp1[i] = temp[temp[i]];
	}


	// элементы
	IO::CSimpleArray<C2DElement3, C2DEL3ARRAY> tmp1;
	tmp1().resize(m_elements().size());
	for (size_t i = 0; i < m_elements().size(); i++)
	{
		m_elements()[i]._0 = temp[m_elements()[i]._0];
		m_elements()[i]._1 = temp[m_elements()[i]._1];
		m_elements()[i]._2 = temp[m_elements()[i]._2];

	}

}

std::vector< std::vector<size_t> > C2DMesh3::GetNodeLinks() const
{
	CMatrix<size_t, short> adj;
	adj.resize(m_nodes().size(), m_nodes().size(), 0);

	std::vector< std::vector<size_t> > ret;
	//ret.resize(m_nodes().GetCount());

	for (size_t i = 0; i < m_elements().size(); i++)
	{
		const C2DElement3& e = m_elements()[i];
		adj.cell(e._0, e._0) =
			adj.cell(e._1, e._1) =
			adj.cell(e._2, e._2) =

			adj.cell(e._0, e._1) =
			adj.cell(e._1, e._0) =

			adj.cell(e._0, e._2) =
			adj.cell(e._2, e._0) =

			adj.cell(e._1, e._2) =
			adj.cell(e._2, e._1) = 1;
	}
	adj.WriteToLog(true);

	for (size_t i = 0; i < adj.rows(); i++) {

		std::vector<size_t> temp;

		temp.clear();
		for (size_t j = 0; j < adj.cols(); j++) {
			if ((i != j) && (adj.cell(i, j))) temp.push_back(j);
		}
		ret.push_back(temp);
	}

	return ret;
}

//! 
void C2DMesh3::WriteToLog() const
{
	DLOG(CString(_T("C2DMesh3<Type:")) + AllToString(GetType()) + CString(_T(">")), log_info);

	DLOG(CString(_T("Node count: ")) + AllToString(m_nodes().size()) + CString(_T(" Element count: ")) + AllToString(m_elements().size()), log_info);

	for (size_t i = 0; i < m_nodes().size(); i++) {
		DLOG(CString(_T("node<")) + AllToString(i) + CString(_T(">[")) + AllToString(i) + CString(_T("]: ")) + AllToString(m_nodes()[i].x) + CString(_T(", ")) + AllToString(m_nodes()[i].y), log_info);
	}

	for (size_t i = 0; i < m_elements().size(); i++) {
		DLOG(CString(_T("element[")) + AllToString(i) + CString(_T("]: ")) + AllToString(m_elements()[i]._0) + CString(_T(", ")) + AllToString(m_elements()[i]._1) + CString(_T(", ")) + AllToString(m_elements()[i]._2), log_info);
	}
}



//! Создание горизонтальных рисок (использовалось)
void C2DMesh3::GenerateMarks(){
/*
for (int i=0; i < m_bordernodes().GetCount(); i++){
//В крайне правое положение
if (m_nodes()[m_bordernodes()[i]].m_x == 5.0 && m_nodes()[m_bordernodes()[i]].m_y == 0.0){
m_bordernodes_marks().Add(m_bordernodes()[i]);
}

if (m_nodes()[m_bordernodes()[i]].m_x == 5.0 && abs(m_nodes()[m_bordernodes()[i]].m_y - 5.0) < 0.1){
m_bordernodes_marks().Add(m_bordernodes()[i]);
}

if (m_nodes()[m_bordernodes()[i]].m_x == 5.0 && abs(m_nodes()[m_bordernodes()[i]].m_y - 10.0) < 0.1){
m_bordernodes_marks().Add(m_bordernodes()[i]);
}

if (m_nodes()[m_bordernodes()[i]].m_x == 5.0 && abs(m_nodes()[m_bordernodes()[i]].m_y - 15.0) < 0.1){
m_bordernodes_marks().Add(m_bordernodes()[i]);
}
}
//*/
} 