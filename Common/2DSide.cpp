#include "stdafx.h"
#include "2DSide.h"

//------------------ C2DVertex ----------------------------

/*
https://www-01.ibm.com/support/knowledgecenter/SSLTBW_2.1.0/com.ibm.zos.v2r1.cbclx01/cplr388.htm
Or a constructor can have an initializer list within the definition but prior to the constructor body:
complx(double r, double i = 0) : re(r), im(i) { ... }
*/

C2DVertex::C2DVertex(double x, double y, double value, bool lvl) : Math::C2DPoint(x, y) {
	level = lvl;
	m_value = value;
};

C2DVertex::C2DVertex(const Math::C2DPoint& point, double value, bool lvl) : Math::C2DPoint(point) {
	level = lvl;
	m_value = value;
};

C2DVertex::C2DVertex(const C2DVertex& vertex) : Math::C2DPoint(static_cast<const Math::C2DPoint&>(vertex)) {
	level = vertex.level;
	m_value = vertex.m_value;
};

C2DVertex& C2DVertex::operator = (const C2DVertex& p)
{
	m_x = p.m_x;
	m_y = p.m_y;
	m_value = p.m_value;
	level = p.level;
	return *this;
}

//----------------------	C2DSide	------------------------

C2DSide::C2DSide(C2DVertex startVertex, C2DVertex endVertex, const C2DScale& scale)
{
	if (startVertex.m_value > endVertex.m_value)
	{
		m_Vertices.push_back(endVertex);
		m_Vertices.push_back(startVertex);
	}
	else
	{
		m_Vertices.push_back(startVertex);
		m_Vertices.push_back(endVertex);
	}

	m_Splitted = false;
	Split(scale);
}

//Разбиение стороны для шкалы
void C2DSide::Split(const C2DScale& scale) {
	const C2DVertex startVertex = GetStartVertex();
	const C2DVertex endVertex = GetEndVertex();

	// Если значения в начальной и конечной точках
	// равны, то разбивать сторону не надо		
	if (m_Splitted || fabs(startVertex.m_value - endVertex.m_value) < EPS) {
		return;
	}

	// Координаты вектора
	C2DVertex sideVector = endVertex - startVertex;

	// Разность значений в конечном и начальном узлах
	double difference = endVertex.m_value - startVertex.m_value;

	// Проходимся по отметкам шкалы
	for (size_t i = 0; i < scale.GetCount(); i++){

		double scaleValue = scale.GetValue(i);

		// Если значение отметки лежит между значениями граничных вершин,
		// то добавляем новую вершину со значением отметки
		if (startVertex.m_value < scaleValue && endVertex.m_value > scaleValue){
			
			// Процент значения отметки от конечного значения
			double percent = (scaleValue - startVertex.m_value) / difference;

			// Новая вершина = Координаты начала стороны + процент * координаты вектора стороны
			C2DVertex vertex = C2DVertex(startVertex + sideVector * percent, scaleValue);
			m_Vertices.insert(m_Vertices.end() - 1, vertex);
		}
	}

	m_Splitted = true;
}

// Возвращает вершину с таким же значением, но другими координатами
C2DVertex* C2DSide::GetMatchValueVertex(const C2DVertex& vertex) const {
	
	for (auto it = m_Vertices.begin(); it != m_Vertices.end(); ++it){
		// Если совпадают значения, но не совпадают координаты
		if (it->m_value == vertex.m_value && (it->x != vertex.x || it->y != vertex.y)) {
			return const_cast<C2DVertex*>(it._Ptr);
		}
	}

	return nullptr;
}

