#pragma once
#include <vector>
#include "2DPoint.h"
#include "2DScale.h"

//------------------ C2DVertex ----------------------------
class C2DVertex : public Math::C2DPoint
{
public:

	double m_value;
	bool level; // если эта вершина соответсвует линии уровня

	C2DVertex(double x = 0, double y = 0, double value = 0, bool lvl = false);

	C2DVertex(const Math::C2DPoint& point, double value = 0, bool lvl = false);

	C2DVertex(const C2DVertex& vertex);

	C2DVertex& operator = (const C2DVertex& p);

};
//----------------------------------------------------------

//----------------------	C2DSide	------------------------
class C2DSide
{
protected:
	bool m_Splitted;
public:
	std::vector<C2DVertex> m_Vertices;

public:
	C2DSide(void){
		m_Splitted = false;
	}
	
	C2DSide(C2DVertex startVertex, C2DVertex endVertex, const C2DScale& scale);
	
	~C2DSide(void) {}

	const C2DVertex &GetStartVertex() const
	{
		return m_Vertices[0];
	}

	const C2DVertex &GetEndVertex() const
	{
		return m_Vertices[m_Vertices.size() - 1];
	}

	// Возвращает вершину с таким же значением, но другими координатами
	C2DVertex* GetMatchValueVertex(const C2DVertex& vertex) const;

	//Разбиение стороны для шкалы
	void Split(const C2DScale& scale);

};
//----------------------------------------------------------
