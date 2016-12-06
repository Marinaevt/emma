#pragma once
#include "2DSide.h"
#include <gl/GL.h>
#include "2DScale.h"

class C2DMesh3;
class C2DElement3;

class C2DPolygonBuilder
{
public:
	C2DPolygonBuilder();

	C2DPolygonBuilder(const C2DPolygonBuilder &builder);

	~C2DPolygonBuilder()
	{
		glDeleteLists(m_DisplayListIndex, 1);
	}

	C2DPolygonBuilder& operator = (const C2DPolygonBuilder& builder);

	void Draw(C2DMesh3 &mesh, int fieldNum, C2DScale &scale);	// in Mesh.cpp
	void ClearDisplayList();

	void UpdateBorderNodes(C2DMesh3 &mesh, int fieldNum, C2DScale &scale, std::vector<C2DVertex>& upd_nodes);
	C2DSide C2DPolygonBuilder::GetBorderSide(C2DMesh3 &mesh, int n1, int n2, int fieldNum, const C2DScale &scale);
	void SetScaleBounds(C2DMesh3 &mesh, int fieldNum, C2DScale &scale);

	void MinMaxOut();

protected:
	GLuint m_DisplayListIndex;
	double m_MaxValue;
	double m_MinValue;
	double m_Norm;

	void DetermineMinMax(C2DMesh3 &mesh, int field_num);
	void DrawElement(C2DMesh3 &mesh, C2DElement3 &e, int fieldNum, const C2DScale &scale);		

	bool IsSimple(const C2DSide &mainSide) const
	{
		return mainSide.m_Vertices.size() == 2;
	}

	void DrawTriangle(const C2DSide &mainSide, const std::vector<C2DSide> &otherSides, const C2DScale &scale) const;
	void DrawPolygons(const C2DSide &mainSide, std::vector<C2DSide> &otherSides, const C2DScale &scale) const;
};