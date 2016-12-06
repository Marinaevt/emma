#include "stdafx.h"
#include "2DPolygonBuilder.h"
#include "2DMesh.h"

C2DPolygonBuilder::C2DPolygonBuilder()
{
	m_DisplayListIndex = 0;
	m_MaxValue = 0;
	m_MinValue = 0;
	m_Norm = 1;
}

C2DPolygonBuilder::C2DPolygonBuilder(const C2DPolygonBuilder &builder)
{
	m_DisplayListIndex = builder.m_DisplayListIndex;
	m_MaxValue = builder.m_MaxValue;
	m_MinValue = builder.m_MinValue;
	m_Norm = builder.m_Norm;
}

C2DPolygonBuilder& C2DPolygonBuilder::operator = (const C2DPolygonBuilder& builder){
		
		m_DisplayListIndex = builder.m_DisplayListIndex;
		m_MaxValue = builder.m_MaxValue;
		m_MinValue = builder.m_MinValue;
		m_Norm = builder.m_Norm;
		return *this;
}

void C2DPolygonBuilder::ClearDisplayList(){
	if (m_DisplayListIndex != 0)
		{
			glDeleteLists(m_DisplayListIndex, 1);
			m_DisplayListIndex = 0;
		}
}

void C2DPolygonBuilder::Draw(C2DMesh3 &mesh, int fieldNum, C2DScale &scale)
{
	if (scale.IsInit() == false)
		return;

	if (m_DisplayListIndex == 0)
	{
		DetermineMinMax(mesh, fieldNum);
		scale.SetBounds(m_MinValue, m_MaxValue);

		double diff = m_MaxValue - m_MinValue;

		if (abs(diff) > EPS)	{
			m_Norm = 1 / diff;
		}else{
			m_Norm = 1;
		}

		m_DisplayListIndex = glGenLists(1);
		glNewList(m_DisplayListIndex, GL_COMPILE);

			for (size_t i = 0; i < mesh.m_elements().size(); i++)
			{
				C2DElement3& e = mesh.m_elements()[i];
				DrawElement(mesh, e, fieldNum, scale);
			}

		glEndList();
	}

	glCallList(m_DisplayListIndex);
}

void C2DPolygonBuilder::DetermineMinMax(C2DMesh3 &mesh, int field_num)
{
	m_MaxValue = m_MinValue = mesh.GetNField(mesh.m_elements()[0].n0, field_num);

	for (size_t i = 0; i < mesh.m_elements().size(); i++)
	{
		C2DElement3& e = mesh.m_elements()[i];

		double n0Value = mesh.GetNField(e.n0, field_num);
		double n1Value = mesh.GetNField(e.n1, field_num);
		double n2Value = mesh.GetNField(e.n2, field_num);

		if (m_MaxValue < n0Value)
			m_MaxValue = n0Value;
		if (m_MinValue > n0Value)
			m_MinValue = n0Value;

		if (m_MaxValue < n1Value)
			m_MaxValue = n1Value;
		if (m_MinValue > n1Value)
			m_MinValue = n1Value;

		if (m_MaxValue < n2Value)
			m_MaxValue = n2Value;
		if (m_MinValue > n2Value)
			m_MinValue = n2Value;
	}	
}

void C2DPolygonBuilder::DrawElement(C2DMesh3 &mesh, C2DElement3& e, int fieldNum, const C2DScale &scale)
{		
	std::vector<C2DSide> sides;

	for (UNINT i = 0; i < 3; i++)
	{
		C2DVertex firstVertex = C2DVertex(mesh.m_nodes()[e.n[i]]);
		firstVertex.m_value = (mesh.GetNField(e.n[i], fieldNum) - m_MinValue) * m_Norm;

		C2DVertex secondVertex = C2DVertex(mesh.m_nodes()[e.n[i == 2 ? 0 : i + 1]]);
		secondVertex.m_value = (mesh.GetNField(e.n[i == 2 ? 0 : i + 1], fieldNum) - m_MinValue) * m_Norm;
	
		sides.push_back(C2DSide(firstVertex, secondVertex, scale));
	}

	// Нахождение главной стороны
	size_t longestIndex = 0;
	for (size_t i = 1; i < 3; i++)
	{
		if (sides[longestIndex].GetEndVertex().m_value - sides[longestIndex].GetStartVertex().m_value
			< sides[i].GetEndVertex().m_value - sides[i].GetStartVertex().m_value)
			longestIndex = i;
	}

	// Удаление главной стороны из списка
	const C2DSide mainSide = sides[longestIndex];
	sides.erase(sides.begin() + longestIndex);

	DrawPolygons(mainSide, sides, scale);	
}

void C2DPolygonBuilder::DrawTriangle(const C2DSide &mainSide, const std::vector<C2DSide> &otherSides, const C2DScale &scale) const{

	const C2DVertex& firstSideStartVertex = mainSide.GetStartVertex();
	const C2DVertex& firstSideEndVertex = mainSide.GetEndVertex();
	const C2DVertex& secondSideStartVertex = otherSides[0].GetStartVertex();
	const C2DVertex& secondSideEndVertex = otherSides[0].GetEndVertex();
	C2DColor color = scale.GetColor(firstSideStartVertex.m_value);
	glColor3d(color.GetR(), color.GetG(), color.GetB());

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glBegin(GL_TRIANGLES);
		
		glVertex2d(firstSideStartVertex.x, firstSideStartVertex.y);

		if (!scale.IsDescrete()) {
			color = scale.GetColor(firstSideEndVertex.m_value);
			glColor3d(color.GetR(), color.GetG(), color.GetB());
		}
				
		glVertex2d(firstSideEndVertex.x, firstSideEndVertex.y);				

		// Если начальная вершина второй стороны не является ни одной из вершин первой
		if (secondSideStartVertex != firstSideStartVertex && secondSideStartVertex != firstSideEndVertex)
		{
			if (!scale.IsDescrete()) {
				color = scale.GetColor(secondSideStartVertex.m_value);
				glColor3d(color.GetR(), color.GetG(), color.GetB());
			}
					
			glVertex2d(secondSideStartVertex.x, secondSideStartVertex.y);					
		}
		else
		{
			if (!scale.IsDescrete()) {
				color = scale.GetColor(secondSideEndVertex.m_value);
				glColor3d(color.GetR(), color.GetG(), color.GetB());
			}
					
			glVertex2d(secondSideEndVertex.x, secondSideEndVertex.y);					
		}

	glEnd();
}

void C2DPolygonBuilder::DrawPolygons(const C2DSide &mainSide, std::vector<C2DSide> &otherSides, const C2DScale &scale) const{

	if (otherSides.size() != 2){
			//throw;	//Catch-a нет
			return;
		}

		if (IsSimple(mainSide))
		{
			DrawTriangle(mainSide, otherSides, scale);
			return;
		}

		// Выбираем такой конец второй стороны, который не совпадает
		// с концами главной стороны. Значение этой вешины может не совпадать
		// ни с одним значением шкалы, и тогда будет образован
		// пятиугольник с этой вершиной.
		const C2DVertex& pentagonVetex = otherSides[0].GetStartVertex() != mainSide.GetStartVertex() && otherSides[0].GetStartVertex() != mainSide.GetEndVertex() ?
			otherSides[0].GetStartVertex() : otherSides[0].GetEndVertex();

		for (size_t i = 0; i < mainSide.m_Vertices.size() - 1; i++)
		{
			C2DVertex* contourLineSecondVertex = nullptr;
			C2DVertex *oppositePoint = nullptr;
			C2DVertex curVertex = mainSide.m_Vertices[i];
			C2DVertex nextVertex = mainSide.m_Vertices[i+1];
			C2DColor color = scale.GetColor(curVertex.m_value);
			
			glColor3d(color.GetR(), color.GetG(), color.GetB());

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

			glBegin(GL_POLYGON);

			// Первая вершина
			glVertex2d(curVertex.x, curVertex.y);					

				// Нахождение вершины с таким же значением на других сторонах для первой вершины
				for (auto it = otherSides.begin(); it != otherSides.end(); ++it) {						
					if ((oppositePoint = it->GetMatchValueVertex(curVertex)) != nullptr)
					{
						// Конечная вершина линии уровня
						contourLineSecondVertex = oppositePoint;

						glVertex2d(oppositePoint->x, oppositePoint->y);
						break;
					}
				}

				// Проверка на пятиугольник
				if (pentagonVetex.m_value > curVertex.m_value &&
					pentagonVetex.m_value < nextVertex.m_value)
				{
					if (!scale.IsDescrete())
					{
						color = scale.GetColor(pentagonVetex.m_value);
						glColor3d(color.GetR(), color.GetG(), color.GetB());							
					}
						
					glVertex2d(pentagonVetex.x, pentagonVetex.y);
				}

				// Нахождение вершины с таким же значением на других сторонах для второй вершины
				for (auto it = otherSides.begin(); it != otherSides.end(); ++it){		
					if ((oppositePoint = it->GetMatchValueVertex(nextVertex)) != nullptr)
					{
						if (!scale.IsDescrete())
						{
							color = scale.GetColor(nextVertex.m_value);
							glColor3d(color.GetR(), color.GetG(), color.GetB());
						}
							
						glVertex2d(oppositePoint->x, oppositePoint->y);
						break;
					}
				}

				if (!scale.IsDescrete()){

					color = scale.GetColor(nextVertex.m_value);
					glColor3d(color.GetR(), color.GetG(), color.GetB());
				}
					
			// Вторая вершина
			glVertex2d(nextVertex.x, nextVertex.y);					

			glEnd();

			// Отрисовка линий уровня
			if (scale.IsDrawLines() && contourLineSecondVertex != nullptr){
				
				glColor3d(0, 0, 0);
				glBegin(GL_LINES);

					glVertex2d(curVertex.x, curVertex.y);
					glVertex2d(contourLineSecondVertex->x, contourLineSecondVertex->y);

				glEnd();
			}
		}
}


C2DSide C2DPolygonBuilder::GetBorderSide(C2DMesh3 &mesh, int n1, int n2, int fieldNum, const C2DScale &scale)
{
	C2DVertex firstVertex = C2DVertex(mesh.m_nodes()[n1], (mesh.GetNField(n1, fieldNum) - m_MinValue) * m_Norm);
	C2DVertex secondVertex = C2DVertex(mesh.m_nodes()[n2], (mesh.GetNField(n2, fieldNum) - m_MinValue) * m_Norm);

	return C2DSide(firstVertex, secondVertex, scale);
}

void C2DPolygonBuilder::UpdateBorderNodes(C2DMesh3 &mesh, int fieldNum, C2DScale &scale, std::vector<C2DVertex>& upd_nodes)
{
	DetermineMinMax(mesh, fieldNum);

	double diff = m_MaxValue - m_MinValue;

	if (abs(diff) > EPS) {
		m_Norm = 1 / diff;
	}
	else {
		m_Norm = 1;
	}

	C2DSide side;
	for (size_t i = 0; i < mesh.m_bordernodes.GetSize(); i++)
	{
		side = GetBorderSide(mesh, mesh.m_bordernodes()[i], mesh.m_bordernodes()[(i + 1) % mesh.m_bordernodes.GetSize()], fieldNum, scale);
		for (size_t j = 0; j < side.m_Vertices.size() - 1; j++)
		{
			if (j > 0)
				side.m_Vertices[j].level = true;

			upd_nodes.push_back(side.m_Vertices[j]);
		}
		//.insert(upd_nodes.end(), side.m_Vertices.begin(), side.m_Vertices.end() - 1);
	}


	upd_nodes.push_back(side.GetEndVertex());
	int i = 0;
	i++;
}

void C2DPolygonBuilder::SetScaleBounds(C2DMesh3 &mesh, int fieldNum, C2DScale &scale)
{
	DetermineMinMax(mesh, fieldNum);
	scale.SetBounds(m_MinValue, m_MaxValue);
}

void C2DPolygonBuilder::MinMaxOut()
{

}