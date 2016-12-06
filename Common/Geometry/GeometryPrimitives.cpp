#include "includes.h"
#include "Utils.h"
#include "GeometryUtils.h"
#include "GeometryPrimitives.h"

namespace FemLibrary {

	size_t Vector2D::GetHash(const Vector2D& p)
	{
		return Utils::hash_combine(std::hash<double>()(p.X), std::hash<double>()(p.Y));
	}

	size_t Vector3D::GetHash(const Vector3D& p)
	{
		auto hash = Utils::hash_combine(std::hash<double>()(p.X), std::hash<double>()(p.Y));
		return Utils::hash_combine(hash, std::hash<double>()(p.Z));
	}

	const Vector2D& SimpleTriangle::NextV(int ind) const
	{
		return V[Utils::NextInd(ind, 3)];
	}

	const Vector2D& SimpleTriangle::PrevV(int ind) const
	{
		return V[Utils::PrevInd(ind, 3)];
	}

	double SimpleTriangle::Cos(int ind) const
	{
		auto p2 = E[Utils::PrevInd(ind, 3)].Vector().LengthSquared();
		auto n2 = E[Utils::NextInd(ind, 3)].Vector().LengthSquared();
		return (p2 + n2 - E[ind].Vector().LengthSquared()) / (2 * sqrt(p2*n2));
	}

	void Triangle::MoveVertex(int vertexInd, const Vector2D& newVal)
	{
		V[vertexInd].X = newVal.X;
		V[vertexInd].Y = newVal.Y;
		E[Utils::PrevInd(vertexInd, 3)].S =
			E[Utils::NextInd(vertexInd, 3)].E = V[vertexInd];
		Recalc();
	}

	std::vector<std::pair<const Triangle*, int>> Triangle::GetStar(int ind) const
	{
		auto center = V[ind];
		std::vector<std::pair<const Triangle*, int>> star = { std::make_pair(this, ind) };
		auto n = Neighbour[Utils::NextInd(ind, 3)];
		auto fromSide = NeighbourInd[Utils::NextInd(ind, 3)];
		auto secondPass = false;
	Again:
		while (n != nullptr && n != this)
		{
			auto next = Utils::NextInd(fromSide, 3);
			auto opps = Utils::PrevInd(fromSide, 3);
			if (n->E[fromSide].E != center) std::swap(next, opps); // todo: remove unnecessary ifs
			star.emplace_back(std::make_pair(n, opps));
			assert(n->E[fromSide].Contains(center) && n->E[next].Contains(center) && !n->E[opps].Contains(center));
			fromSide = n->NeighbourInd[next];
			n = n->Neighbour[next];
		}
		if (n == this || secondPass) return star;
		fromSide = NeighbourInd[Utils::PrevInd(ind, 3)];
		n = Neighbour[Utils::PrevInd(ind, 3)];
		secondPass = true;
		goto Again;
	}

	void Triangle::Recalc()
	{
		auto SignedArea = GeometryUtils::TriangleOrientedArea(V[0], V[1], V[2]);
		Orientation = GeometryUtils::Orientation(SignedArea);
		if (Orientation == PointsOrder::Collinear)
			throw std::exception("Triangle.Ctor: collinear points");

		// в некоторых случаях возможно появление почти вырожденных треугольников на границе
		// но они будут удалены, как внешние треугольники
		//assert(abs(cosinus) < 1);
		ObtuseVertexInd = -1;
		auto len20 = E[0].Vector().LengthSquared();
		auto len21 = E[1].Vector().LengthSquared();
		auto len22 = E[2].Vector().LengthSquared(); // прямой не считается тупым !!!
		if (len20 > len21 + len22) ObtuseVertexInd = 0;
		if (len21 > len20 + len22) ObtuseVertexInd = 1;
		if (len22 > len20 + len21) ObtuseVertexInd = 2;

		AlphaQuality = 4 * sqrt(3) * SignedArea / (E[0].Vector().LengthSquared() + E[1].Vector().LengthSquared() + E[2].Vector().LengthSquared());

		auto res = GeometryUtils::Circumcenter(V[0], V[1], V[2], Circumcenter);
		assert(res);
		// при построении выпуклой триангуляции методом D&C возможен нулевой центр

		//Assert(Circumcenter != null); 
		CircleRadius = (V[0] - Circumcenter).Length();

		//var shortestLength = E[0].Vector().Length(); ShortestEdgeInd = 0;
		//if (E[1].Vector.Length < shortestLength) { shortestLength = E[1].Vector.Length; ShortestEdgeInd = 1; }
		//if (E[2].Vector.Length < shortestLength) { shortestLength = E[2].Vector.Length; ShortestEdgeInd = 2; }
		//RadiusEdgeRatio = CircleRadius / shortestLength;
	}
}