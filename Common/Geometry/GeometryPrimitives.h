#pragma once
#include <algorithm>
#include <vector>
#include "FemLibrary.h"


namespace FemLibrary {

	enum class PointsOrder { Counterclockwise, Clockwise, Collinear };

	class Vector2D
	{
	public: double X, Y;
			Vector2D operator -()                  const { return Vector2D(-X, -Y); }
			Vector2D operator -(const Vector2D& p) const { return Vector2D(X - p.X, Y - p.Y); }
			Vector2D operator +(const Vector2D& p) const { return Vector2D(X + p.X, Y + p.Y); }
			double   operator *(const Vector2D& p) const { return X * p.X + Y * p.Y; }
			Vector2D& operator+=(const Vector2D& rhs) { X += rhs.X; Y += rhs.Y; return *this; }
			Vector2D& operator-=(const Vector2D& rhs) { X -= rhs.X; Y -= rhs.Y; return *this; }
			bool operator==(const Vector2D& rhs) const { return X == rhs.X && Y == rhs.Y; }
			bool operator!=(const Vector2D& rhs) const { return !this->operator==(rhs); }
			Vector2D operator *(double k) const { return Vector2D(X * k, Y * k); }
			Vector2D operator /(double k) const { return Vector2D(X / k, Y / k); }

			Vector2D() :X(std::numeric_limits<double>::quiet_NaN()), Y(std::numeric_limits<double>::quiet_NaN()) {}
			Vector2D(double x, double y) :X(x), Y(y) {}
			double LengthSquared() const { return X * X + Y * Y; }
			double Length() const { return sqrt(LengthSquared()); }
			bool IsZero() const { return X == .0 && Y == .0; }

			static double CrossProduct(const Vector2D& p1, const Vector2D& p2) { return p1.X * p2.Y - p1.Y * p2.X; }
			Vector2D Norm() const { return Vector2D(-Y, X); }

			static size_t GetHash(const Vector2D& p);

			std::string ToString() const { return std::to_string(X) + " " + std::to_string(Y); }
			const double* ToArray() const { return &X; } // may be dangerous
	};
	inline Vector2D operator*(double k, const Vector2D& v) { return v*k; }
}

namespace std {
	template <> struct hash<FemLibrary::Vector2D> {
		size_t operator()(const FemLibrary::Vector2D& p) const { return FemLibrary::Vector2D::GetHash(p); }
	};
}

namespace FemLibrary {
	class Vector3D
	{
	public: double X, Y, Z;
			Vector3D operator -()                  const { return Vector3D(-X, -Y, -Z); }
			Vector3D operator -(const Vector3D& p) const { return Vector3D(X - p.X, Y - p.Y, Z - p.Z); }
			Vector3D operator +(const Vector3D& p) const { return Vector3D(X + p.X, Y + p.Y, Z + p.Z); }
			double   operator *(const Vector3D& p) const { return X * p.X + Y * p.Y + Z*p.Z; }
			bool operator==(const Vector3D& rhs) const { return X == rhs.X && Y == rhs.Y && Z == rhs.Z; }
			bool operator!=(const Vector3D& rhs) const { return !this->operator==(rhs); }
			Vector3D operator *(double k) const { return Vector3D(X * k, Y * k, Z * k); }
			Vector3D operator /(double k) const { return Vector3D(X / k, Y / k, Z / k); }
			Vector3D() :X(std::numeric_limits<double>::quiet_NaN()), Y(std::numeric_limits<double>::quiet_NaN()), Z(std::numeric_limits<double>::quiet_NaN()) {}
			Vector3D(double x, double y, double z) :X(x), Y(y), Z(z) {}
			double LengthSquared() const { return X * X + Y * Y + Z * Z; }
			double Length() const { return sqrt(LengthSquared()); }
			bool IsZero() const { return X == .0 && Y == .0 && Z == .0; }
			//static double CrossProduct(const Vector3D& p1, const Vector3D& p2) { return p1.X * p2.Y - p1.Y * p2.X; }
			//Vector3D Norm() const { return Vector3D(-Y, X); }

			static size_t GetHash(const Vector3D& p);

			std::string ToString() const { return std::to_string(X) + " " + std::to_string(Y) + " " + std::to_string(Z); }
			const double* ToArray() const { return &X; } // may be dangerous
	};
	inline Vector3D operator*(double k, const Vector3D& v) { return v*k; }
}

namespace std {
	template <> struct hash<FemLibrary::Vector3D> {
		size_t operator()(const FemLibrary::Vector3D& p) const { return FemLibrary::Vector3D::GetHash(p); }
	};
}

namespace FemLibrary{

	class BoundingBox
	{
	public:
		double MinX = DBL_MAX, MinY = DBL_MAX;
		double MaxX = -DBL_MAX, MaxY = -DBL_MAX;
		Vector2D Min() const { return Vector2D(MinX, MinY); }
		Vector2D Max() const { return Vector2D(MaxX, MaxY); }
		Vector2D MinMax() const { return Vector2D(MinX, MaxY); }
		Vector2D MaxMin() const { return Vector2D(MaxX, MinY); }
		double DeltaX() const { return MaxX - MinX; }
		double DeltaY() const { return MaxY - MinY; }
		void Union(const Vector2D& p)
		{
			MinX = (std::min)(MinX, p.X);
			MinY = (std::min)(MinY, p.Y);
			MaxX = (std::max)(MaxX, p.X);
			MaxY = (std::max)(MaxY, p.Y);
		}
		BoundingBox() :MinX(.0), MinY(.0), MaxX(.0), MaxY(.0) {}
		BoundingBox(const std::initializer_list<Vector2D>& points) { for (auto point : points) Union(point); }
		BoundingBox(const std::vector          <Vector2D>& points) { for (auto point : points) Union(point); }
		//BoundingBox(const std::vector<Vector2D>& points) { for (auto point : points) Union(point); }
		bool IntersectsWith(const BoundingBox& bbox) const { return !(MaxX < bbox.MinX || MinX > bbox.MaxX || MaxY < bbox.MinY || MinY > bbox.MaxY); }
		std::string ToString() const { return "[" + Min().ToString() + " - " + Max().ToString() + "]"; }
		bool Contains(const Vector2D& p) const { return MinX <= p.X && p.X <= MaxX && MinY <= p.Y && p.Y <= MaxY; }
	};

	class LineSeg
	{
	public: Vector2D S, E;
			Vector2D Vector() const { return E - S; }
			double Length() const { return Vector().Length(); }
			LineSeg() {}
			LineSeg(const Vector2D& s, const Vector2D& e) : S(s), E(e) {}
			LineSeg(double s1, double e1, double s2, double e2) : S(Vector2D(s1, e1)), E(Vector2D(s2, e2)) {}
			bool Contains(const Vector2D& p) const { return S == p || E == p; }
			bool Inside(const Vector2D& p)const {
				return Vector2D::CrossProduct(Vector(), p - S) == 0 &&
					(S.X != E.X ? (S.X <= p.X && p.X <= E.X || E.X <= p.X && p.X <= S.X)
						: (S.Y <= p.Y && p.Y <= E.Y || E.Y <= p.Y && p.Y <= S.Y));
			}
			//LineSeg& operator=(const LineSeg& arg) { // Caution: опасный код, но иного способа использовать ссылки внутри класса с копированием нет (?)
			//	if (this == &arg) return *this;
			//	//this->~LineSeg();
			//	new(this) LineSeg(arg);
			//	return *this;
			//}
			Vector2D Middle() const { return (S + E) / 2; }
			BoundingBox GetBoundingBox() const { return BoundingBox({ S, E }); }
			bool operator !=(const LineSeg& rhs) const { return !this->operator==(rhs); };
			bool operator ==(const LineSeg& rhs) const {
				return S == rhs.S && E == rhs.E ||
					S == rhs.E && E == rhs.S;
			}
			static size_t GetHash(const LineSeg& s) { return Vector2D::GetHash(s.S) ^ Vector2D::GetHash(s.E); }
	};
}

namespace std {
	template <> struct hash<FemLibrary::LineSeg> {
		size_t operator()(const FemLibrary::LineSeg& p) const { return FemLibrary::LineSeg::GetHash(p); }
	};
}

namespace FemLibrary{

	class SimpleTriangle
	{
	public:
		LineSeg  E[3];

		Vector2D Centroid;

		Vector2D V[3];
		const Vector2D& V0() const { return V[0]; }
		const Vector2D& V1() const { return V[1]; }
		const Vector2D& V2() const { return V[2]; }

		const Vector2D& NextV(int ind) const;

		const Vector2D& PrevV(int ind) const;

		int VertexInd(const Vector2D& p) const
		{
			if (V[0] == p) return 0;
			if (V[1] == p) return 1;
			if (V[2] == p) return 2;
			return -1;
		}

		SimpleTriangle() {}
		SimpleTriangle(const Vector2D& v0, const Vector2D& v1, const Vector2D& v2)
		{
			Centroid = (v0 + v1 + v2) / 3;		//    1        нумерация ребер и соседей:
			V[0] = v0; E[0] = LineSeg(v1, v2);	//   /  \       E[0]: 1->2
			V[1] = v1; E[1] = LineSeg(v2, v0);	//  /    \      E[1]: 2->0
			V[2] = v2; E[2] = LineSeg(v0, v1);	// 0 ---- 2     E[2]: 0->1
		}

		double Cos(int ind) const;
	};

	class Triangle : public SimpleTriangle
	{
	public:
		Triangle* Neighbour[3] = { nullptr,nullptr,nullptr };
		int       NeighbourInd[3] = { -1, -1, -1 };

		bool IsExternal = false;
		bool isStructureEdge[3] = { false, false, false };
		bool IsStructureEdge(int ind) const {
			return isStructureEdge[ind] ||
				Neighbour[ind] == nullptr ||
				Neighbour[ind]->IsExternal ||
				IsExternal;
		}
		void SetStructureEdge(int ind, bool value = true)
		{
			isStructureEdge[ind] = value;
			auto neighbour = Neighbour[ind];
			if (neighbour != nullptr)
				neighbour->isStructureEdge[NeighbourInd[ind]] = value;
		}

		Triangle(const Vector2D& v0, const Vector2D& v1, const Vector2D& v2) : SimpleTriangle(v0, v1, v2)
		{
			Neighbour[0] = Neighbour[1] = Neighbour[2] = nullptr;
			NeighbourInd[0] = NeighbourInd[1] = NeighbourInd[2] = -1;

			//Cos = new double[3];
			Recalc();
		}
		SimpleTriangle ToSimpleTriangle() { return *static_cast<SimpleTriangle*>(this); }


		//public bool IsValidMove(int ind, Vector2D newPos) = > Orientation == OrderOfPoints(PrevV(ind), newPos, NextV(ind));
		//
		void MoveVertex(int vertexInd, const Vector2D& newVal);

		//public bool IsPointInside(Vector2D p) = > E.All(edge = > Orientation == OrderOfPoints(edge.S, edge.E, p));
		//
		std::vector<std::pair<const Triangle*, int>> GetStar(int ind) const;
		//public override string ToString() = > $"({V0}, {V1}, {V2})";
		//public string LogOut() = > $"( {V0.LogOut()} , {V1.LogOut()} , {V2.LogOut()} )";

		//==============================
		//       Additional fields
		//==============================
	public:
		double AlphaQuality;
		Vector2D Circumcenter;
		double CircleRadius;
		//double RadiusEdgeRatio{ get; private set; }
		//int ShortestEdgeInd{ get; private set; }
		int ObtuseVertexInd;
		//double SignedArea{ get; private set; }
		//double Area = > Math.Abs(SignedArea);

		PointsOrder Orientation;

		//public double VertexSin = > Area * 2;
		//public double[] Cos;
		bool IsObtuse() const { return ObtuseVertexInd >= 0; } // прямой тоже считается тупым !!!

	private:
		void Recalc();
	};

	struct FragmentGeometry
	{
		std::vector<Vector2D> ContourPoints;
		std::vector<FragmentGeometry> Holes;
	};
}