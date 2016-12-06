#pragma once
#include <vector>
#include "FemLibrary.h"

namespace FemLibrary {

	class Vector2D;
	class LineSeg;
	struct FragmentGeometry;


	class FEMLIBRARY_API NonIntersectedCurvesPool
	{
	public:
		std::vector<LineSeg> Lines;
		void AddLine(LineSeg line);
	};

	enum class PointsOrder;

	class FEMLIBRARY_API GeometryUtils
	{
	public:
		static std::vector<FragmentGeometry> FindContours(const std::vector<LineSeg>& lines);
		static bool IsPointInsideContour(const Vector2D& p, const std::vector<Vector2D>& contour, bool isBoundaryInside = true);
		static double ContourSignedArea(const std::vector<Vector2D>& points);
		static std::tuple<std::vector<LineSeg>, std::vector<LineSeg>> DivideLines(const LineSeg& line1, const LineSeg& line2);

		// robust predicates
		static bool InCircle(const Vector2D& a, const Vector2D& b, const Vector2D& c, const Vector2D& d, bool clockwise = false);
		static double TriangleOrientedArea(const Vector2D& a, const Vector2D& b, const Vector2D& c);
		static double TriangleAlphaQuality(const Vector2D& a, const Vector2D& b, const Vector2D& c);
		static PointsOrder Orientation(const Vector2D& a, const Vector2D& b, const Vector2D& c);
		static PointsOrder Orientation(double triangleOrientedArea);
		static bool Circumcenter(const Vector2D& p1, const Vector2D& p2, const Vector2D& p3, Vector2D& center);
		static bool Circumcenter(const Vector2D& p1, const Vector2D& p2, double radius, bool clockwise, bool largeArc, Vector2D& center);

		static bool IntersectionOfNonCollinearLines(const Vector2D& start1, const Vector2D& end1, const Vector2D& start2, const Vector2D& end2, Vector2D& intersection);
	};

}