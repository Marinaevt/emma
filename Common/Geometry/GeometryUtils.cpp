#include "includes.h"
#include "Utils.h"
#include "GeometryUtils.h"
#include "GeometryPrimitives.h"
#include "Mesher.h"

namespace FemLibrary
{

	std::tuple<std::vector<LineSeg>, std::vector<LineSeg>> GeometryUtils::DivideLines(const LineSeg& line1, const LineSeg& line2)
	{
		auto list1 = std::vector<LineSeg>(), list2 = std::vector<LineSeg>();
		list1.reserve(3);
		list2.reserve(3);
		const Vector2D* intersection;
		//
		auto dir1 = line1.Vector();
		auto dir2 = line2.Vector();

		// обработка вырожденных случаев
		if (dir1.IsZero())
		{
			if (!dir2.IsZero())
			{
				if (!line2.Inside(line1.S))
					list2.push_back(line2);
				else
				{
					list2.push_back(LineSeg(line1.S, line2.S));
					list2.push_back(LineSeg(line1.S, line2.E));
				}
			}
			return make_tuple(list1, list2);
		}
		if (dir2.IsZero())
		{
			if (!line1.Inside(line2.S))
				list1.push_back(line1);
			else
			{
				list1.push_back(LineSeg(line2.S, line1.S));
				list1.push_back(LineSeg(line2.S, line1.E));
			}
			return make_tuple(list1, list2);
		}
		// line1 совпадает с line2 или их ограничивающие коробки не пересекаются
		if (line1 == line2 || !line1.GetBoundingBox().IntersectsWith(line2.GetBoundingBox())) {
			list1.push_back(line1);
			list2.push_back(line2);
			return make_tuple(list1, list2);
		}

		//считаем уравнения прямых проходящих через отрезки
		auto a1 = -dir1.Y;
		auto b1 = +dir1.X;
		auto d1 = -(a1 * line1.S.X + b1 * line1.S.Y);

		auto a2 = -dir2.Y;
		auto b2 = +dir2.X;
		auto d2 = -(a2 * line2.S.X + b2 * line2.S.Y);

		//подставляем концы отрезков, для выяснения в каких полуплоскотях они
		auto seg1_line2_start = a2 * line1.S.X + b2 * line1.S.Y + d2;
		auto seg1_line2_end = a2 * line1.E.X + b2 * line1.E.Y + d2;
		auto seg2_line1_start = a1 * line2.S.X + b1 * line2.S.Y + d1;
		auto seg2_line1_end = a1 * line2.E.X + b1 * line2.E.Y + d1;
		const double Eps = 1e-8;                                                               // !!! epsilon !!!
		if (abs(seg1_line2_start) < Eps) seg1_line2_start = .0;
		if (abs(seg1_line2_end) < Eps) seg1_line2_end = .0;
		if (abs(seg2_line1_start) < Eps) seg2_line1_start = .0;
		if (abs(seg2_line1_end) < Eps) seg2_line1_end = .0;

		//если концы одного отрезка имеют один знак, значит он в одной полуплоскости и пересечения нет.
		if (seg1_line2_start * seg1_line2_end > 0 ||
			seg2_line1_start * seg2_line1_end > 0)
		{
			list1.push_back(line1);
			list2.push_back(line2);
			return make_tuple(list1, list2);
		}

		auto isCollinear1 = Vector2D::CrossProduct(dir1, dir2) == 0.0;
		auto isCollinear2 = seg1_line2_start == .0 && seg1_line2_end == .0;
		auto isCollinear3 = seg2_line1_start == .0 && seg2_line1_end == .0;
		//if(!(isCollinear3 == isCollinear2 && isCollinear2 == isCollinear1))
		//    throw new Exception(Util.ContextString("Collinearity mismatch", MethodBase.GetCurrentMethod(), line1, line2));
		//Assert(isCollinear3 == isCollinear2 && isCollinear2 == isCollinear1, "Collinearity mismatch");
		if (isCollinear1 || isCollinear2 || isCollinear3)
		{
			auto minP = line2.S, maxP = line2.E;
			auto s = line2.S - line1.S;
			auto e = line2.E - line1.S;
			auto min = abs(dir1.X) < abs(dir1.Y) ? s.Y / dir1.Y : s.X / dir1.X;
			auto max = abs(dir1.X) < abs(dir1.Y) ? e.Y / dir1.Y : e.X / dir1.X;
			if (max < min) {
				std::swap(max, min);
				std::swap(minP, maxP);
			}
			if (minP == line1.S) min = 0.0;
			else if (minP == line1.E) min = 1.0;
			if (maxP == line1.S) max = 0.0;
			else if (maxP == line1.E) max = 1.0;

			if (0 < min && min < 1) {
				list1.push_back(LineSeg(line1.S, minP));
				if (max <= 1.0) {
					list1.push_back(line2);
					if (max < 1.0)
						list1.push_back(LineSeg(maxP, line1.E));
				}
				else list1.push_back(LineSeg(minP, line1.E));
			}
			else if (0 < max && max <= 1) {
				list1.push_back(LineSeg(line1.S, maxP));
				if (max < 1.0)
					list1.push_back(LineSeg(maxP, line1.E));
			}
			if (min < 0 && 0 < max) {
				list2.push_back(LineSeg(minP, line1.S));
				if (1 <= max) {
					list2.push_back(line1);
					if (1 < max)
						list2.push_back(LineSeg(line1.E, maxP));
				}
				else list2.push_back(LineSeg(line1.S, maxP));
			}
			else if (min < 1 && 1 <= max) {
				list2.push_back(LineSeg(minP, line1.E));
				if (1 < max)
					list2.push_back(LineSeg(line1.E, maxP));
			}
			if (list1.size() == 0) list1.push_back(line1);
			if (list2.size() == 0) list2.push_back(line2);
			return make_tuple(list1, list2);
		}

		// неколлинеарны
		// пересечение по двум вершинам
		if (line1.S == line2.S ||
			line1.S == line2.E ||
			line1.E == line2.S ||
			line1.E == line2.E)
		{
			list1.push_back(line1);
			list2.push_back(line2);
			return make_tuple(list1, list2);
		}

		assert(seg1_line2_start != seg1_line2_end && "Collinearity mismatch2");
		assert(seg2_line1_start != seg2_line1_end && "Collinearity mismatch3");
		auto u1 = seg1_line2_start / (seg1_line2_start - seg1_line2_end);
		auto u2 = seg2_line1_start / (seg2_line1_start - seg2_line1_end);
		if (u1 == .0 || abs(u1) == 1.0)
		{
			// пересечение по стороне line2 и вершине line1
			intersection = u1 < .5 ? &line1.S : &line1.E;
			list1.push_back(line1);
			list2.push_back(LineSeg(line2.S, *intersection));
			list2.push_back(LineSeg(line2.E, *intersection));
		}
		else if (u2 == .0 || abs(u2) == 1.0)
		{   // пересечение по стороне line1 и вершине line2
			intersection = u2 < .5 ? &line2.S : &line2.E;
			list2.push_back(line2);
			list1.push_back(LineSeg(line1.S, *intersection));
			list1.push_back(LineSeg(line1.E, *intersection));
		}
		else
		{   // пересечение посередине
			intersection = new Vector2D(line1.S + dir1 * u1);
			list1.push_back(LineSeg(line1.S, *intersection));
			list1.push_back(LineSeg(line1.E, *intersection));
			list2.push_back(LineSeg(line2.S, *intersection));
			list2.push_back(LineSeg(line2.E, *intersection));
		}
		return make_tuple(list1, list2);
	}
	std::vector<FragmentGeometry> GeometryUtils::FindContours(const std::vector<LineSeg>& nonIntersectedLines)
	{
		auto capacity = nonIntersectedLines.size() * 2;
		auto result = std::vector<FragmentGeometry>();									result.reserve(capacity);
		auto g = std::vector<std::vector<std::tuple<double, int, double>>>();			     g.reserve(capacity);
		auto points = std::vector<Vector2D>();											points.reserve(capacity);
		auto t = &Vector2D::GetHash;
		std::unordered_map<Vector2D, int, decltype(t) > map(0, t); map.reserve(capacity);
		auto minX = DBL_MAX;
		auto mostLeftPointIndex = -1, counter = 0;
		auto addPoint = [&](Vector2D p) -> void
		{
			if (map.count(p) > 0) return;
			map.emplace(p, counter);
			g.push_back(std::vector< std::tuple < double, int, double>>());
			points.push_back(p);
			if (p.X < minX)
			{
				mostLeftPointIndex = counter;
				minX = p.X;
			}
			counter++;
		};
		for (auto line : nonIntersectedLines)
		{
			addPoint(line.S);
			addPoint(line.E);
		}
		auto n = points.size();
		if (counter < 3) return result;
		{   // начинаем с самой крайней точки, чтобы первый контур получился внешним
			// todo: сортировка по координате X
			map[points[mostLeftPointIndex]] = 0;
			map[points[0]] = mostLeftPointIndex;
			std::swap(points[0], points[mostLeftPointIndex]);
		}

		for (auto line : nonIntersectedLines)
		{
			auto i1 = map[line.S];
			auto i2 = map[line.E];
			auto v = line.Vector();
			auto value1 = atan2(v.Y, v.X);
			auto value2 = atan2(-v.Y, -v.X); //value1 + (value1 <= .0 ? Math.PI : -Math.PI);
			auto length = line.Length();
			g[i1].push_back(std::make_tuple(value1, i2, length));
			g[i2].push_back(std::make_tuple(value2, i1, length));
		}

		// для каждой вершины сортируем исходящие ребра по углу
		auto g2 = std::vector<std::vector<int>>(n);
		for (size_t i = 0; i < n; i++) {
			auto& list = g[i];
			std::sort(list.begin(), list.end(), [](std::tuple < double, int, double> x, std::tuple < double, int, double> y)
			{
				if (std::get<0>(x) < std::get<0>(y)) return true;
				if (std::get<0>(x) > std::get<0>(y)) return  false;
				if (std::get<2>(x) < std::get<2>(y)) return true;
				if (std::get<2>(x) > std::get<2>(y)) return  false;
				return true;
			});
			g2[i].resize(list.size());
			std::transform(list.begin(), list.end(), g2[i].begin(), [](std::tuple < double, int, double> t) { return std::get<1>(t); });
		}

		std::vector<bool> used(n*n, false);
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < g2[i].size(); ++j)
				if (!used[i*n + j])
				{
					used[i*n + j] = true;
					int v = g2[i][j], pv = i;
					auto facet = FragmentGeometry();
					for (;;)
					{
						facet.ContourPoints.push_back(points[v]);

						// ищем индекс предыдущей вершины
						size_t it;
						for (it = 0; it < g2[v].size(); it++)
							if (g2[v][it] == pv)
								break;
						assert(it < g2[v].size());
						it = (it == g2[v].size() - 1) ? 0 : it + 1;

						if (used[v*n + it]) break;
						used[v*n + it] = true;
						pv = v;
						v = g2[v][it];
					}
					result.push_back(facet);
				}
		//
		auto newCounter = 0;
		for (size_t j = 0; j < result.size(); j++)
		{
			auto& contourPoints = result[j].ContourPoints;
			for (auto oldCount = -1; oldCount != contourPoints.size() && contourPoints.size() > 2;)
			{
				oldCount = contourPoints.size();
				for (size_t i = 0; i < contourPoints.size() && contourPoints.size() > 2; i++)
				{
					auto prevInd = Utils::PrevInd(i, contourPoints.size());
					auto nextInd = Utils::NextInd(i, contourPoints.size());
					if (contourPoints[prevInd] == contourPoints[nextInd])
					{
						contourPoints.erase(contourPoints.begin() + std::max(i, nextInd)); // не эффективно
						contourPoints.erase(contourPoints.begin() + std::min(i, nextInd)); // не эффективно
						i--;
					}
				}
			}
			if (contourPoints.size() > 2 && ContourSignedArea(contourPoints) < 0)
				result[newCounter++] = result[j];
		}
		result.resize(newCounter);
		return result;
	}
	bool GeometryUtils::IsPointInsideContour(const Vector2D& p, const std::vector<Vector2D>& contour, bool isBoundaryInside)
	{
		assert(contour.size() > 2);
		auto parity = false;
		auto start = contour.back();
		for (auto& end : contour)
		{
			auto a = end - start;
			auto b = p - start;
			auto sa = Vector2D::CrossProduct(a, b);
			if ((sa > 0 && start.Y < p.Y && p.Y <= end.Y) ||
				(sa < 0 && end.Y < p.Y && p.Y <= start.Y))
				parity = !parity;
			else if (sa == 0 && a.X * b.X >= 0 && a.Y * b.Y >= 0 && a.LengthSquared() >= b.LengthSquared())
				return isBoundaryInside; // on edge
			start = end;
		}
		return parity;
	}
	double GeometryUtils::ContourSignedArea(const std::vector<Vector2D>& points)
	{
		auto area = 0.0;
		auto li = points.size() - 1;
		for (size_t i = 0; i < points.size(); li = i++)
			area += (points[i].X + points[li].X) * (points[i].Y - points[li].Y);
		return area / 2;
	}

	void NonIntersectedCurvesPool::AddLine(LineSeg _line)
	{
		auto queue = std::queue<std::tuple<LineSeg, unsigned int>>();
		queue.push(std::make_tuple(_line, 0));
		while (queue.size() > 0)
		{
			auto linet = queue.front(); queue.pop();
			auto line = std::get<0>(linet);
			if (line.S == line.E) continue;
			for (auto j = std::get<1>(linet); j < Lines.size(); j++)
			{
				auto r = GeometryUtils::DivideLines(Lines[j], line);
				auto buf1 = std::get<0>(r), buf2 = std::get<1>(r);
				assert(buf1.size() > 0 && buf2.size() > 0);
				if (buf1.size() > 1)
				{
					Lines[j] = buf1[0];
					for (size_t i = 1; i < buf1.size(); i++)
						Lines.push_back(buf1[i]);
				}
				if (buf2.size() == 1) continue;
				for (auto subline : buf2)
					queue.push(std::make_tuple(subline, j));
				goto End;
			}
			// когда появляются дубликаты?
			auto dublicate = false;
			for (auto l : Lines)
				if (dublicate = (l == line))
					break;
			if (!dublicate) Lines.push_back(line);
		End:
			;
		}
	}
	bool GeometryUtils::Circumcenter(const Vector2D& p1, const Vector2D& p2, const Vector2D& p3, Vector2D& center)
	{	// Shewchuk's circumcenter location
		// todo: is there is a way to calculate circumcenter and radius simultaneously?
		auto v1 = p2 - p1;
		auto v2 = p3 - p1;
		auto l1 = v1.LengthSquared();
		auto l2 = v2.LengthSquared();
		auto d = 4 * TriangleOrientedArea(p1, p2, p3);
		if (d == 0) return false;
		center = p1 + (Vector2D(v2.Y*l1 - v1.Y*l2, v1.X*l2 - v2.X*l1) / d);
		return true;
	}

	inline bool IsLargeArc(const Vector2D& p1, const Vector2D& p2, const Vector2D& c, bool clockwise)
	{
		auto crossProduct = Vector2D::CrossProduct(p1 - c, p2 - c) > 0;
		return (clockwise && crossProduct) || (!clockwise && !crossProduct);
	}
	bool GeometryUtils::Circumcenter(const Vector2D& p1, const Vector2D& p2, double radius, bool clockwise, bool largeArc, Vector2D& center)
	{
		auto v = p2 - p1;
		auto d = v.Length();
		auto h = radius*radius - d*d / 4;
		if (h <= 0 || d == 0) return false;
		h = sqrt(h);
		Vector2D n(v.Y, -v.X);
		auto c1 = 0.5 * (p1 + p2) + h / d * n;
		auto c2 = 0.5 * (p1 + p2) - h / d * n;
		center = (largeArc == IsLargeArc(p1, p2, c1, clockwise) ? c1 : c2);
		return true;
	}

	//------------------------------------------------------------
	//					robust predicates
	//------------------------------------------------------------
	// todo: exactinit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	extern "C" double counterclockwise(const double* pa, const double* pb, const double* pc);
	extern "C" double incircle(const double* pa, const double* pb, const double* pc, const double* pd);
	extern "C" void exactinit();

	static bool robust_predicates_initialized = false;

	bool GeometryUtils::InCircle(const Vector2D& a, const Vector2D& b, const Vector2D& c, const Vector2D& d, bool clockwise)
	{
		assert(Orientation(a, b, c) == (clockwise ? PointsOrder::Clockwise : PointsOrder::Counterclockwise));
		if (!robust_predicates_initialized) {
			exactinit();
			robust_predicates_initialized = true;
		}
		auto res = incircle(a.ToArray(), b.ToArray(), c.ToArray(), d.ToArray());
		return (res > 0 && !clockwise) || (res < 0 && clockwise);
	}

	double GeometryUtils::TriangleOrientedArea(const Vector2D& a, const Vector2D& b, const Vector2D& c)
	{
		return Vector2D::CrossProduct(b - a, c - a) / 2;
		if (!robust_predicates_initialized) {
			exactinit();
			robust_predicates_initialized = true;
		}
		return counterclockwise(a.ToArray(), b.ToArray(), c.ToArray()) / 2;
	}

	double GeometryUtils::TriangleAlphaQuality(const Vector2D& a, const Vector2D& b, const Vector2D& c)
	{
		return 4 * sqrt(3)*TriangleOrientedArea(a, b, c) / ((a - b).LengthSquared() + (b - c).LengthSquared() + (a - c).LengthSquared());
	}

	PointsOrder GeometryUtils::Orientation(const Vector2D& a, const Vector2D& b, const Vector2D& c)
	{
		return Orientation(TriangleOrientedArea(a, b, c));
	}

	PointsOrder GeometryUtils::Orientation(double triangleOrientedArea)
	{
		// todo: unexact version
		const double eps = 1e-8;                                    // !!! epsilon !!!
		if (abs(triangleOrientedArea) < eps) return PointsOrder::Collinear;
		return triangleOrientedArea < 0 ? PointsOrder::Clockwise : PointsOrder::Counterclockwise;
	}

	// очень специфичная функция
	// отрезки должны быть неколлинеарны
	// концы первого могут лежать на втором, но не наоборот
	bool GeometryUtils::IntersectionOfNonCollinearLines(const Vector2D& start1, const Vector2D& end1, const Vector2D& start2, const Vector2D& end2, Vector2D& intersection)
	{
		if (start1 == start2) { intersection = start1; return true; }
		if (start1 == end2) { intersection = start1; return true; }
		if (end1 == start2) { intersection = end1; return true; }
		if (end1 == end2) { intersection = end1; return true; }

		//подставляем концы отрезков, для выяснения в каких полуплоскотях они
		auto seg1_line2_start = TriangleOrientedArea(start2, end2, start1);
		auto seg1_line2_end = TriangleOrientedArea(start2, end2, end1);
		auto seg2_line1_start = TriangleOrientedArea(start1, end1, start2);
		auto seg2_line1_end = TriangleOrientedArea(start1, end1, end2);

		//если концы одного отрезка имеют один знак, значит он в одной полуплоскости и пересечения нет.
		if (seg1_line2_start * seg1_line2_end > 0 ||
			seg2_line1_start * seg2_line1_end > 0)
			return false;

		if (seg1_line2_start == 0 && seg1_line2_end == 0) return false; // а это для чего? todo: заменить на assert 

		auto u1 = seg1_line2_start / (seg1_line2_start - seg1_line2_end);
		auto i = start1 + u1 * (end1 - start1);
		if (u1 < 0.5)
		{
			if (PointsOrder::Collinear == Orientation(start2, start1, i) ||
				PointsOrder::Collinear == Orientation(end2, start1, i))
			{
				intersection = start1; return true;
			}
		}
		else if (PointsOrder::Collinear == Orientation(start2, end1, i) ||
			PointsOrder::Collinear == Orientation(end2, end1, i))
		{
			intersection = end1; return true;
		}
		intersection = i;
		return true;
	}

}