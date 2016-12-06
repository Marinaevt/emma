#include "includes.h"
#include "Utils.h"
#include "GeometryUtils.h"
#include "GeometryPrimitives.h"
#include "Mesher.h"
#include <iterator>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;


namespace FemLibrary {

	bool MovePointMSW(const Vector2D& old_point, const vector<pair<Triangle*, int>>& vertex_triangles, Vector2D& new_point);
	bool MovePointNewton(const Vector2D& old_point, const vector<pair<Triangle*, int>>& vertex_triangles, Vector2D& new_point);

	class HalfEdge
	{
	public:
		static int creation_counter; // todo: remove after unit testd were completed

		string ToString() const { return Start.ToString() + " -> " + Sym->Start.ToString(); }
		Vector2D Start;
		HalfEdge* Sym = nullptr;
		HalfEdge* Prev = nullptr;
		HalfEdge* Next = nullptr;
		const Vector2D& getEnd() const { return Sym->Start; }
		void setEnd(const Vector2D& value) const { Sym->Start = value; }

		explicit HalfEdge(const Vector2D& point) : Start(point) { }

		static pair<HalfEdge*, HalfEdge*> CreatePair(const Vector2D& start, const Vector2D& end)
		{
			auto a = new HalfEdge(start);
			auto b = new HalfEdge(end);

			assert(creation_counter += 2);

			a->Sym = b;
			b->Sym = a;
			return make_pair(a, b);
		}

		void ConnectTo(HalfEdge* he)
		{
			assert(Next == nullptr && Sym->Prev == nullptr && he->Prev == nullptr && he->Sym->Next == nullptr);
			Next = he;
			Sym->Prev = he->Sym;
			he->Prev = this;
			he->Sym->Next = Sym;
		}
		pair<HalfEdge*, HalfEdge*> AddNext(const Vector2D& next)
		{
			auto res = CreatePair(getEnd(), next);
			ConnectTo(res.first);
			return res;
		}
	private: // use delete edge instead
		void Delete() const           // Next\   / Sym.Prev
		{                             //      \e/
			Next->Prev = Sym->Prev;   //       |       |
			Prev->Next = Sym->Next;   //       ^   Sym |
			Sym->Next->Prev = Prev;   //       |       V
			Sym->Prev->Next = Next;   //      /s\
									  								  // Prev/   \ Sym.Next
		}
	public:
		static void DeleteEdge(HalfEdge* he)
		{
			he->Delete();
			delete he->Sym;
			delete he;
			assert(creation_counter -= 2);
		}
		void ConnectToEdgeFromLeft(HalfEdge* e)
		{                           //             /|\|
			Next = e;               //              | 
			Sym->Prev = e->Prev;    //   Next       | e
			e->Prev->Next = Sym;    // -->-->-->-->-(
			e->Prev = this;         //   SymPrev   / \
																	//      e.Prev/   \ |
		}

		bool operator==(const HalfEdge& rhs) const { return Start == rhs.Start && getEnd() == rhs.getEnd(); }
		bool operator!=(const HalfEdge& rhs) const { return !this->operator==(rhs); }
		static size_t GetHash(const HalfEdge& p) { return Utils::hash_combine(Vector2D::GetHash(p.Start), Vector2D::GetHash(p.getEnd())); }
	};
	int HalfEdge::creation_counter = 0;
}

namespace std {
	template <> struct hash<FemLibrary::HalfEdge> {
		size_t operator()(const FemLibrary::HalfEdge& p) const { return FemLibrary::HalfEdge::GetHash(p); }
	};
}

namespace FemLibrary{

class TriangleSearchTable
{
	BoundingBox _bounds;
	vector<Triangle*> _table;
	vector<Triangle*> _deleteList;
	unordered_set<Triangle*> _aliveTriangles;
    int KoefResizingTriTable = 4;
	int _tableSize = 1;
public:
	void Clear()
	{
		for (auto triangle : _aliveTriangles)
			delete triangle;
		for (auto triangle : _deleteList)
			delete triangle;
		_aliveTriangles.clear();
		_deleteList.clear(); _deleteList.shrink_to_fit();
		_table.resize(_tableSize = 1); _table.shrink_to_fit();
	}
	bool IsAlive(Triangle* t) const { assert(t != nullptr); return _aliveTriangles.find(t) != _aliveTriangles.end(); }
	void Set(const Vector2D& p, Triangle* t)
	{
		assert(t != nullptr);
		_table[PointToTriTableKey(p)] = t;
	}
	Triangle* Get(const Vector2D& p) const
	{
		auto res = _table[PointToTriTableKey(p)];
		assert(res != nullptr);
		return res;
	}
	const unordered_set<Triangle*>& GetSet() const { assert(CheckBoundedness()); return _aliveTriangles; }
	int GetNumberOfAliveTriangles() const { return static_cast<int>(_aliveTriangles.size()); }
	TriangleSearchTable& Add(Triangle* triangle)
    { 
		assert(triangle != nullptr);
		assert(_bounds.Contains(triangle->V[0]) && _bounds.Contains(triangle->V[1]) && _bounds.Contains(triangle->V[2]));
		auto res = _aliveTriangles.emplace(triangle);
		assert(res.second);
		_table[PointToTriTableKey(triangle->Centroid)] = triangle;
		while (GetNumberOfAliveTriangles() > _tableSize * _tableSize * KoefResizingTriTable)
			Resize();
		return *this;
    }
	TriangleSearchTable& Replace(Triangle* oldTri, Triangle* newTri)
    {
		ReplaceWithExistent(oldTri, newTri);
		Add(newTri);
		return *this;
    }
	TriangleSearchTable& ReplaceWithExistent(Triangle* oldTri, Triangle* existentTri)
	{
		assert(oldTri != nullptr && existentTri != nullptr);
		auto res = _aliveTriangles.erase(oldTri);
		assert(res == 1);

		for (auto i = 0; i < _table.size(); ++i)
			if (_table[i] == oldTri)
				_table[i] = existentTri;
		// мы не применяем физическое удаление, т.к. адреса треугольников используются в качестве ключей хеш-таблиц и внутри очереди с приоритетом
		// и при не удачном расположении звезд новый треугольник возникнет на месте старого
		//delete oldTri;
		_deleteList.emplace_back(oldTri);
		return *this;
	}
	//TriangleSearchTable& ReplaceIfExists(Triangle* oldTri, Triangle* newTri)
    //{
	//	assert(oldTri != nullptr && newTri != nullptr);
	//	auto res = _aliveTriangles.erase(oldTri);
	//	if(res == 1)
	//	{
	//		for (auto i = 0; i < _table.size(); ++i)
	//			if (_table[i] == oldTri)
	//				_table[i] = newTri;
	//		delete oldTri;
	//	}
	//	Add(newTri);
	//	return *this;
    //}

    void BuildFromTriangles(const vector<Triangle*>& triangles, const BoundingBox& bounds)
    {
		Clear();
        
		if (triangles.size() <= 0) throw exception("empty triangles vector");

		_bounds = bounds;
		_table[0] = triangles[0];

		for(auto tri: triangles)
			Add(tri);
    }

	//=============================================================================================
	//									Searching
	//=============================================================================================
public:
	enum class PointPositionInTriangle { InsideTriangle, TriangleVertex, TriangleSide, OutsideTriangle };

	tuple<const Triangle*, PointPositionInTriangle, int> GetPointInTriangleInfo(const Vector2D& p)
	{
		if (!_bounds.Contains(p))
			return make_tuple(static_cast<Triangle*>(nullptr), PointPositionInTriangle::OutsideTriangle, -1);
			//throw exception(("GetPointInTriangleInfo: point " + p.ToString() + " is outside triangulation bounds").c_str());

		auto tri = Get(p);
		auto initialTri = tri;
#ifdef _DEBUG
		unordered_set<Triangle*> path = { tri };
#endif
		PointsOrder orders[3];
		pair<PointPositionInTriangle, int> pos;
		while (true)
		{
			pos = GetPointPositionInTriangle(p, tri, orders);
			if (pos.first != PointPositionInTriangle::OutsideTriangle) break;
			// если точка не находится в текущем треугольнике, то ищем в смежном треугольнике, 
			// противоположная вершина для которого лежит по другую сторону смежного ребра относительно нашей точки
			Triangle* newTri = nullptr;
			for (auto i = 0; i < 3 && newTri == nullptr; i++)
				if (tri->Neighbour[i] != nullptr && tri->Orientation != orders[i] && orders[i] != PointsOrder::Collinear)
					newTri = tri->Neighbour[i];
			if (newTri == nullptr)
				return make_tuple(tri, PointPositionInTriangle::OutsideTriangle, -1);
			tri = newTri;
#ifdef _DEBUG
			if (!path.emplace(tri).second)
				throw exception(("GetPointInTriangleInfo: path recursion, point:" + p.ToString()).c_str());
#endif
		}
		if (initialTri != tri) Set(p, tri);
		return make_tuple(tri, pos.first, pos.second);
	}

private:
	static pair<PointPositionInTriangle, int> GetPointPositionInTriangle(const Vector2D& p, const Triangle* t, PointsOrder(&orders)[3])
	{
		auto ind = t->VertexInd(p);
		if (ind >= 0) return make_pair(PointPositionInTriangle::TriangleVertex, ind);

		// todo: use InCircle test
		//if (Math.Abs((t.Center - p).LengthSquared - t.BoundingCircleSquareRadius) > Eps)
		//    return Tuple<PointPositionInTriangle, int>(PointPositionInTriangle.OutsideTriangle, -1);
		orders[0] = GeometryUtils::Orientation(t->E[0].S, t->E[0].E, p);
		orders[1] = GeometryUtils::Orientation(t->E[1].S, t->E[1].E, p);
		orders[2] = GeometryUtils::Orientation(t->E[2].S, t->E[2].E, p);
		if (orders[0] == orders[1] && orders[1] == orders[2])
		{
			assert(orders[0] != PointsOrder::Collinear && orders[0] == t->Orientation);
			return make_pair(PointPositionInTriangle::InsideTriangle, -1);
		}
		if (orders[0] == PointsOrder::Collinear && orders[1] == orders[2]) return make_pair(PointPositionInTriangle::TriangleSide, 0);
		if (orders[1] == PointsOrder::Collinear && orders[2] == orders[0]) return make_pair(PointPositionInTriangle::TriangleSide, 1);
		if (orders[2] == PointsOrder::Collinear && orders[0] == orders[1]) return make_pair(PointPositionInTriangle::TriangleSide, 2);
		return make_pair(PointPositionInTriangle::OutsideTriangle, -1);
	}

private:
	bool CheckBoundedness() const
	{
		unordered_set<Triangle*> triangleSet;
		for (auto triangleQueue = _table; triangleQueue.size() > 0; )
		{
			auto tri = triangleQueue.back();
			triangleQueue.pop_back();

			if (!triangleSet.emplace(tri).second) continue;

			for (auto nTri : tri->Neighbour)
				if (nTri != nullptr)
					triangleQueue.push_back(nTri);
		}
		for (auto tri : _aliveTriangles) { // Symmetric except with _aliveTriangles
			auto res = triangleSet.emplace(tri);
			if (!res.second) triangleSet.erase(res.first);
		}
		return triangleSet.size() == 0;
	}
    int PointToTriTableKey(const Vector2D& p) const
    {
        if(!_bounds.Contains(p))
            throw exception(("point {"+p.ToString()+"} not in triangle table at all: "+_bounds.ToString()).c_str());
		auto o = p - _bounds.Min();
        auto x = min(_tableSize - 1, static_cast<int>(floor(o.X * _tableSize / _bounds.DeltaX())));
        auto y = min(_tableSize - 1, static_cast<int>(floor(o.Y * _tableSize / _bounds.DeltaY())));
        return x + y * _tableSize;
    }
    void Resize()
    {
		auto tbuffer = _table;
		auto newSize = _tableSize * 2;
		_table.resize(newSize * newSize);
		for (auto i = 0; i < _tableSize;++i)
			for (auto j = 0; j < _tableSize;++j)
			{
				auto i2 = i * 2, j2 = j * 2;
				_table[ i2    *newSize + j2    ] =
				_table[ i2    *newSize + j2 + 1] =
				_table[(i2+1) *newSize + j2    ] =
				_table[(i2+1) *newSize + j2 + 1] = tbuffer[i*_tableSize + j];
			}
        _tableSize = newSize;
#ifdef _DEBUG
		for (auto t : _table)
			assert(t != nullptr && t->Orientation != PointsOrder::Collinear);
#endif
    }
};

class Mesh::MeshImpl
{
public:
	MeshingStatistics stat;
	MeshingParameters params;
	FragmentGeometry cg;
	MeshingResult result;
	vector<IntermediateStep> intermediateSteps;

	BoundingBox bounds;
	TriangleSearchTable triTable;

	// result
	vector<SimpleTriangle> triangles;

	bool Generate()
	{
		try
		{
			result.Clear();
			triangles.clear();
			intermediateSteps.clear();
			processStatisticsStart();

			// 1 этап: деление структурных ребер в зависимости от условия
			unordered_set<Vector2D> structurePoints;
			unordered_set<LineSeg> structureSegments;
			process_input2(structurePoints, structureSegments); // todo: добавить поддержку дырок
			processStatisticsAfterInputProcessing();

			// 2 этап: триангуляция Делоне
			auto delaunayTriangles = DivideAndConquer(structurePoints);
			triTable.BuildFromTriangles(delaunayTriangles, bounds);
			if (params.saveIntermediateSteps) SaveIntermediateStep();
			processStatisticsAfterConvexDelaunay();

			assert(CheckNeighboursMatch(triTable.GetSet()));
			assert(CheckValidDelaunay(triTable.GetSet()));

			if(params.ConvexDelaynayOnly)
				return process_result();

			// 3 этап: ограниченная триангуляция Делоне
			unordered_set<LineSeg> structureSegmentsToInsert(structureSegments);
			for(auto tri: triTable.GetSet())
			{
				for(auto i=0;i<3;++i)
				{
					auto it = structureSegmentsToInsert.find(tri->E[i]);
					if (it != structureSegmentsToInsert.end()) {
						tri->SetStructureEdge(i);
						structureSegmentsToInsert.erase(it);
					}
				}
			}
			for (auto& structureSegment : structureSegmentsToInsert)
			{
				auto newTriangles = InsertStructureEdgeAnglada(structureSegment.S, structureSegment.E);
				FlipAllInconditionalTriangles(newTriangles);
			}
			assert(CheckStructureEdges(triTable.GetSet(), structureSegments));
			if (params.saveIntermediateSteps) SaveIntermediateStep();

			// помечаем внешние треугольники
			for (auto tri : triTable.GetSet())
				tri->IsExternal = !GeometryUtils::IsPointInsideContour(tri->Centroid, cg.ContourPoints); // todo: обработка дырок
			if (params.saveIntermediateSteps) SaveIntermediateStep();
			processStatisticsAfterConstrainedDelaunay();

			assert(CheckNeighboursMatch(triTable.GetSet()));
			assert(CheckValidDelaunay(triTable.GetSet()));

			// 
			//process_boundary(structurePoints, structureSegments);
			//if (saveIntermediateSteps) SaveIntermediateStep();
			//assert(CheckNeighboursMatch(triTable.GetSet()));
			//assert(CheckValidDelaunay(triTable.GetSet()));

			// 4 этап: распространение фронта
			unordered_set<Triangle*> acceptedTriangles;
			auto compareFunc = [](const pair<double, Triangle*>& x, const pair<double, Triangle*>& y) -> bool { return x.first > y.first; };
			priority_queue<pair<double, Triangle*>, vector<pair<double, Triangle*>>, decltype(compareFunc)> activeTrianglesQueue(compareFunc);
			ProcessNewTriangles(triTable.GetSet(), acceptedTriangles, activeTrianglesQueue);
			while (activeTrianglesQueue.size() > 0)
			{
				// проверка превышения количества треугольников
				if (params.MaxNumberOfTriangles > 0 && triTable.GetNumberOfAliveTriangles() > params.MaxNumberOfTriangles) {
					auto numberOfInternalTriangles = 0;
					for (auto t : triTable.GetSet())
						if (!t->IsExternal)
							numberOfInternalTriangles++;
					if (numberOfInternalTriangles > params.MaxNumberOfTriangles)
						throw exception(("Exceeded number of triangles " + to_string(params.MaxNumberOfTriangles)).c_str());
				}

				auto pair = activeTrianglesQueue.top(); activeTrianglesQueue.pop();
				auto tri = pair.second;
				if (!triTable.IsAlive(tri)) continue;

				// find shortest edge of this triangle located on the front
				auto ind = -1;
				double minLength = DBL_MAX;
				for (auto i = 0; i < 3; ++i)
					if (tri->IsStructureEdge(i) || acceptedTriangles.find(tri->Neighbour[i]) != acceptedTriangles.end()) {
						auto len = tri->E[i].Length();
						if (len < minLength) {
							ind = i;
							minLength = len;
						}
					}
				assert(ind >= 0);

				auto insertionPoint = GetInsertionPoint(tri, ind);
				vector<Triangle*> newTriangles;
				if (!InsertPoint(insertionPoint, newTriangles))
				{
					// надо исследовать, к чему это приводит
					//throw exception(("Can't insert point " + insertionPoint.ToString()).c_str());
					acceptedTriangles.emplace(tri);
					continue;
				}

				assert(CheckNeighboursMatch(triTable.GetSet()));
				assert(CheckValidDelaunay(triTable.GetSet()));

				ProcessNewTriangles(newTriangles, acceptedTriangles, activeTrianglesQueue);

				if (params.saveIntermediateSteps) SaveIntermediateStep(&insertionPoint);
				//_intermediateStates.Add(Tuple.Create(insertionPoint, acceptedTriangles.Where(t => _triTable.AliveTriangles.Contains(t)).ToList()));
			}
			processStatisticsAfterAFM();

			// 5 этап: сглаживание
			if (params.smoothing_params.smoothing)
				Smooth(structurePoints, structurePoints);
				//Smooth(params.smoothing_params.smoothingIterationsMultiplier, structurePoints);
			processStatisticsEnd();

			// 6 этап: обработка ответа
#ifdef _DEBUG
			// check that all alive triangles are acccepted or external
			for (auto triangle : triTable.GetSet()) 
				if (acceptedTriangles.end() == acceptedTriangles.find(triangle)) 
					throw exception("Not accepted internal triangle found!");
#endif
			return process_result();
		}
		catch(exception &exc)
		{
			result.errorMsg = string(exc.what());
			bool bRes = process_result();
			return result.success = false;
			return result.success = true; // for debugging
		}
	}
private:
	void process_input2(unordered_set<Vector2D>& structurePoints, unordered_set<LineSeg>& structureEdges)
	{
		if (cg.ContourPoints.size() < 3)
			throw exception("Not enought points to preform triangulation");

		bounds = BoundingBox(cg.ContourPoints);
		if (bounds.DeltaX() * bounds.DeltaY() <= 0.0)
			throw exception("Invalid points to preform triangulation");

		structurePoints.clear();
		structureEdges.clear();

		// todo: process holes
		auto prevPoint = cg.ContourPoints.back();
		for (auto const& point : cg.ContourPoints) {
			auto division = divide_edge(LineSeg(prevPoint, point));

			for (auto& p: division) structurePoints.emplace(p);
			for (size_t i = 1; i < division.size();++i) structureEdges.emplace(LineSeg(division[i-1], division[i]));
			prevPoint = point;
		}
	}
	double rib_length(const Vector2D& pos) const
	{
		return params.RibLength;
	}
	vector<Vector2D> divide_edge(const LineSeg& edge) const
	{
		const auto number_of_integration_steps = 100;
		const auto number_of_interations = 10;
		auto v = edge.Vector();

		// 
		auto T = 0.0;
		auto h = v / number_of_integration_steps; //Шаг сетки
		auto a = edge.S + h / 2;
		for (auto i = 0; i < number_of_integration_steps; ++i)
			T += 1.0 / rib_length(a + h * i); 
		auto I = ceil(T * h.Length());

		auto delta = v / I;
		vector<Vector2D> s; 
		s.reserve(I + 2);
		s.emplace_back(edge.S);
		
		for(auto i=1;i<I;++i)
			s.emplace_back(edge.S + i*delta);
		s.emplace_back(edge.E);

		//
		for (auto iteration = 0; iteration < number_of_interations; ++iteration)
		{
			int ind;
			for (ind = 1; ind < I; ++ind)
			{
				auto dmi  = rib_length((s[ind    ] + s[ind - 1]) / 2);
				auto dmi1 = rib_length((s[ind + 1] + s[ind    ]) / 2);

				s[ind] = ((s[ind + 1] + s[ind - 1]) + (s[ind + 1] - s[ind - 1]) * (dmi - dmi1) / (dmi + dmi1)) / 2;
			}
			for (ind = I - 1; ind > 0; --ind)
			{
				auto dmi  = rib_length((s[ind    ] + s[ind - 1]) / 2);
				auto dmi1 = rib_length((s[ind + 1] + s[ind    ]) / 2);

				s[ind] = ((s[ind + 1] + s[ind - 1]) + (s[ind + 1] - s[ind - 1]) * (dmi - dmi1) / (dmi + dmi1)) / 2;
			}
		}
		return s;
	}
	

	bool process_result()
	{
		auto allTriangles = triTable.GetSet();
		result.Clear();
		triangles             .clear();
		triangles		      .reserve(allTriangles.size());
		result.points         .reserve(allTriangles.size());
		result.triangleIndexes.reserve(allTriangles.size() * 3);
		unordered_map<Vector2D, int> points_indexes;
		for (auto triangle : allTriangles) {
			if (triangle->IsExternal) continue;

			triangles.emplace_back(triangle->ToSimpleTriangle());
			for(auto& point: triangle->V)
			{
				auto res = points_indexes.emplace(point, static_cast<int>(result.points.size()));
				if(res.second)
					result.points.emplace_back(point);
				auto ind = (*res.first).second;
				result.triangleIndexes.emplace_back(ind);
			}
		}
		triTable.Clear();
		return result.success = true;
	}
	//---------------------------------------------
	//				Common
	//---------------------------------------------
private:
	static void UpdateReference(Triangle* triOld, Triangle* triNew, int indOld, int indNew)
	{
		auto neighbour = triOld->Neighbour   [indOld];
		auto nind      = triOld->NeighbourInd[indOld];
		if (neighbour != nullptr) // для случая граничных треугольников
		{
			assert(triNew->E[indNew] == neighbour->E[nind]);
			neighbour->Neighbour      [nind] = triNew;
			neighbour->NeighbourInd   [nind] = indNew;
			neighbour->isStructureEdge[nind] = triOld->isStructureEdge[indOld];
		}
		triNew->Neighbour      [indNew] = neighbour;
		triNew->NeighbourInd   [indNew] = nind;
		triNew->isStructureEdge[indNew] = triOld->isStructureEdge[indOld];
		assert(CheckNeighbourMatch(triNew));
		assert(neighbour == nullptr || CheckNeighbourMatch(neighbour));
	}
	pair<Triangle*, Triangle*> FlipTriangles(Triangle* tri, int ind)
	{
		//Assert(CheckNeighbourPointMatch(tri), "FlipTriangles: CheckNeighbourMatch 0");
		auto nTri = tri->Neighbour[ind];
		assert(!tri->IsExternal && !nTri->IsExternal);
		auto nind = tri->NeighbourInd[ind];
		auto&  Points =  tri->V;
		auto& nPoints = nTri->V;
		auto IndS  = Utils::NextInd( ind, 3);
		auto IndE  = Utils::PrevInd( ind, 3);
		auto nIndS = Utils::NextInd(nind, 3);
		auto nIndE = Utils::PrevInd(nind, 3);
		if (nPoints[nIndS] != Points[IndS]) swap(nIndS, nIndE);
		if (nPoints[nIndS] != Points[IndS] || nPoints[nIndE] != Points[IndE])
			throw exception("DelaunayTriangulation.FlipTriangles: can't match points");

		// todo: replace creation with modification
		auto t1 = new Triangle(Points[ind], Points[IndS], nPoints[nind]);//       S
		auto t2 = new Triangle(Points[ind], Points[IndE], nPoints[nind]);//     / t1 \
		// общим вершинам 0 и 2 соответствует ребро 1                    // ind --|--- nind
		t1->Neighbour   [1] = t2;                                        //     \ t2 /
		t2->Neighbour   [1] = t1;                                        //       E
		t1->NeighbourInd[1] = t2->NeighbourInd[1] = 1;

		UpdateReference(tri, t1, IndE, 2);
		UpdateReference(tri, t2, IndS, 2);
		UpdateReference(nTri, t1, nIndE, 0);
		UpdateReference(nTri, t2, nIndS, 0);
		triTable.Replace(tri, t1).Replace(nTri, t2);
		assert(CheckNeighbourMatch(t1));
		assert(CheckNeighbourMatch(t2));
		return make_pair(t1, t2);
	}
	void FlipAllInconditionalTriangles(const vector<Triangle*>& triList)
	{
		auto counter = 0;
		vector<Triangle*> triQueue(triList); // 
		//vector<Triangle*> newTriangles;
		while (triQueue.size() > 0)
		{
			auto tri = triQueue.back(); triQueue.pop_back();
			if (!triTable.IsAlive(tri))
				continue;
			for (auto i = 0; i < 3; ++i) // ищем индекс первого треугольника-соседа не прошедшего проверку Делоне
				if (!DelaunayCheck(tri, i))
				{
					auto pair = FlipTriangles(tri, i);
					triQueue.emplace_back(pair.first ); // newTriangles.emplace_back(pair.first );
					triQueue.emplace_back(pair.second); // newTriangles.emplace_back(pair.second);
					break;
				}
			if (counter++ > 100000) // todo: remove this check after testing
				throw exception("FlipAllInconditionalTriangles: endless recursion");
		}
		//return newTriangles;
	}
	//=============================================
	//          Intermediate states
	//=============================================
private:
	void SaveIntermediateStep(const Vector2D* point = nullptr) // todo: save only internal triangles
	{
		auto triangleSet = triTable.GetSet();
		intermediateSteps.emplace_back(IntermediateStep());
		auto& step = intermediateSteps.back();

		//step.triangles.resize(triangleSet.size());	
		//transform(triangleSet.begin(), triangleSet.end(), step.triangles.begin(), [](Triangle* t) { return t->ToSimpleTriangle(); });
		for (auto t : triangleSet)
			if (!t->IsExternal)
				step.triangles.emplace_back(t->ToSimpleTriangle());
		
		step.newPoint = point != nullptr ? *point : Vector2D(0, 0);
	}
	//=============================================
	//          advancing front refinement
	//=============================================
private:
	double DesiredCircumscribedRadius(const Vector2D& point) const
	{
		return rib_length(point) / sqrt(3);
	}
	bool IsAcceptedTriangle(const Triangle* tri) const
	{
		if (tri->IsExternal || (tri->IsObtuse() && tri->IsStructureEdge(tri->ObtuseVertexInd)))
			return true;
		return tri->CircleRadius < DesiredCircumscribedRadius(tri->Centroid) * params.ToleranceRatio;
		//foreach(var edge in tri.E)
		//	if (edge.Length > maxRibLength)
		//		return false;
		//return true;
	}
	Vector2D GetInsertionPoint(const Triangle* tri, int ind) const
	{
		// ищем точку
		//![image alt text](img\\new_pointh.jpg "title")

		auto middle = tri->E[ind].Middle();
		auto normal = tri->Circumcenter - middle;
		auto normalLength = normal.Length();
		auto p = tri->E[ind].Length() / 2;
		auto radius = min(max(DesiredCircumscribedRadius(middle), p), 0.5*(p*p + normalLength*normalLength) / normalLength);
		auto offset = radius + sqrt(radius*radius - p*p);
		auto insertionPoint = middle + offset / normalLength*normal;
		//var offset = Math.Sqrt(3) / 2 * maxRibLength;
		//var insertionPoint = normalLength < offset ? tri.Circumcenter : middle + offset / normalLength * normal;
		return insertionPoint;
	}
	template<class T, class Q>
	void ProcessNewTriangles(T& newTriangles, unordered_set<Triangle*>& acceptedTriangles, Q& activeTrianglesQueue)
	{
		// find all accepted triangles
		for(Triangle* tri: newTriangles)
			if (IsAcceptedTriangle(tri)) {
				auto res = acceptedTriangles.emplace(tri);
				assert(res.second);
			}
		// find all active triangles
		for(Triangle* tri: newTriangles)
			if (!IsAcceptedTriangle(tri))
			{
				for (auto i = 0; i < 3; ++i)
					if (tri->IsStructureEdge(i) || acceptedTriangles.find(tri->Neighbour[i]) != acceptedTriangles.end())
					{
						activeTrianglesQueue.push(make_pair(tri->CircleRadius, tri));//tri.E[tri.ShortestEdgeInd].Length
						break;
					}
			}
	}
	bool InsertPoint(const Vector2D& insertionPoint, vector<Triangle*>& newTriangles)
	{	// весьма хреновая реализация — куча не нужных действий и структур
		auto pointInTriangleInfo = triTable.GetPointInTriangleInfo(insertionPoint);
		auto pointPosition = get<1>(pointInTriangleInfo);
		if (pointPosition == TriangleSearchTable::PointPositionInTriangle::TriangleVertex ) return true;
		if (pointPosition == TriangleSearchTable::PointPositionInTriangle::OutsideTriangle) return false;

		unordered_set<Triangle*> visitedTriangles;
		vector<pair<Triangle*, int>> borderEdges;
		{
			vector<pair<Triangle*, int>> triangleQueue;
			triangleQueue.emplace_back(make_pair(const_cast<Triangle*>(get<0>(pointInTriangleInfo)), -1));  // const cast
			while (triangleQueue.size() > 0)
			{
				auto pair = triangleQueue.back(); triangleQueue.pop_back();
				auto curTri = pair.first;
				if (!visitedTriangles.emplace(curTri).second) continue;

				for (auto i = 0; i < 3; ++i)
				{
					if (i == pair.second) continue;
					auto neighbour = curTri->Neighbour[i];
					if (curTri->IsStructureEdge(i) || !InCircle(neighbour, insertionPoint))
						//(insertionPoint - neighbour->Circumcenter).Length() > neighbour->CircleRadius) // !GeometryUtils.InCircle(neighbour, insertionPoint))
						borderEdges.emplace_back(make_pair(curTri, i));
					else triangleQueue.emplace_back(make_pair(neighbour, curTri->NeighbourInd[i]));
				}
			}
		}

		assert(borderEdges.size() > 0);
#ifdef _DEBUG
		for(auto t: borderEdges)
			assert(visitedTriangles.find(t.first) != visitedTriangles.end());
#endif
		vector<pair<Triangle*, Triangle*>> replaceList;
		for(auto borderEdgePair: borderEdges)
		{
			auto borderTriangle = borderEdgePair.first ;
			auto borderInd      = borderEdgePair.second;
			auto newTriangle = new Triangle(insertionPoint, borderTriangle->NextV(borderInd),
															borderTriangle->PrevV(borderInd));
			newTriangles.emplace_back(newTriangle);

			if (visitedTriangles.erase(borderTriangle) == 1) replaceList.emplace_back(make_pair(borderTriangle, newTriangle));
			else											 triTable.Add(newTriangle);
		}
		FillAllNeighbours(newTriangles);
		for (size_t i = 0; i < newTriangles.size(); ++i)
			UpdateReference(borderEdges[i].first, newTriangles[i], borderEdges[i].second, 0);

		for (auto& p : replaceList) triTable.Replace(p.first, p.second);
		for(auto visitedTriangle: visitedTriangles)
			triTable.ReplaceWithExistent(visitedTriangle, newTriangles[0]);
		return true;
	}
	void Smooth(double iterationMultiplier, const unordered_set<Vector2D>& structurePointSet) const
	{
		vector<bool> isStructurePoint;
		vector<Vector2D> vertexes;
		unordered_map<Vector2D,int> vertexIndexes;
		vector<vector<pair<Triangle*,int>>> vertexTriangles;
		for(auto triangle: triTable.GetSet())
		{
			if (triangle->IsExternal) continue;
			for (auto index = 0; index < 3; ++index)
			{
				auto& vertex = triangle->V[index];
				auto p = vertexIndexes.emplace(vertex, static_cast<int>(vertexes.size()));
				if (p.second) // если добавление прошло успешно
				{
					vertexes        .emplace_back(vertex);
					vertexTriangles .emplace_back(vector<pair<Triangle*, int>>());
					isStructurePoint.emplace_back(structurePointSet.find(vertex) != structurePointSet.end());
				}
				vertexTriangles[(*p.first).second].emplace_back(make_pair(triangle, index));
			}
		}
		auto nVertexes = vertexes.size();


		auto compareFunc = [](const pair<double, int>& x, const pair<double, int>& y) -> bool { return x.first > y.first; };
		priority_queue<pair<double, int>, vector<pair<double, int>>, decltype(compareFunc)> pointQueue(compareFunc);
		vector<vector<int>> vertexNeighbours    (nVertexes);
		vector<double>      vertexActualQuality (nVertexes, 0.);
		auto recalcQuality = [&vertexTriangles](int i) -> double {
			auto& star = vertexTriangles[i];

			auto meanQuality = accumulate(star.begin(), star.end(), 1.0, [](double acc, const pair<Triangle*, int>& p)
			{
				return acc*p.first->AlphaQuality;
			});
			return pow(abs(meanQuality), 1.0 / star.size());
		};
		for (auto i = 0; i < nVertexes; ++i)
		{
			if (isStructurePoint[i]) continue; // пропускаем структурные вершины
				
			pointQueue.push(make_pair(vertexActualQuality[i] = recalcQuality(i), i));

			// заполняем списки соседних вершин
			for(auto& pair: vertexTriangles[i])
			{
				auto starVertex = pair.first->Orientation == PointsOrder::Clockwise ? pair.first->NextV(pair.second) : pair.first->PrevV(pair.second);
				vertexNeighbours[i].emplace_back(vertexIndexes.at(starVertex));
			}
		}

		auto smoothCounter = 0;
		while (pointQueue.size() > 0)
		{
			auto qualityIndPair = pointQueue.top(); pointQueue.pop();
			auto ind = qualityIndPair.second;
			auto oldMeanQuality = qualityIndPair.first;
			if (oldMeanQuality != vertexActualQuality[ind])
				continue;
			assert(!isStructurePoint[ind]);

			auto& star       = vertexTriangles [ind];
			auto& neighbours = vertexNeighbours[ind];

			// находим центроид
			auto newPoint = accumulate(neighbours.begin(), neighbours.end(), Vector2D(0, 0), [&vertexes](const Vector2D& acc, int i)
			{
				return acc + vertexes[i];
			}) / static_cast<double>(neighbours.size());

			auto isValidMove = true;
			auto calcQuality = [&newPoint](const Triangle* t, int i, const Vector2D& point) -> double {
				auto& p = t->PrevV(i);
				auto& n = t->NextV(i);
				auto signedArea = GeometryUtils::TriangleOrientedArea(p, point, n);
				if (GeometryUtils::Orientation(signedArea) != t->Orientation) return 0.0;
				return 4 * sqrt(3) * signedArea / (t->E[i].Vector().LengthSquared() + (p - newPoint).LengthSquared() + (n - newPoint).LengthSquared());
			};
			auto meanQuality = accumulate(star.begin(), star.end(), 1.0, [&newPoint, &calcQuality](double acc, const pair<Triangle*, int>& p)
			{
				return acc * calcQuality(p.first, p.second, newPoint);
			});
			if (isValidMove && oldMeanQuality < pow(abs(meanQuality), 1.0 / star.size()))
			{
				vertexes[ind] = newPoint;
				for(auto& pair: star)
					pair.first->MoveVertex(pair.second, newPoint);
				for(auto neighbourInd: neighbours)
					if(!isStructurePoint[neighbourInd])
						pointQueue.push(make_pair(vertexActualQuality[neighbourInd] = recalcQuality(neighbourInd), neighbourInd));
			}
			smoothCounter++;
			if (smoothCounter > nVertexes*iterationMultiplier)
				return;
		}
	}
	struct SmoothingVertex
	{
		bool isStructurePoint = false, isBorderPoint = false;
		Vector2D vertex;
		vector<pair<Triangle*, int>> star;
		unordered_set<Vector2D> neighbourVertexSet;
		vector<SmoothingVertex*> neighbourVertexes;
		double meanQ = 1.0, minQ = 1.0;
		explicit SmoothingVertex(const Vector2D& vertex):vertex(vertex){}
		void AddStarTriangle(Triangle* tri, int ind)
		{
			star.emplace_back(tri, ind);
			neighbourVertexSet.emplace(tri->NextV(ind));
			neighbourVertexSet.emplace(tri->PrevV(ind));
		}
		void ProcessNeighbours(const unordered_map<Vector2D, SmoothingVertex*>& vertexDict)
		{
			// заполняем списки соседних вершин
			for (auto& vertex : neighbourVertexSet)
				neighbourVertexes.emplace_back(vertexDict.at(vertex));
		}
		void RecalcQuality()
		{
			// рассчитываем среднее и минимальное качество
			meanQ = 1.0;
			minQ = DBL_MAX;
			for (auto& pair : star)
			{
				meanQ *= pair.first->AlphaQuality;
				minQ = min(minQ, pair.first->AlphaQuality);
			}
			meanQ = pow(abs(meanQ), 1.0 / star.size());
		}
		bool TryToMoveVertex(const Vector2D& newPos, bool checkMeanQuality)
		{
			auto oldMinQ = minQ;
			auto oldMeanQ = meanQ;
			auto invertedTriangles = false;
			for (auto& starPair : star)
			{
				auto t = starPair.first;
				auto i = starPair.second;
				t->MoveVertex(i, newPos);
				if (t->AlphaQuality <= 0)
					invertedTriangles = true;
			}
			RecalcQuality();
			if(invertedTriangles || minQ < oldMinQ || (checkMeanQuality && meanQ < oldMeanQ))
			{
				for (auto& starPair : star)
					starPair.first->MoveVertex(starPair.second, vertex);
				minQ = oldMinQ;
				meanQ = oldMeanQ;
				return false;
			}
			vertex = newPos;
			for (auto neighbour : neighbourVertexes)
				neighbour->RecalcQuality();
			return true;
		}

	};
	bool Smooth(const unordered_set<Vector2D>& structurePointSet, const unordered_set<Vector2D>& borderPointSet)
	{
		if (params.smoothing_params.smoothingIterationsMultiplier < 0 || 
			params.smoothing_params.minQualitySmoothingRatio < 0 || 
			params.smoothing_params.minQualitySmoothingRatio > 1 || 
			params.smoothing_params.moveThreshold <= 0)
			return false;

		auto triangleSet = triTable.GetSet();

		// меняем ориентацию всех треугольников на PointsOrder::Counterclockwise, меняя местами вершины 0 и 1
		for (auto triangle : triangleSet)
		{
			assert(triangle->Orientation != PointsOrder::Collinear);
			if (triangle->Orientation == PointsOrder::Counterclockwise) continue;
			// else PointsOrder::Clockwise
			triangle->Neighbour[0]->NeighbourInd[triangle->NeighbourInd[0]] = 1;
			triangle->Neighbour[1]->NeighbourInd[triangle->NeighbourInd[1]] = 0;
			swap(triangle->Neighbour   [0], triangle->Neighbour   [1]);
			swap(triangle->NeighbourInd[0], triangle->NeighbourInd[1]);
			swap(triangle->E[0], triangle->E[1]);
			swap(triangle->V[0], triangle->V[1]);
			swap(triangle->E[0].S, triangle->E[0].E);
			swap(triangle->E[1].S, triangle->E[1].E);
			swap(triangle->E[2].S, triangle->E[2].E);
			swap(triangle->isStructureEdge[0], triangle->isStructureEdge[1]);
			triangle->Orientation = PointsOrder::Counterclockwise;
			     if (triangle->ObtuseVertexInd == 1) triangle->ObtuseVertexInd = 0;
			else if (triangle->ObtuseVertexInd == 0) triangle->ObtuseVertexInd = 1;
		}
		assert(CheckNeighboursMatch(triangleSet));

		unordered_map<Vector2D, SmoothingVertex*> vertexDict;
		for (auto triangle : triangleSet)
		{
			if (triangle->IsExternal) continue;
			for (auto index = 0; index < 3; ++index)
			{
				auto& vertex = triangle->V[index];
				auto p = vertexDict.emplace(vertex, nullptr);
				if (p.second) // если добавление прошло успешно
				{
					auto smoothingVertex = new SmoothingVertex(vertex);
					p.first->second = smoothingVertex;
					smoothingVertex->isStructurePoint = structurePointSet.find(vertex) != structurePointSet.end();
					smoothingVertex->isBorderPoint = false; // todo
				}
				p.first->second->AddStarTriangle(triangle, index);
			}
		}
		for (auto vertexPair : vertexDict) {
			vertexPair.second->ProcessNeighbours(vertexDict);
			vertexPair.second->RecalcQuality();
		}

		auto iteration = 0; // счетчик количества операций сглаживания (успешных или всех?)
		auto totalSmoothingOperations = static_cast<int>(vertexDict.size() * params.smoothing_params.smoothingIterationsMultiplier);
		stat.minQualityIterations = static_cast<int>(totalSmoothingOperations * params.smoothing_params.minQualitySmoothingRatio);
		stat.meanQualityIterations = totalSmoothingOperations - stat.minQualityIterations;
		assert(stat.minQualityIterations >= 0 && stat.meanQualityIterations >= 0);

		auto compareFunc  = [](const pair<double, SmoothingVertex*>& x, const pair<double, SmoothingVertex*>& y) -> bool { return x.first > y.first; };

		// todo: решить, как поступать в случае невозможности выполнить шаг оптимизации
		//=========================================================
		//					min quality smoothing
		//=========================================================
		{
			priority_queue<pair<double, SmoothingVertex*>, vector<pair<double, SmoothingVertex*>>, decltype(compareFunc)> pointQueue(compareFunc);

			for (auto& vertexPair : vertexDict)
				if (!vertexPair.second->isStructurePoint)
					pointQueue.emplace(vertexPair.second->minQ, vertexPair.second);

			stat.minQualitySuccessfulMoveCounter = stat.minQualityTooSmalllMoveCounter = 0;
			for (; pointQueue.size() > 0 && iteration < stat.minQualityIterations; ++iteration)
			{
				auto qualityIndPair = pointQueue.top(); pointQueue.pop();
				auto smoothingVertex = qualityIndPair.second;
				auto oldQuality = qualityIndPair.first;		// если качество вершины не совпадает с ключом очереди приоритетов,
				if (oldQuality != smoothingVertex->minQ) {	// значит эта вершина или ее соседка уже перемещалась, значит трогать ее больше нельзя
					iteration--;							// счетчик итераций при этом изменяться не должен, ибо ничего не считалось
					continue;								
				}
				assert(!smoothingVertex->isStructurePoint);

				Vector2D newPoint;
				if (!MovePointMSW(smoothingVertex->vertex, smoothingVertex->star, newPoint))
					continue;

				auto minimalMove = rib_length(smoothingVertex->vertex) * params.smoothing_params.moveThreshold;
				if (minimalMove * minimalMove > (newPoint - smoothingVertex->vertex).LengthSquared()) {
					stat.minQualityTooSmalllMoveCounter++;
					continue;
				}

				if (smoothingVertex->TryToMoveVertex(newPoint, false))
				{
					stat.minQualitySuccessfulMoveCounter++;

					pointQueue.emplace(smoothingVertex->minQ, smoothingVertex);
					for (auto neighbour : smoothingVertex->neighbourVertexes)
						if (!neighbour->isStructurePoint)
							pointQueue.emplace(neighbour->minQ, neighbour);

					if (stat.isExtended && iteration*stat.displaySize >= stat.minQualityDuringOptimization.size()*totalSmoothingOperations)
						processStatisticsDuringOptimization(iteration);
				}
			}
			stat.minQualityMoveFailureCounter = stat.minQualityIterations - stat.minQualitySuccessfulMoveCounter - stat.minQualityTooSmalllMoveCounter;
		}
		processStatisticsDuringOptimization(iteration);
		processStatisticsAfterSmoothingMin();
		//=========================================================
		//					mean quality smoothing
		//=========================================================
		{
			priority_queue<pair<double, SmoothingVertex*>, vector<pair<double, SmoothingVertex*>>, decltype(compareFunc)> pointQueue(compareFunc);

			for (auto& vertexPair : vertexDict)
				if (!vertexPair.second->isStructurePoint)
					pointQueue.emplace(vertexPair.second->meanQ, vertexPair.second);

			stat.meanQualitySuccessfulMoveCounter = stat.meanQualityTooSmalllMoveCounter = 0;
			for (; pointQueue.size() > 0 && iteration < totalSmoothingOperations; ++iteration)
			{
				auto qualityIndPair = pointQueue.top(); pointQueue.pop();
				auto smoothingVertex = qualityIndPair.second;
				auto oldQuality = qualityIndPair.first;		// если качество вершины не совпадает с ключом очереди приоритетов,
				if (oldQuality != smoothingVertex->meanQ) {	// значит эта вершина или ее соседка уже перемещалась, значит трогать ее больше нельзя
					iteration--;							// счетчик итераций при этом изменяться не должен, ибо ничего не считалось
					continue;
				}
				assert(!smoothingVertex->isStructurePoint);

				Vector2D newPoint;
				if (!MovePointNewton(smoothingVertex->vertex, smoothingVertex->star, newPoint))
					continue;

				auto minimalMove = rib_length(smoothingVertex->vertex) * params.smoothing_params.moveThreshold;
				if (minimalMove * minimalMove > (newPoint - smoothingVertex->vertex).LengthSquared()) {
					stat.meanQualityTooSmalllMoveCounter++;
					continue;
				}

				if (smoothingVertex->TryToMoveVertex(newPoint, true))
				{
					stat.meanQualitySuccessfulMoveCounter++;

					pointQueue.emplace(smoothingVertex->meanQ, smoothingVertex);
					for (auto neighbour : smoothingVertex->neighbourVertexes)
						if (!neighbour->isStructurePoint)
							pointQueue.emplace(neighbour->meanQ, neighbour);

					if (stat.isExtended && iteration*stat.displaySize >= stat.minQualityDuringOptimization.size()*totalSmoothingOperations)
						processStatisticsDuringOptimization(iteration);
				}
			}
			processStatisticsDuringOptimization(iteration);
			stat.meanQualityMoveFailureCounter = stat.meanQualityIterations - stat.meanQualitySuccessfulMoveCounter - stat.meanQualityTooSmalllMoveCounter;
		}
		processStatisticsAfterSmoothingMean();

		// освобождаем память VertexDict
		for (auto vertexPair : vertexDict)
			delete vertexPair.second;
		return true;
	}
	//=============================================
    //          Divide and conquer
    //=============================================
	// todo: надо что-нибудь придумать с ориентацией, например заставить всегда использовать неточную версию
	// или при делении структурных ребер добавить функцию выравнивания по линии,
	// ибо при обработке вертикальных линий на границе контура возможно попадание большого количества коллинеарных точек в один "букет"
	// в таком случае алгоритм может "сойти с ума"
private:
	static void FillAllNeighbours(const vector<Triangle*>& triangles)
	{
		// заполняем всех соседей в треугольниках
		unordered_map<LineSeg, pair<Triangle*, int>> mapEdge2Triangles;
		for (auto triangle : triangles)
			for (auto i = 0; i < 3; ++i)
			{
				auto res = mapEdge2Triangles.emplace(triangle->E[i], make_pair(triangle, i));
				if (!res.second)
				{
					auto pair = (*res.first).second;
					auto nTri = pair.first;
					auto nInd = pair.second;
					assert(nTri->Neighbour[nInd] == nullptr && triangle->Neighbour[i] == nullptr);

					triangle->Neighbour[i] = nTri;
					triangle->NeighbourInd[i] = nInd;
					nTri->Neighbour[nInd] = triangle;
					nTri->NeighbourInd[nInd] = i;
				}
			}
	}
    static pair<HalfEdge*, HalfEdge*> DivideAndConquer(const vector<Vector2D>& points, int start, int n)
    {
		pair<HalfEdge*, HalfEdge*> res;
        if (n > 3)
        {
            auto i = (n / 2);
            auto L = DivideAndConquer(points, start, i);
            auto R = DivideAndConquer(points, start + i, n - i);

            HalfEdge *ldi = L.second, *rdi = R.first ;
            HalfEdge *ldo = L.first , *rdo = R.second;
            while (true)
            {
                if (PointsOrder::Counterclockwise == GeometryUtils::Orientation(ldi->Start, ldi->getEnd(), rdi->Start))
                    ldi = ldi->Next;
                else if (PointsOrder::Clockwise == GeometryUtils::Orientation(rdi->Start, rdi->getEnd(), ldi->Start))
                    rdi = rdi->Sym->Prev->Sym;                                      // ldo                 rdo
                else break;                                                        //   \   \ /        \ /
            }                                                                      //    ————*==<==<==<=*
            auto basePair = HalfEdge::CreatePair(rdi->Start, ldi->Start);              //     <-ldi        rdi->
            auto basel = basePair.first;
            basel->ConnectToEdgeFromLeft(ldi);
            basel->Sym->ConnectToEdgeFromLeft(rdi->Sym->Next);

            if (ldo->Start == ldi->Start) ldo = basel->Sym;
            if (rdo->Start == rdi->Start) rdo = basel;

            // merge loop
            while (true)
            {
                auto rvalid = false, lvalid = false;
                auto lcand = basel->Sym->Prev->Sym; // lcand <- basel.Sym.Onext; не так!!!
                //IF Valid[lcand] THEN           Valid[e]=RightOf[e.Dest,basel] = CCW[e.Dest, basel.Dest, basel.Org].
                if (PointsOrder::Counterclockwise == GeometryUtils::Orientation(lcand->getEnd(), basel->getEnd(), basel->Start))
                {
                    //    WHILE GeometryUtils::InCircle[basel.Dest, basel.Org, lcand.Dest, lcand.Onext.Dest] DO
                    while (GeometryUtils::InCircle(basel->getEnd(), basel->Start, lcand->getEnd(), lcand->Next->getEnd()))
                    {
                        auto t = lcand->Prev->Sym;  //  t <- lcand.Onext;
                        HalfEdge::DeleteEdge(lcand);//   DeleteEdge[lcand];
                        lcand = t;                  //    lcand <- t;
                    } // todo: optimize
                    lvalid = PointsOrder::Counterclockwise == GeometryUtils::Orientation(lcand->getEnd(), basel->getEnd(), basel->Start);
                }

                auto rcand = basel->Sym->Next; // rcand <- basel.Oprev;
                // IF Valid[rcand] THEN        Valid[e]=RightOf[e.Dest,basel] = CCW[e.Dest, basel.Dest, basel.Org].
                if (PointsOrder::Counterclockwise == GeometryUtils::Orientation(rcand->getEnd(), basel->getEnd(), basel->Start))
                {
                    //     WHILE GeometryUtils::InCircle[basel.Dest, basel.Org, rcand.Dest, rcand.Oprev.Dest] DO
                    while (GeometryUtils::InCircle(basel->getEnd(), basel->Start, rcand->getEnd(), rcand->Sym->Next->getEnd()))
                    {
                        auto t = rcand->Sym->Next;  // t <- rcand.Oprev;
						HalfEdge::DeleteEdge(rcand);//  DeleteEdge[rcand];
                        rcand = t;                  //   rcand <- t;
                    }
                    rvalid = PointsOrder::Counterclockwise == GeometryUtils::Orientation(rcand->getEnd(), basel->getEnd(), basel->Start);
                }

                if (!rvalid && !lvalid) break;

                // IF NOT Valid[lcand] OR(Valid[rcand] AND GeometryUtils::InCircle[lcand.Dest, lcand.Org, rcand.Org, rcand.Dest]) THEN
                if (!lvalid || (rvalid && GeometryUtils::InCircle(lcand->getEnd(), lcand->Start, rcand->Start, rcand->getEnd())))
                {
                    // Add cross edge base1 from rcand.Dest to basel.Dest:)
                    auto pair = HalfEdge::CreatePair(rcand->getEnd(), basel->getEnd()); // basel < -Connect[rcand, basel.Sym];
                    auto nbasel = pair.first;                          //          /* rcand.End
                    nbasel->ConnectToEdgeFromLeft(basel->Sym);         //  \     /  |  /
                    nbasel->Sym->ConnectToEdgeFromLeft(rcand->Next);   //   \  /    | /
                    basel = nbasel;                                    // --- E=<=<=S
                }
                else
                { // Add cross edge base1 from basel.Org to lcand.Dest
                    auto pair = HalfEdge::CreatePair(basel->Start, lcand->getEnd()); // base1 <- Connect[basel.Sym, lcand.Sym]
                    auto nbasel = pair.first;                               //     *\lcand.End
                    nbasel->ConnectToEdgeFromLeft(lcand->Sym);				//  \  |  \    /
                    nbasel->Sym->ConnectToEdgeFromLeft(basel->Sym->Next);   //   \ |    \ /
                    basel = nbasel;                                         //     E=<=<=S----
                }
            }
            res = make_pair(ldo, rdo);
        }
        else if (n == 3)
        {
            auto pair1 = HalfEdge::CreatePair(points[start], points[start + 1]);
            auto pair2 = pair1.first->AddNext(points[start + 2]);
            pair<HalfEdge*, HalfEdge*> pair3;

            auto order = GeometryUtils::Orientation(points[start], points[start + 1], points[start + 2]);
            if (order != PointsOrder::Collinear)
            {
                pair3 = pair2.first->AddNext(points[start]);
                pair3.first->ConnectTo(pair1.first);
            }
            else
            {
                pair2.first->Next = pair2.second;
                pair2.second->Prev = pair2.first;
                pair1.second->Next = pair1.first;
                pair1.first->Prev = pair1.second;
            }
            res = order == PointsOrder::Clockwise
                ? make_pair(pair3.second, pair3.first)
                : make_pair(pair1.first, pair2.second);
        }
        else if (n == 2)
        {
            res = HalfEdge::CreatePair(points[start], points[start + 1]);
            res.first ->Next = res.first ->Prev = res.second;
            res.second->Next = res.second->Prev = res.first ;
        }
        else assert(false);
        return res;
    }
	static vector<Triangle*> DivideAndConquer(const unordered_set<Vector2D>& pointSet)
    {
		HalfEdge::creation_counter = 0;

		vector<Vector2D> pointsUnique(pointSet.begin(), pointSet.end());
		sort(pointsUnique.begin(), pointsUnique.end(),[](const Vector2D& a, const Vector2D& b)->bool
		{
			if (a.X < b.X) return true;
			if (a.X > b.X) return false;
			if (a.Y < b.Y) return true;
			if (a.Y > b.Y) return false;
			assert(false); // сортировка так реализована, что иногда одно значение сравнивает с самим собой
			return false;
		});
        auto he = DivideAndConquer(pointsUnique, 0, static_cast<int>(pointsUnique.size())).second;
		assert(he != nullptr);

        // извлекаем все полуребра
        unordered_set<HalfEdge*> allEdges;
        {

            vector<HalfEdge*> queue;
            for (queue.push_back(he); queue.size() > 0;)
            {
                auto edge = queue.back();
				queue.pop_back();

				if (!allEdges.emplace(edge).second) continue;

                if (edge->Next != nullptr) queue.push_back(edge->Next);
				if (edge->Prev != nullptr) queue.push_back(edge->Prev);
				if (edge->Sym  != nullptr) queue.push_back(edge->Sym );
            }
        }

		vector<Vector2D> contour;;
		unordered_set<const HalfEdge*> processedEdges;
		auto findContour = [&contour, &processedEdges](const HalfEdge* edge) {
			contour.clear();
			for (auto cur = edge; processedEdges.emplace(cur).second; cur = cur->Next)
				contour.push_back(cur->Start);
		};

		// ищем внешний контур
        findContour(he);

        // ищем все треугольники
        vector<Triangle*> triangles;
        for (auto edge: allEdges)
        {
			findContour(edge);
			if (contour.size() <= 0) continue;
            if (contour.size() != 3)
                throw exception("invalid triangle found");
			assert(PointsOrder::Collinear != GeometryUtils::Orientation(contour[0], contour[1], contour[2]));
            triangles.push_back(new Triangle(contour[0], contour[1], contour[2]));
        }

		// удаляем все полуребра
		assert(HalfEdge::creation_counter == allEdges.size());
		for (auto halfEdge : allEdges)
			delete halfEdge;

        FillAllNeighbours(triangles);
		return triangles;
    }
	//--------------------------------------
	//				Constrained
	//--------------------------------------
private:
	//              /\			   	 рекурсивная перетриангуляция полупустоты, полученной удалением треугольников, пересекающих структурное ребро
	//         /\  /  \			     принцип: 1. надуваем окружность на ребре ab до первой встреченной вершины
	//        /  \/    \    /\				  2. строим полученный треугольник и рекурсивно выполняем перетриангуляцию для двух других сторон
	//     /\/          \  /  \					
	//    /             с\/    \
	// a *----------------------* b
	tuple<Triangle*, int, bool> TriangulatePseudopolygonDelaunay(
		const vector<Vector2D>& vertexList,
		int sInd, int size,
		const Vector2D& a, const Vector2D& b,
		vector<Triangle*>& newTriangles,
		const unordered_map<LineSeg, tuple<Triangle*, int, bool>>& externalNeighbours,
		bool orientation) const
	{
		if (size == 0) return externalNeighbours.at(LineSeg(a, b)); // throws exception if fails
		
		// находим индекс первой встреченной вершины
		auto cInd = sInd;
		for (auto vInd = sInd+1; vInd < sInd+size; ++vInd)
		{
			auto r = orientation ? GeometryUtils::InCircle(a, b, vertexList[cInd], vertexList[vInd]) // todo: replace with one InCircle call
				                 : GeometryUtils::InCircle(b, a, vertexList[cInd], vertexList[vInd]);
			if (r) cInd = vInd;
		}

		// выполняем рекурсивно триангуляцию для левой и правой сторон
		auto c = vertexList[cInd];
		auto left  = TriangulatePseudopolygonDelaunay(vertexList,     sInd,            cInd - sInd, a, c, newTriangles, externalNeighbours, orientation);
		auto right = TriangulatePseudopolygonDelaunay(vertexList, cInd + 1, size - 1 - cInd + sInd, c, b, newTriangles, externalNeighbours, orientation);

		auto tri = new Triangle(a, c, b);
		auto rTri = get<0>(right);
		auto lTri = get<0>( left);
		auto rInd = get<1>(right);
		auto lInd = get<1>(left);
		tri->Neighbour[0] = rTri;
		tri->Neighbour[2] = lTri;
		tri->NeighbourInd[0] = rInd;
		tri->NeighbourInd[2] = lInd;
		tri->isStructureEdge[0] = get<2>(right);
		tri->isStructureEdge[2] = get<2>( left);
		if (lTri != nullptr)
		{
			lTri->Neighbour   [lInd] = tri;
			lTri->NeighbourInd[lInd] = 2;
			tri->isStructureEdge[2] = lTri->isStructureEdge[lInd];
		}
		if (rTri != nullptr)
		{
			rTri->Neighbour   [rInd] = tri;
			rTri->NeighbourInd[rInd] = 0;
			tri->isStructureEdge[0] = rTri->isStructureEdge[rInd];
		}
		newTriangles.emplace_back(tri);
		return make_tuple(tri, 1, false);
	}
	Vector2D FindCavityHull(
		Triangle* initialTriangle,		                          // треугольник, который содержит первый участок структурного ребра
		int indEdge,					                          // индекс вершины и ребра, через которые проходит структурное ребро
		const Vector2D& startPoint, const Vector2D& endPoint,     // концы структурного ребра
		vector<Vector2D>& upperPart,                              // список вершин над структурным ребром
		vector<Vector2D>& lowerPart,						      // список вершин под структурным ребром
		unordered_map<LineSeg, tuple<Triangle*, int, bool>>& neighbours, // треугольники, внешние по отношению к найденной пустоте, требуются для перетриангуляции
		vector<Triangle*>& trianglesToDelete) const				  // треугольники внутри пустоты
	{
		auto curTriangle = initialTriangle;
		upperPart.clear(); lowerPart.clear(); neighbours.clear(); trianglesToDelete.clear();
		trianglesToDelete.emplace_back(curTriangle);

		// строим списки верхних и нижних вершин
		auto nextInd = Utils::NextInd(indEdge, 3);
		auto prevInd = Utils::PrevInd(indEdge, 3);
		auto oS = GeometryUtils::Orientation(startPoint, endPoint, curTriangle->V[prevInd]);
		auto oE = GeometryUtils::Orientation(startPoint, endPoint, curTriangle->V[nextInd]);
		assert(oS != PointsOrder::Collinear && oE != PointsOrder::Collinear);
		(oS == PointsOrder::Clockwise ? lowerPart : upperPart).emplace_back(curTriangle->V[prevInd]);
		(oE == PointsOrder::Clockwise ? lowerPart : upperPart).emplace_back(curTriangle->V[nextInd]);
		assert(upperPart.size() == 1 && lowerPart.size() == 1);
		neighbours.emplace(curTriangle->E[nextInd], make_tuple(curTriangle->Neighbour[nextInd], curTriangle->NeighbourInd[nextInd], curTriangle->isStructureEdge[nextInd]));
		neighbours.emplace(curTriangle->E[prevInd], make_tuple(curTriangle->Neighbour[prevInd], curTriangle->NeighbourInd[prevInd], curTriangle->isStructureEdge[prevInd]));
		PointsOrder o;
		Vector2D curP;
		do
		{
			//     uInd *
			//         / \ 
			//   -----/---\---->
			//       /     \
			// lInd *-------* nInd
			auto nind = curTriangle->NeighbourInd[indEdge];
			trianglesToDelete.emplace_back(curTriangle = curTriangle->Neighbour[indEdge]);
			assert(curTriangle != nullptr);
			auto lInd = Utils::PrevInd(nind, 3);
			auto uInd = Utils::NextInd(nind, 3);
			   if (upperPart.back() == curTriangle->V[lInd]) swap(lInd, uInd);
			assert(lowerPart.back() == curTriangle->V[lInd]);

			curP = curTriangle->V[nind];
			o = GeometryUtils::Orientation(startPoint, endPoint, curP);
			if (o != PointsOrder::Clockwise)
			{
				if (o != PointsOrder::Collinear) upperPart.emplace_back(curP);
				indEdge = uInd;
				neighbours.emplace(curTriangle->E[lInd], make_tuple(curTriangle->Neighbour[lInd], curTriangle->NeighbourInd[lInd], curTriangle->isStructureEdge[lInd]));
			}
			if (o != PointsOrder::Counterclockwise)
			{
				if (o != PointsOrder::Collinear) lowerPart.emplace_back(curP);
				indEdge = lInd;
				neighbours.emplace(curTriangle->E[uInd], make_tuple(curTriangle->Neighbour[uInd], curTriangle->NeighbourInd[uInd], curTriangle->isStructureEdge[uInd]));
			}

		} while (o != PointsOrder::Collinear);
		return curP;
	}
	vector<Triangle*> InsertStructureEdgeAnglada(const Vector2D& startPoint, const Vector2D& endPoint)
	{
		vector<Triangle*> totalNewTriangles;

		auto curPoint = startPoint;
		auto info = triTable.GetPointInTriangleInfo(curPoint);
		assert(get<1>(info) == TriangleSearchTable::PointPositionInTriangle::TriangleVertex);
		auto curTriangle = const_cast<Triangle*>(get<0>(info)); // kludge: const_cast
		while (curPoint != endPoint)
		{
			// ищем по всем соседям трегольники, содержащие данную точку
			auto ind2 = curTriangle->VertexInd(curPoint);
			assert(ind2 >= 0);
			auto curTriList = curTriangle->GetStar(ind2);

			// находим список всех возможных пересечений
			// todo: заменить на нахождение первого пересечения, сейчас я нахожу список всех пересечений в целях отладки
			vector<tuple<const Triangle*, Vector2D, int, bool>> resList;
			for(auto& pair: curTriList)
			{ // todo: заменить IntersectionOfNonCollinearLines на более пристойный вариант
				auto triangle = pair.first;
				auto ind      = pair.second;
				assert(triangle->V[ind] == curPoint);
				Vector2D intersection;
				if(!GeometryUtils::IntersectionOfNonCollinearLines(triangle->E[ind].S, triangle->E[ind].E, curPoint, endPoint, intersection))
					continue;
				auto nextInd = Utils::NextInd(ind, 3);
				auto prevInd = Utils::PrevInd(ind, 3);
				     if (intersection == triangle->V[nextInd]) resList.emplace_back(make_tuple(triangle, intersection, prevInd, true ));
				else if (intersection == triangle->V[prevInd]) resList.emplace_back(make_tuple(triangle, intersection, nextInd, true ));
				else                                           resList.emplace_back(make_tuple(triangle, intersection,     ind, false));
			}
			if (resList.size() == 0)
				throw exception("Can't find intersection in all triangles");
			if ((resList.size() == 2 && !get<3>(resList[0])) || resList.size() > 2)
				throw exception("Multiple intersections in all triangles");

			curTriangle  = const_cast<Triangle*>(get<0>(resList[0])); // kludge: const_cast
			auto indEdge = get<2>(resList[0]);
			if (!get<3>(resList[0])) // структурное ребро пересекает существующее ребро
			{
				vector<Vector2D> uList;
				vector<Vector2D> lList;
				unordered_map<LineSeg, tuple<Triangle*, int, bool>> externalNeighbours;
				vector<Triangle*> newTriangles;
				vector<Triangle*> trianglesToDelete;
				auto curBegPoint = curTriangle->V[indEdge];
				auto curEndPoint = FindCavityHull(curTriangle, indEdge, startPoint, endPoint, uList, lList, externalNeighbours, trianglesToDelete);
				auto upper = TriangulatePseudopolygonDelaunay(uList, 0, static_cast<int>(uList.size()), curBegPoint, curEndPoint, newTriangles, externalNeighbours, true );
				auto lower = TriangulatePseudopolygonDelaunay(lList, 0, static_cast<int>(lList.size()), curBegPoint, curEndPoint, newTriangles, externalNeighbours, false);
				auto uTri = get<0>(upper);
				auto lTri = get<0>(lower);
				auto uInd = get<1>(upper);
				auto lInd = get<1>(lower);
				uTri->Neighbour   [uInd] = lTri;
				lTri->Neighbour   [lInd] = uTri;
				uTri->NeighbourInd[uInd] = lInd;
				lTri->NeighbourInd[lInd] = uInd;
				curTriangle = uTri;
				indEdge     = uInd;
				curPoint = curEndPoint;
				for (auto triangle : newTriangles) {
					triTable.Add(triangle);
					totalNewTriangles.push_back(triangle);
				}
				for (auto triangle : trianglesToDelete)
					triTable.ReplaceWithExistent(triangle, newTriangles[0]);
			}
			else curPoint = get<1>(resList[0]);
			curTriangle->SetStructureEdge(indEdge);
		}
		assert(CheckNeighboursMatch(totalNewTriangles));
		return totalNewTriangles;
	}
	//=============================================
	//					Statistics
	//=============================================
private:
	vector<double> angles;
	vector<double> angles1;
	vector<double> angles2;
	vector<double> qualities;
	vector<double> qualities1;
	vector<double> qualities2;

	void processStatisticsStart() {
		stat.isExtended = params.extendedStatistics;
		stat.startTime = stat.endTime = clock();
		angles .clear();
		angles1.clear();
		angles2.clear();
		qualities .clear();
		qualities1.clear();
		qualities2.clear();
		stat.minQualityDuringOptimization .clear();
		stat.meanQualityDuringOptimization.clear();
	}
	void processStatisticsAfterInputProcessing() {
		auto timeNow = clock();
		stat.timeInputProcessing = static_cast<double>(timeNow - stat.endTime) / CLOCKS_PER_SEC;
		stat.endTime = timeNow;
	}
	void processStatisticsAfterConvexDelaunay() {
		auto timeNow = clock();
		stat.timeConvexDelaunay = static_cast<double>(timeNow - stat.endTime) / CLOCKS_PER_SEC;
		stat.endTime = timeNow;

	}
	void processStatisticsAfterConstrainedDelaunay() {
		auto timeNow = clock();
		stat.timeConstainedDelaunay = static_cast<double>(timeNow - stat.endTime) / CLOCKS_PER_SEC;
		stat.endTime = timeNow;

	}
	void processStatisticsAfterAFM() {
		auto timeNow = clock();
		stat.timeAFM = static_cast<double>(timeNow - stat.endTime) / CLOCKS_PER_SEC;
		stat.endTime = timeNow;

		if (stat.isExtended) {
			stat.numTriangles = 0;

			unordered_set<Vector2D> vertexSet;
			for (auto triangle : triTable.GetSet())
			{
				if (triangle->IsExternal) continue;
				stat.numTriangles++;
				for (auto index = 0; index < 3; ++index)
					vertexSet.emplace(triangle->V[index]);
			}
			stat.numPoints = vertexSet.size();
			collectHistogramsData(angles, qualities);
		}
	}
	void processStatisticsAfterSmoothingMin() {
		auto timeNow = clock();
		stat.timeSmoothingMin = static_cast<double>(timeNow - stat.endTime) / CLOCKS_PER_SEC;
		stat.endTime = timeNow;

		if (stat.isExtended) collectHistogramsData(angles1, qualities1);
	}
	void processStatisticsAfterSmoothingMean() {
		auto timeNow = clock();
		stat.timeSmoothingMean = static_cast<double>(timeNow - stat.endTime) / CLOCKS_PER_SEC;
		stat.endTime = timeNow;

		if (stat.isExtended) collectHistogramsData(angles2, qualities2);
	}
	void processStatisticsEnd() {
		stat.endTime = clock();
		stat.timeTotal = static_cast<double>(stat.endTime - stat.startTime) / CLOCKS_PER_SEC;
		
		if (stat.isExtended){
			auto minAngle   = DBL_MAX, maxAngle   = DBL_MIN;
			auto minQuality = DBL_MAX, maxQuality = DBL_MIN;

			for (auto angle : angles) {
				minAngle = min(minAngle, angle);
				maxAngle = max(maxAngle, angle);
			}
			for (auto angle : angles1) {
				minAngle = min(minAngle, angle);
				maxAngle = max(maxAngle, angle);
			}
			for (auto angle : angles2) {
				minAngle = min(minAngle, angle);
				maxAngle = max(maxAngle, angle);
			}
			stat.minAngle = minAngle;
			stat.maxAngle = maxAngle;

			for (auto quality : qualities) {
				minQuality = min(minQuality, quality);
				maxQuality = max(maxQuality, quality);
			}
			for (auto quality : qualities1) {
				minQuality = min(minQuality, quality);
				maxQuality = max(maxQuality, quality);
			}
			for (auto quality : qualities2) {
				minQuality = min(minQuality, quality);
				maxQuality = max(maxQuality, quality);
			}
			stat.minQuality = minQuality;
			stat.maxQuality = maxQuality;

			stat.histTriangleAngles   .clear();
			stat.histTriangleAnglesAO1.clear();
			stat.histTriangleAnglesAO2.clear();
			auto anglesK = stat.displaySize / (maxAngle - minAngle);
			for (auto i = 0; i < stat.displaySize; ++i) {
				Vector2D p = { (i + 0.5) / anglesK + minAngle, 0 };
				stat.histTriangleAngles   .emplace_back(p);
				stat.histTriangleAnglesAO1.emplace_back(p);
				stat.histTriangleAnglesAO2.emplace_back(p);
			}
			for (auto angle : angles ) stat.histTriangleAngles   [min(stat.displaySize - 1, static_cast<int>((angle - minAngle)*anglesK))].Y++;
			for (auto angle : angles1) stat.histTriangleAnglesAO1[min(stat.displaySize - 1, static_cast<int>((angle - minAngle)*anglesK))].Y++;
			for (auto angle : angles2) stat.histTriangleAnglesAO2[min(stat.displaySize - 1, static_cast<int>((angle - minAngle)*anglesK))].Y++;

			stat.histTriangleQualities   .clear();
			stat.histTriangleQualitiesAO1.clear();
			stat.histTriangleQualitiesAO2.clear();
			auto qualitiesK = stat.displaySize / (maxQuality - minQuality);
			for (auto i = 0; i < stat.displaySize; ++i) {
				Vector2D p = { (i + 0.5) / qualitiesK + minQuality, 0 };
				stat.histTriangleQualities   .emplace_back(p);
				stat.histTriangleQualitiesAO1.emplace_back(p);
				stat.histTriangleQualitiesAO2.emplace_back(p);
			}
			for (auto quality : qualities ) stat.histTriangleQualities   [min(stat.displaySize - 1, static_cast<int>((quality - minQuality)*qualitiesK))].Y++;
			for (auto quality : qualities1) stat.histTriangleQualitiesAO1[min(stat.displaySize - 1, static_cast<int>((quality - minQuality)*qualitiesK))].Y++;
			for (auto quality : qualities2) stat.histTriangleQualitiesAO2[min(stat.displaySize - 1, static_cast<int>((quality - minQuality)*qualitiesK))].Y++;
		}
	}
	void collectHistogramsData(vector<double>& anglesHistograms, vector<double>& qualitiesHistograms) const {
		for (auto triangle : triTable.GetSet())
		{
			if (triangle->IsExternal) continue;
			if (triangle->IsStructureEdge(0) && triangle->IsStructureEdge(1) ||
				triangle->IsStructureEdge(1) && triangle->IsStructureEdge(2) ||
				triangle->IsStructureEdge(2) && triangle->IsStructureEdge(0))
				continue;
			qualitiesHistograms.emplace_back(abs(triangle->AlphaQuality));
			for (auto index = 0; index < 3; ++index)
				anglesHistograms.emplace_back(acos(triangle->Cos(index))*180.0/M_PI);
		}
	}
	void processStatisticsDuringOptimization(int ind) {
		// todo: valid processing of structure triangles
		auto minQuality = DBL_MAX;
		auto meanQuality = 1.0;
		auto triCounter = 0;
		for (auto triangle : triTable.GetSet())
		{
			if (triangle->IsExternal) continue;
			if (triangle->IsStructureEdge(0) && triangle->IsStructureEdge(1) ||
				triangle->IsStructureEdge(1) && triangle->IsStructureEdge(2) ||
				triangle->IsStructureEdge(2) && triangle->IsStructureEdge(0))
				continue;
			auto q = abs(triangle->AlphaQuality);
			minQuality = min(minQuality, q);
			meanQuality *= q;
			triCounter++;
		}
		meanQuality = pow(meanQuality, 1.0 / triCounter);
		stat. minQualityDuringOptimization.emplace_back(Vector2D{ static_cast<double>(ind),  minQuality });
		stat.meanQualityDuringOptimization.emplace_back(Vector2D{ static_cast<double>(ind), meanQuality });
	}
	//--------------------------------------
	//				Checking
	//--------------------------------------
public: 
	template<class T>
	static bool CheckValidDelaunay(const T& triangles)
	{
		for(auto triangle: triangles)
			for (auto i = 0; i < 3; ++i)
				if (!DelaunayCheck(triangle, i))
					return false;
		return true;
	}
private:
	static bool InCircle(const Triangle* tri, const Vector2D& point)
	{
		return GeometryUtils::InCircle(tri->V[0], tri->V[1], tri->V[2], point, tri->Orientation == PointsOrder::Clockwise);
	}
	static bool DelaunayCheck(const Triangle* tri, int ind)
	{
		if (tri->IsStructureEdge(ind)) return true; // структурные ребра не участвуют в проверке условия делоне
		auto nTri = tri->Neighbour[ind];
		return !InCircle(tri, nTri->V[tri->NeighbourInd[ind]]) && !InCircle(nTri, tri->V[ind]);
	}
	template <class T>
	static bool CheckNeighboursMatch(const T& triangles)
	{
		for (const Triangle* triangle : triangles)
			if (!CheckNeighbourMatch(triangle))
				return false;
		return true;
	}
	static bool CheckNeighbourMatch(const Triangle* tri)
	{
		if (tri == nullptr) return true;
		for (auto ind = 0; ind < 3; ++ind)
		{
			auto nTri = tri->Neighbour[ind];
			auto nind = tri->NeighbourInd[ind];
			if (nTri == nullptr) continue;
			if (nTri->Neighbour[nind] != tri) 
				return false;
			if (tri->Neighbour[nTri->NeighbourInd[nind]] != nTri) 
				return false;
			if (tri->E[ind] != nTri->E[nind])
				return false;
			if (tri->E[ind].S != tri->NextV(ind) ||
				tri->E[ind].E != tri->PrevV(ind))
				return false;
		}
		return true;
	}
	template <class T>
	static bool CheckStructureEdges(const T& triangles, const unordered_set<LineSeg>& structureEdges)
	{
		unordered_set<LineSeg> structure_segment_set;
		for (auto t : triangles)
			for (auto ind = 0; ind<3; ++ind)
				if (t->isStructureEdge[ind]) {
					if (structureEdges.find(t->E[ind]) == structureEdges.end())
						return false;
					structure_segment_set.emplace(t->E[ind]);
				}
		return structure_segment_set.size() == structureEdges.size();
	}
};


Mesh::~Mesh() {}
Mesh::Mesh(): d(make_unique<MeshImpl>()){}

MeshingResult Mesh::generateMesh(const FragmentGeometry& cg, const MeshingParameters& params) const
{
	d->params = params;
	d->cg = cg;
	d->Generate();
	return d->result;
}

vector<IntermediateStep> Mesh::getIntermediateSteps() const { return d->intermediateSteps; }
MeshingStatistics Mesh::getStatistics() const { 
	// todo: debug only
	SaveMeshingStatisticsToTxtFile(d->stat, "statitics.txt");
	return d->stat; 
}

void SaveMeshingStatisticsToTxtFile(const MeshingStatistics & stat, const char * filename)
{
	ofstream ofs(filename);
	ofs << "Min quality during optimization:\n";
	for (auto& p : stat.minQualityDuringOptimization)
		ofs << "  " << p.X << "; " << p.Y << "\n";
	ofs << "Mean quality during optimization:\n";
	for (auto& p : stat.meanQualityDuringOptimization)
		ofs << "  " << p.X << "; " << p.Y << "\n";

	ofs << "Angles before optimization:\n";
	for (auto& p : stat.histTriangleAngles)
		ofs << "  " << p.X << "; " << p.Y << "\n";
	ofs << "Angles after optimization(I):\n";
	for (auto& p : stat.histTriangleAnglesAO1)
		ofs << "  " << p.X << "; " << p.Y << "\n";
	ofs << "Angles after optimization(II):\n";
	for (auto& p : stat.histTriangleAnglesAO2)
		ofs << "  " << p.X << "; " << p.Y << "\n";

	ofs << "Qualities before optimization:\n";
	for (auto& p : stat.histTriangleQualities)
		ofs << "  " << p.X << "; " << p.Y << "\n";
	ofs << "Qualities after optimization(I):\n";
	for (auto& p : stat.histTriangleQualitiesAO1)
		ofs << "  " << p.X << "; " << p.Y << "\n";
	ofs << "Qualities after optimization(II):\n";
	for (auto& p : stat.histTriangleQualitiesAO2)
		ofs << "  " << p.X << "; " << p.Y << "\n";

	ofs << "Optimization(I) iterations: "  << stat.minQualityIterations << "\n";
	ofs << "Optimization(II) iterations: " << stat.meanQualityIterations << "\n";

	ofs << "Number of points: " << stat.numPoints << "\n";
	ofs << "Number of triangles: " << stat.numTriangles << "\n";
}

};


//void process_input(unordered_set<Vector2D>& structurePoints, unordered_set<LineSeg>& structureEdges)
//{
//	if (boundary_points.size() < 3)
//		throw exception("Not enought points to preform triangulation");
//
//	bounds = BoundingBox(boundary_points);
//	if (bounds.DeltaX() * bounds.DeltaY() <= 0.0)
//		throw exception("Invalid points to preform triangulation");
//
//	structurePoints.clear();
//	structureEdges.clear();
//	//vector<LineSeg> edgesQueue;
//
//	// todo: process holes
//	auto prevPoint = boundary_points.back();
//	for (auto const& point : boundary_points) {
//		structurePoints.emplace(point);
//		structureEdges.emplace(LineSeg(prevPoint, point));
//		//edgesQueue.emplace_back(LineSeg(prevPoint, point));
//		prevPoint = point;
//	}
//}
//bool process_boundary(unordered_set<Vector2D>& structurePoints, const unordered_set<LineSeg>& structure_segment_set)
//{
//	std::vector<pair<Triangle*, int>> queue;
//	for (auto t : triTable.GetSet())
//		for (auto ind = 0; ind<3; ++ind)
//			if (t->isStructureEdge[ind])
//				queue.emplace_back(t, ind);
//
//	//  for debugging purposes
//	std::sort(queue.begin(), queue.end(), [](const pair<Triangle*, int>& a, const pair<Triangle*, int>& b) { return a.first->AlphaQuality < b.first->AlphaQuality; });
//
//	std::vector<Triangle*> newTriangles;
//	while (queue.size()>0)
//	{
//		auto pairTI = queue.back(); queue.pop_back();
//		auto t = pairTI.first;
//		auto i = pairTI.second;
//		if (!triTable.IsAlive(t)) continue;
//		if (structure_segment_set.find(t->E[i]) == structure_segment_set.end()) continue;
//
//		auto internal_points = divide_edge(t->E[i]);
//		if (internal_points.size() == 2) continue;
//
//		assert(internal_points[0] == t->NextV(i) && internal_points.size() > 2);
//
//		if (t->Neighbour[i] == nullptr)
//		{
//			for (size_t ind = 1; ind<internal_points.size(); ++ind)
//			{
//				auto t1 = new Triangle(internal_points[ind - 1], internal_points[ind], t->V[i]);
//				t1->isStructureEdge[2] = true;
//
//				if (ind == 1)
//				{
//					auto ci = Utils::PrevInd(i, 3);
//					assert(t->E[ci].Contains(t->E[i].S));
//					UpdateReference(t, t1, ci, 1);
//					if (t->IsStructureEdge(ci)) queue.emplace_back(t1, 1);
//				}
//				else
//				{
//					auto l1 = newTriangles.back();
//					l1->Neighbour[0] = t1;
//					t1->Neighbour[1] = l1;
//					l1->NeighbourInd[0] = 1;
//					t1->NeighbourInd[1] = 0;
//				}
//
//				triTable.Add(t1);
//				newTriangles.emplace_back(t1);
//				structurePoints.emplace(internal_points[ind]);
//			}
//			auto fi = Utils::NextInd(i, 3);
//			assert(t->E[fi].Contains(t->E[i].E));
//			UpdateReference(t, newTriangles.back(), fi, 0);
//			if (t->IsStructureEdge(fi)) queue.emplace_back(newTriangles.back(), 0);
//
//			// окончательно удаляем все сведения об исходных треугольниках
//			triTable.ReplaceWithExistent(t, newTriangles.back());
//
//			if (params.saveIntermediateSteps) SaveIntermediateStep();
//		}
//		else
//		{
//			auto nt = t->Neighbour[i];
//			auto ni = t->NeighbourInd[i];
//			for (size_t ind = 1; ind<internal_points.size(); ++ind)
//			{
//				auto t1 = new Triangle(internal_points[ind - 1], internal_points[ind], t->V[i]);
//				auto t2 = new Triangle(internal_points[ind - 1], internal_points[ind], nt->V[ni]);
//				t1->IsExternal = t->IsExternal;
//				t2->IsExternal = nt->IsExternal;
//
//				// обновляем ссылки друг на друга
//				t1->Neighbour[2] = t2;
//				t2->Neighbour[2] = t1;
//				t1->NeighbourInd[2] = t2->NeighbourInd[2] = 2;
//				t1->isStructureEdge[2] = t2->isStructureEdge[2] = true;
//
//				if (ind == 1)
//				{
//					auto ci = Utils::PrevInd(i, 3);
//					assert(t->E[ci].Contains(t->E[i].S));
//					UpdateReference(t, t1, ci, 1);
//					if (t->IsStructureEdge(ci)) queue.emplace_back(t1, 1);
//
//					auto nci = t->Orientation == nt->Orientation ? Utils::NextInd(ni, 3) : Utils::PrevInd(ni, 3);
//					assert(nt->E[nci].Contains(t->E[i].S));
//					UpdateReference(nt, t2, nci, 1);
//					if (nt->IsStructureEdge(nci)) queue.emplace_back(t2, 1);
//				}
//				else
//				{
//					auto l1 = newTriangles[newTriangles.size() - 2];
//					auto l2 = newTriangles[newTriangles.size() - 1];
//					l1->Neighbour[0] = t1;
//					t1->Neighbour[1] = l1;
//					l2->Neighbour[0] = t2;
//					t2->Neighbour[1] = l2;
//					l1->NeighbourInd[0] = l2->NeighbourInd[0] = 1;
//					t1->NeighbourInd[1] = t2->NeighbourInd[1] = 0;
//				}
//
//				triTable.Add(t1);
//				triTable.Add(t2);
//				newTriangles.emplace_back(t1);
//				newTriangles.emplace_back(t2);
//				structurePoints.emplace(internal_points[ind]);
//			}
//			auto fi = Utils::NextInd(i, 3);
//			auto nfi = t->Orientation == nt->Orientation ? Utils::PrevInd(ni, 3) : Utils::NextInd(ni, 3);
//			assert(t->E[fi].Contains(t->E[i].E));
//			assert(nt->E[nfi].Contains(t->E[i].E));
//			UpdateReference(t, newTriangles[newTriangles.size() - 2], fi, 0);
//			UpdateReference(nt, newTriangles[newTriangles.size() - 1], nfi, 0);
//			if (t->IsStructureEdge(fi)) queue.emplace_back(newTriangles[newTriangles.size() - 2], 0);
//			if (nt->IsStructureEdge(nfi)) queue.emplace_back(newTriangles[newTriangles.size() - 1], 0);
//
//			// окончательно удаляем все сведения об исходных треугольниках
//			triTable.ReplaceWithExistent(t, newTriangles[newTriangles.size() - 2]);
//			triTable.ReplaceWithExistent(nt, newTriangles[newTriangles.size() - 1]);
//
//			if (params.saveIntermediateSteps) SaveIntermediateStep();
//		}
//	}
//
//	assert(CheckNeighboursMatch(triTable.GetSet()));
//	FlipAllInconditionalTriangles(newTriangles);
//	return true;
//}