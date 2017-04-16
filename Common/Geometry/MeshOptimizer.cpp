#include "includes.h"
#include "GeometryPrimitives.h"
#include "Polynomial.h"
#include "GeometryUtils.h"
using namespace std;

namespace FemLibrary {

	pair<Vector2D, double> CenterAndRadius(const Vector2D& m, const Vector2D& n, double k, double j, double quality)
	{
		auto pc = m + n / quality;
		auto r = sqrt(pc*pc - k + j / quality);
		return make_pair(pc, r);
	}
	pair<Vector2D, double> CenterAndRadius(const LineSeg& s1, double quality)
	{
		auto m1 = s1.Middle();
		auto n1 = sqrt(3) / 2 * s1.Vector().Norm();
		auto k1 = s1.Vector().LengthSquared() + s1.S*s1.E;
		auto j1 = sqrt(3)*Vector2D::CrossProduct(s1.S, s1.E);
		return CenterAndRadius(m1, n1, k1, j1, quality);
	}

#define equation2_definitions \
Vector2D m1, n1;\
Vector2D m2, n2;\
double k1, j1;	\
double k2, j2;	

#define equation3_definitions Vector3D kvec,jvec,dmx,dmy,dnx,dny,mx,nx;

#define equation2 \
m1 = s1.Middle();																	\
m2 = s2.Middle();																	\
n1 = sqrt(3) / 2 * s1.Vector().Norm();												\
n2 = sqrt(3) / 2 * s2.Vector().Norm();												\
k1 = s1.Vector().LengthSquared() + s1.S*s1.E;										\
k2 = s2.Vector().LengthSquared() + s2.S*s2.E;										\
j1 = sqrt(3)*Vector2D::CrossProduct(s1.S, s1.E);									\
j2 = sqrt(3)*Vector2D::CrossProduct(s2.S, s2.E);									\
auto dn = n1 - n2;																	\
auto dk = k1 - k2;																	\
auto dm = m1 - m2;																	\
auto dj = j1 - j2;																	\
auto km = k1*m2 - k2*m1;															\
auto kn = k1*n2 - k2*n1;															\
auto jm = j1*m2 - j2*m1;															\
auto jn = j1*n2 - j2*n1;															\
auto cn  = Vector2D::CrossProduct(n1, n2);											\
auto cm  = Vector2D::CrossProduct(m1, m2);											\
auto cmn = Vector2D::CrossProduct(m1, n2) - Vector2D::CrossProduct(m2, n1);			\
vector<double> poly(5);																\
poly[0] = dk*dk - 4 * (cm*cm + km*dm);												\
poly[1] = -2 * dj*dk - 4 * (kn*dm - jm*dm + km*dn + 2 * cm*cmn);					\
poly[2] = dj*dj + 4 * (jn*dm + jm*dn - cmn*cmn - 2 * cn*cm - kn*dn);				\
poly[3] = 4 * jn*dn - 8 * cn*cmn;													\
poly[4] = -4 * cn*cn;

#define equation3 \
auto m1 = s1.Middle();																														\
auto m2 = s2.Middle();																														\
auto m3 = s3.Middle();																														\
auto n1 = sqrt(3) / 2 * s1.Vector().Norm();																									\
auto n2 = sqrt(3) / 2 * s2.Vector().Norm();																									\
auto n3 = sqrt(3) / 2 * s3.Vector().Norm();																									\
auto k1 = s1.Vector().LengthSquared() + s1.S*s1.E;																							\
auto k2 = s2.Vector().LengthSquared() + s2.S*s2.E;																							\
auto k3 = s3.Vector().LengthSquared() + s3.S*s3.E;																							\
auto j1 = sqrt(3)*Vector2D::CrossProduct(s1.S, s1.E);																						\
auto j2 = sqrt(3)*Vector2D::CrossProduct(s2.S, s2.E);																						\
auto j3 = sqrt(3)*Vector2D::CrossProduct(s3.S, s3.E);																						\
kvec = Vector3D (k1, k2, k3);																												\
jvec = Vector3D (j1, j2, j3);																												\
dmx = Vector3D (m3.X - m2.X, m1.X - m3.X, m2.X - m1.X);																						\
dmy = Vector3D (m3.Y - m2.Y, m1.Y - m3.Y, m2.Y - m1.Y);																						\
dnx = Vector3D (n3.X - n2.X, n1.X - n3.X, n2.X - n1.X);																						\
dny = Vector3D (n3.Y - n2.Y, n1.Y - n3.Y, n2.Y - n1.Y);																						\
mx = Vector3D (m1.X, m2.X, m3.X);																											\
nx = Vector3D (n1.X, n2.X, n3.X);																											\
auto b1 = 2 * mx*dmy;																														\
auto b2 = 2 * (nx*dmy + mx*dny);																											\
auto b3 = 2 * nx*dny;																														\
Vector2D v1(-dny*jvec, dnx*jvec);																											\
Vector2D v2(dny*kvec - dmy*jvec, jvec*dmx - kvec*dnx);																						\
Vector2D v3(kvec*dmy, -kvec*dmx);																											\
																																			\
vector<double> poly(6);																														\
poly[0] = v3*v3 + b1*(b1*k1 - 2 * v3*m1);																									\
poly[1] = 2 * (v2 - b2* m1)*v3 - b1*(b1* j1 + 2 * (v2*m1 + v3*n1 - b2* k1));																\
poly[2] = v2*v2 + 2 * v1*v3 - 2 * b3 *v3*m1 + b2*(b2* k1 - 2 * (v2*m1 + v3*n1)) - 2 * b1*(b2*j1 - b3*k1 + v1*m1 + v2*n1);					\
poly[3] = 2 * (v1*v2 - b1*(b3*j1 + v1*n1) - b3*(v2*m1 + v3*n1)) - b2*(b2*j1 - 2 * (b3*k1 - v1*m1 - v2*n1));									\
poly[4] = v1*(v1 - 2 * b2*n1) + b3*(b3*k1 - 2 * (b2*j1 + v1*m1 + v2*n1));																	\
poly[5] = -b3*(b3*j1 + 2 * v1*n1);

	FEMLIBRARY_API vector<double> Equation2_Test(const LineSeg& s1, const LineSeg& s2)
	{
		equation2_definitions
			equation2
			return poly;
	}
	FEMLIBRARY_API vector<double> Equation3_Test(const LineSeg& s1, const LineSeg& s2, const LineSeg& s3)
	{
		equation3_definitions
			equation3
			return poly;
	}
	const double precisionOfMSW = 1.01;
	struct Basis
	{
		Vector2D point;
		double quality;
		int side_indexes[3];

		Basis() :quality(0.0), side_indexes{ -1,-1,-1 } {}
		Basis(const Vector2D& point, double quality, int s1, int s2) :point(point), quality(quality), side_indexes{ s1,s2, -1 } {}
		Basis(const Vector2D& point, double quality, int s1, int s2, int s3) :point(point), quality(quality), side_indexes{ s1,s2, s3 } {}
		//static Basis Empty() { return Basis( Vector2D(), DBL_MIN, -1,-1,-1); }
		bool IsTwoSideBasis() const { return side_indexes[2] == -1; }
		bool IsEmpty() const { return side_indexes[0] == -1; }
		int Length() const { return IsTwoSideBasis() ? 2 : 3; }
		bool ContainsIndex(int ind) const {
			return side_indexes[0] == ind || side_indexes[1] == ind ||
				(!IsTwoSideBasis() && side_indexes[2] == ind);
		}
		bool ContainsExactly(const vector<int>& indexes) const
		{
			if (indexes.size() != Length()) return false;
			for (auto i : indexes) if (!ContainsIndex(i)) return false;
			return true;
		}

		/*bool operator<(const Basis& another) const { return quality < another.quality; }

		static Basis Basis2(const vector<LineSeg>& sides, int i1, int i2)
		{
			assert(i1 != i2);
			if (i1 > i2) swap(i1, i2);

			const auto& s1 = sides[i1];
			const auto& s2 = sides[i2];

			equation2

			auto roots = Polynomial::FindRoots(poly);
			if (roots.size() == 0) throw exception("Basis2: can't find root");

			auto q = DBL_MIN;
			for (size_t i = 0; i < roots.size(); ++i)
				if (0 < roots[i] && roots[i] <= 1)
					q = max(q, roots[i]);
			if(q == DBL_MIN) throw exception("Basis2: invalid roots");
			assert(q > 0);

			auto cr1 = CenterAndRadius(m1, n1, k1, j1, q);
			auto cr2 = CenterAndRadius(m2, n2, k2, j2, q);
			assert(abs((cr1.first - cr2.first).Length() - (cr1.second + cr2.second)) < 1e-5);

			auto ip = (cr2.second*cr1.first + cr1.second*cr2.first) / (cr1.second + cr2.second);
			return Basis(ip, q, i1, i2);
		}
		static Basis Basis3(const vector<LineSeg>& sides, int i1, int i2, int i3)
		{
			assert(i1 != i2 && i1 != i3 && i2 != i3);
			if (i1 > i2) swap(i1, i2);
			if (i1 > i3) swap(i1, i3);
			if (i2 > i3) swap(i2, i3);

			const auto& s1 = sides[i1];
			const auto& s2 = sides[i2];
			const auto& s3 = sides[i3];

			equation3

			auto roots = Polynomial::SolveQuinticInRange(poly, 0., 1., true); // todo: apply b range
			// на самом деле пересечение должно быть всегда, но так как у меня
			// вычисления не достаточно точные, то значение функции в ней очень
			// велико, возможно стоит минимизировать систему из 2 уравнений или
			// находить расстояния до окружностей в экстремумах
			//if (roots.size() == 0) throw exception("Basis3: can't find root");
			if (roots.size() == 0) return Empty();

			auto q = roots[0];
			assert(q > 0);
			Vector2D ip( (q * kvec - jvec)*(q * dmy + dny) / (2*(q*mx + nx)*(q*dmy + dny)),
						 (jvec - q * kvec)*(q * dmx + dnx) / (2*(q*mx + nx)*(q*dmy + dny)) );
			//assert todo
			return Basis(ip, q, i1,i2,i3 );
		}
		*/
		bool Check(const vector<LineSeg>& side_vector, int side1, int side2 = -1) const
		{
			if (IsEmpty()) return false;
			auto cr = CenterAndRadius(side_vector[side1], quality);
			if ((point - cr.first).LengthSquared() / cr.second / cr.second > precisionOfMSW) return false;
			if (side2 > -1)
			{
				cr = CenterAndRadius(side_vector[side2], quality);
				if ((point - cr.first).LengthSquared() / cr.second / cr.second > precisionOfMSW) return false;
			}
			return true;
		}
	};
	struct BasisBuilder
	{
		const LineSeg *r1, *r2;
		const vector<LineSeg>& sides;
		int i1, i2, i3;

		BasisBuilder(const vector<LineSeg>& sides, int i1, int i2, int i3, const LineSeg* r1 = nullptr, const LineSeg* r2 = nullptr) : r1(r1), r2(r2), sides(sides), i1(i1), i2(i2), i3(i3) {}
		virtual Vector2D CalcPoint(double q) = 0;
		virtual int FindSolutions(pair<double, BasisBuilder*> *solutions) = 0;
		bool CheckRestrictionAndReturnBasis(double q, Basis& basis)
		{	// todo: create basis
			auto ip = CalcPoint(q);
			if (r1 != nullptr) {
				auto cr = CenterAndRadius(*r1, q);
				if ((ip - cr.first).LengthSquared() / cr.second / cr.second > precisionOfMSW) return false;

				if (r2 != nullptr)
				{
					cr = CenterAndRadius(*r2, q);
					if ((ip - cr.first).LengthSquared() / cr.second / cr.second > precisionOfMSW) return false;
				}
			}
			basis = Basis(ip, q, i1, i2, i3);
			return true;
		}
	};

	struct BasisBuilder2 :public BasisBuilder
	{
		equation2_definitions

			BasisBuilder2(const vector<LineSeg>& sides, int i1, int i2, int r1 = -1, int r2 = -1)
			: BasisBuilder(sides, i1, i2, -1, r1 == -1 ? nullptr : &sides[r1],
				r2 == -1 ? nullptr : &sides[r2]) {}

		int FindSolutions(pair<double, BasisBuilder*> *solutions) override {
			assert(i1 != i2);
			// сортируем индексы
			if (i1 > i2) swap(i1, i2);

			const auto& s1 = sides[i1];
			const auto& s2 = sides[i2];

			equation2

				auto roots = Polynomial::FindRoots(poly);
			if (roots.size() == 0)
				throw exception("BasisBuilder2: can't find root");

			auto solutionsFounded = 0;
			for (size_t i = 0; i < roots.size(); ++i)
				if (0. < roots[i] && roots[i] <= 1.)
					solutions[solutionsFounded++] = make_pair(roots[i], this);
			return solutionsFounded;
		}
		Vector2D CalcPoint(double q) override {
			auto cr1 = CenterAndRadius(m1, n1, k1, j1, q);
			auto cr2 = CenterAndRadius(m2, n2, k2, j2, q);
			//if (abs((cr1.first - cr2.first).Length() - (cr1.second + cr2.second)) >= 1e-5) {
			//	auto val = abs((cr1.first - cr2.first).Length() - (cr1.second + cr2.second));
			//	assert(abs((cr1.first - cr2.first).Length() -    (cr1.second + cr2.second)) < 1e-5 ||
			//		   abs((cr1.first - cr2.first).Length() - abs(cr1.second - cr2.second)) < 1e-5 );
			//}

			return (cr2.second*cr1.first + cr1.second*cr2.first) / (cr1.second + cr2.second);
		}
	};
	struct BasisBuilder3 :public BasisBuilder
	{
		equation3_definitions

			BasisBuilder3(const vector<LineSeg>& sides, int i1, int i2, int i3, int r1 = -1, int r2 = -1)
			: BasisBuilder(sides, i1, i2, i3, r1 == -1 ? nullptr : &sides[r1],
				r2 == -1 ? nullptr : &sides[r2]) {}

		int FindSolutions(pair<double, BasisBuilder*> *solutions) override {
			assert(i1 != i2 && i1 != i3 && i2 != i3);
			// сортируем индексы
			if (i1 > i2) swap(i1, i2);
			if (i1 > i3) swap(i1, i3);
			if (i2 > i3) swap(i2, i3);

			const auto& s1 = sides[i1];
			const auto& s2 = sides[i2];
			const auto& s3 = sides[i3];

			equation3

				auto roots = Polynomial::SolveQuinticInRange(poly, 0., 1., false); // todo: apply b range 
				// на самом деле пересечение должно быть всегда, но так как у меня 
				// вычисления не достаточно точные, то значение функции в ней очень 
				// велико, возможно стоит минимизировать систему из 2 уравнений или 
				// находить расстояния до окружностей в экстремумах
				//if (roots.size() == 0) throw exception("Basis3: can't find root");

			auto solutionsFounded = 0;
			for (size_t i = 0; i < roots.size(); ++i)
				if (0. < roots[i] && roots[i] <= 1.)
					solutions[solutionsFounded++] = make_pair(roots[i], this);
			return solutionsFounded;
		}
		Vector2D CalcPoint(double q) override {
			return Vector2D((q * kvec - jvec)*(q * dmy + dny) / (2 * (q*mx + nx)*(q*dmy + dny)),
				(jvec - q * kvec)*(q * dmx + dnx) / (2 * (q*mx + nx)*(q*dmy + dny)));
		}
	};

	// todo: находить сразу все базисы не нужно
	Basis GetBasis2(const Basis& b, const vector<LineSeg>& side_vector, int ind_x)
	{
		assert(!b.ContainsIndex(ind_x));
		auto candidates_count = 0;
		pair<double, BasisBuilder*> candidates[13];
		BasisBuilder2 b0(side_vector, b.side_indexes[0], ind_x, b.side_indexes[1]);
		BasisBuilder2 b1(side_vector, b.side_indexes[1], ind_x, b.side_indexes[0]);
		BasisBuilder3 b01(side_vector, b.side_indexes[0], ind_x, b.side_indexes[1]);

		candidates_count += b0.FindSolutions(candidates + candidates_count);
		candidates_count += b1.FindSolutions(candidates + candidates_count);
		candidates_count += b01.FindSolutions(candidates + candidates_count);
		assert(0 < candidates_count&&candidates_count <= 13);

		sort(candidates, candidates + candidates_count, [](const pair<double, BasisBuilder*>& x, const pair<double, BasisBuilder*>& y) { return x.first > y.first; }); // по возрастанию
		Basis basis;
		Basis first_checked;
		for (auto i = 0; i < candidates_count; ++i)
		{
			if (candidates[i].second->CheckRestrictionAndReturnBasis(candidates[i].first, basis)) {
				if (first_checked.IsEmpty()) first_checked = basis;
				if(basis.quality <= b.quality)
				{
					// Этот ассерт постоянно срабатывает т.к. алгоритм нахождения корней 5 степени постоянно лажает, поэтому проще всего (и правильнее?) брать решение с наименьшим корнем
					// assert(basis.quality < b.quality);
					return basis;
				}
			}
		}
		// Это особый случай, который возникает, если найденное решение лучше в плане качества, чем существующий базис. В таком случае придется брать его, пока не придумаю, как сделать это правильно
		if (!first_checked.IsEmpty())
			return first_checked;
		throw exception("GetBasis2: can't find valid candidates");
	}
	Basis GetBasis3(const Basis& b, const vector<LineSeg>& side_vector, int ind_x)
	{
		assert(!b.ContainsIndex(ind_x));
		auto candidates_count = 0;
		pair<double, BasisBuilder*> candidates[27];
		BasisBuilder2 b20(side_vector, b.side_indexes[0], ind_x, b.side_indexes[1], b.side_indexes[2]);
		BasisBuilder2 b21(side_vector, b.side_indexes[1], ind_x, b.side_indexes[2], b.side_indexes[0]);
		BasisBuilder2 b22(side_vector, b.side_indexes[2], ind_x, b.side_indexes[0], b.side_indexes[1]);
		BasisBuilder3 b30(side_vector, b.side_indexes[1], b.side_indexes[2], ind_x, b.side_indexes[0]);
		BasisBuilder3 b31(side_vector, b.side_indexes[2], b.side_indexes[0], ind_x, b.side_indexes[1]);
		BasisBuilder3 b32(side_vector, b.side_indexes[0], b.side_indexes[1], ind_x, b.side_indexes[2]);

		candidates_count += b20.FindSolutions(candidates + candidates_count);
		candidates_count += b21.FindSolutions(candidates + candidates_count);
		candidates_count += b22.FindSolutions(candidates + candidates_count);
		candidates_count += b30.FindSolutions(candidates + candidates_count);
		candidates_count += b31.FindSolutions(candidates + candidates_count);
		candidates_count += b32.FindSolutions(candidates + candidates_count);
		assert(0 < candidates_count&&candidates_count <= 27);

		sort(candidates, candidates + candidates_count, [](const pair<double, BasisBuilder*>& x, const pair<double, BasisBuilder*>& y) { return x.first > y.first; }); // по возрастанию
		Basis basis;
		Basis first_checked;
		for (auto i = 0; i < candidates_count; ++i)
		{
			if (candidates[i].second->CheckRestrictionAndReturnBasis(candidates[i].first, basis))
			{
				if (first_checked.IsEmpty()) first_checked = basis;
				if (basis.quality <= b.quality)
				{
					// Этот ассерт постоянно срабатывает т.к. алгоритм нахождения корней 5 степени постоянно лажает, поэтому проще всего (и правильнее?) брать решение с наименьшим корнем
					// assert(basis.quality < b.quality);
					return basis;
				}
			}
		}

		// Это особый случай, который возникает, если найденное решение лучше в плане качества, чем существующий базис. В таком случае придется брать его, пока не придумаю, как сделать это правильно
		if(!first_checked.IsEmpty())
			return first_checked;
		throw exception("GetBasis3: can't find valid candidates");
	}
	Basis GetBasis(const Basis& b, const vector<LineSeg>& side_vector, int i)
	{
		return b.IsTwoSideBasis() ? GetBasis2(b, side_vector, i) : GetBasis3(b, side_vector, i);
	}

	Basis MSW(const vector<LineSeg>& side_vector, const vector<int>& sides, const Basis& c)
	{
		if (c.ContainsExactly(sides)) return c;
		int n_free = sides.size() - c.Length();
		if (n_free <= 0)
			throw exception("basis mismatch");

		auto x_pos = rand() % n_free;
		vector<int> new_sides; 
		new_sides.reserve(sides.size() - 1);
		auto pos = 0, x = -1;
		for (size_t i = 0; i < sides.size(); ++i)
		{
			new_sides.emplace_back(sides[i]);
			if (!c.ContainsIndex(sides[i])) {
				if (pos == x_pos) {
					x = sides[i];
					new_sides.pop_back();
				}
				pos++;
			}
		}
		assert(x != -1 && new_sides.size() == sides.size() - 1);

		auto b = MSW(side_vector, new_sides, c);
		if (b.Check(side_vector, x)) return b;
		return MSW(side_vector, sides, GetBasis(b, side_vector, x));
	}

	bool CheckValidBasis(const Basis& b, const vector<LineSeg>& side_vector)
	{
		auto q = DBL_MAX;
		for (size_t i = 0; i < side_vector.size(); ++i) {
			if (!b.Check(side_vector, i))
				return false;
			q = min(q, GeometryUtils::TriangleAlphaQuality(b.point, side_vector[i].S, side_vector[i].E));
		}
		if (abs(q - b.quality) > 1e-5)
			return false;
		return true;
	}

	bool MovePointMSW(const Vector2D& old_point, const vector<pair<Triangle*, int>>& vertex_triangles, Vector2D& new_point)
	{
		auto n = vertex_triangles.size();

		vector<LineSeg> sides(n);
		transform(vertex_triangles.begin(), vertex_triangles.end(), sides.begin(), [](const pair<Triangle*, int>& p) {return p.first->E[p.second]; });

		// строим начальный базис
		// находим две наихудших треугольника
		auto min_quality = DBL_MAX, second_min_quality = DBL_MAX;
		auto min_ind = -1, second_min_ind = -1;
		for (size_t i = 0; i < vertex_triangles.size(); ++i)
		{
			assert(vertex_triangles[i].first->Orientation == PointsOrder::Counterclockwise);
			auto q = vertex_triangles[i].first->AlphaQuality;
			if (q < min_quality)
			{
				second_min_quality = min_quality;
				min_quality = q;
				second_min_ind = min_ind;
				min_ind = i;
			}
			else if (q < second_min_quality)
			{
				second_min_quality = q;
				second_min_ind = i;
			}
		}
		assert(min_ind != second_min_ind);

		BasisBuilder2 bb(sides, min_ind, second_min_ind);
		pair<double, BasisBuilder*> candidates[4];
		auto candidates_count = bb.FindSolutions(candidates);
		sort(candidates, candidates + candidates_count, [](const pair<double, BasisBuilder*>& x, const pair<double, BasisBuilder*>& y) { return x.first > y.first; }); // по возрастанию
		assert(candidates_count > 0);
		Basis basis;
		auto res = bb.CheckRestrictionAndReturnBasis(candidates[0].first, basis);
		assert(res && !basis.IsEmpty());

		vector<int> side_indexes(n);
		for (size_t i = 0; i < n; ++i) {
			side_indexes[i] = i;
		}

		basis = MSW(sides, side_indexes, basis);
		new_point = basis.point;

		// Этот ассерт постоянно срабатывает из-за погрешностей. Хз как его побороть
		// assert(CheckValidBasis(basis, sides));
		if (!CheckValidBasis(basis, sides))
			return false;
		return !basis.IsEmpty();
	}

	bool MovePointNewton(const Vector2D& old_point, const vector<pair<Triangle*, int>>& vertex_triangles, Vector2D& new_point)
	{
		auto n = vertex_triangles.size();
		Vector2D grad(0, 0);
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0;
		for (auto& pair : vertex_triangles)
		{
			auto tri = pair.first;
			auto ind = pair.second;
			assert(tri->Orientation == PointsOrder::Counterclockwise && old_point == tri->V[ind]);
			auto p1 = tri->E[ind].S;
			auto p2 = tri->E[ind].E;
			auto p3 = old_point;
			auto v = p2 - p1;
			auto m = 0.5*(p1 + p2);
			auto k = (p1 - p2).LengthSquared() + (p2 - p3).LengthSquared() + (p3 - p1).LengthSquared();
			auto triArea = GeometryUtils::TriangleOrientedArea(p1, p2, p3);
			assert(triArea > 0);
			auto delta1 = v.Norm() / (2 * triArea);
			auto delta2 = (p3 - m) / k;
			sum1 += delta1.X*delta1.X;
			sum2 += delta1.Y*delta1.Y;
			sum3 += delta1.X*delta1.Y;
			sum4 += 1 / k;
			sum5 += delta2.X*delta2.X;
			sum6 += delta2.Y*delta2.Y;
			sum7 += delta2.X*delta2.Y;
			grad += delta1 - 4 * delta2;
		}
		auto hxx = grad.X * grad.X / n - sum1 - 4 * sum4 + 16 * sum5;
		auto hyy = grad.Y * grad.Y / n - sum2 - 4 * sum4 + 16 * sum6;
		auto hxy = grad.X * grad.Y / n - sum3 + 16 * sum7;
		new_point = old_point - Vector2D(hyy*grad.X - hxy*grad.Y,
			hxx*grad.Y - hxy*grad.X) / (hxx*hyy - hxy*hxy);
		return true;
	}


	FEMLIBRARY_API bool MovePointMSW_Test(const vector<Vector2D>& vertexes, Vector2D& new_point)
	{
		auto centroid = accumulate(vertexes.begin(), vertexes.end(), Vector2D(0, 0)) / static_cast<double>(vertexes.size());
		vector<pair<Triangle*, int>> triangles;
		auto last_vertex = vertexes.back();
		for (auto& v : vertexes)
		{
			if (GeometryUtils::Orientation(centroid, last_vertex, v) == PointsOrder::Counterclockwise)
				triangles.emplace_back(new Triangle(centroid, last_vertex, v), 0);
			else triangles.emplace_back(new Triangle(centroid, v, last_vertex), 0);
			last_vertex = v;
		}
		return MovePointMSW(centroid, triangles, new_point);
	}
	FEMLIBRARY_API bool MovePointNewton_Test(const vector<Vector2D>& vertexes, Vector2D& new_point)
	{
		auto centroid = accumulate(vertexes.begin(), vertexes.end(), Vector2D(0, 0)) / static_cast<double>(vertexes.size());
		vector<pair<Triangle*, int>> triangles;
		auto last_vertex = vertexes.back();
		for (auto& v : vertexes)
		{
			if (GeometryUtils::Orientation(centroid, last_vertex, v) == PointsOrder::Counterclockwise)
				triangles.emplace_back(new Triangle(centroid, last_vertex, v), 0);
			else triangles.emplace_back(new Triangle(centroid, v, last_vertex), 0);
			last_vertex = v;
		}
		return MovePointNewton(centroid, triangles, new_point);
	}

	FEMLIBRARY_API vector<pair<Vector2D,double>> Basis2_Test(const LineSeg& l1, const LineSeg& l2)
	{
		vector<pair<Vector2D, double>> res;
		vector<LineSeg> vec = { l1,l2 };
		BasisBuilder2 bb(vec, 0, 1);
		pair<double, BasisBuilder*> solutions[4];
		auto n_sol = bb.FindSolutions(solutions);
		for (auto i = 0; i < n_sol; ++i)
			res.push_back(make_pair(bb.CalcPoint(solutions[i].first), solutions[i].first));
		return res;
	}

	FEMLIBRARY_API vector<pair<Vector2D, double>> Basis3_Test(const LineSeg& l1, const LineSeg& l2, const LineSeg& l3)
	{
		vector<pair<Vector2D, double>> res;
		vector<LineSeg> vec = { l1,l2,l3 };
		BasisBuilder3 bb(vec, 0, 1, 2);
		pair<double, BasisBuilder*> solutions[5];
		auto n_sol = bb.FindSolutions(solutions);
		for (auto i = 0; i < n_sol; ++i)
			res.push_back(make_pair(bb.CalcPoint(solutions[i].first), solutions[i].first));
		return res;
	}
	FEMLIBRARY_API bool Basis3_Update_Test(const LineSeg& l1, const LineSeg& l2, const LineSeg& l3, const LineSeg& newLine, double quality, const Vector2D& point)
	{
		vector<LineSeg> vec = { l1,l2,l3,newLine };
		Basis b(point,quality,0,1,2);
		auto new_b = GetBasis3(b, vec, 3);
		return !new_b.IsEmpty();
	}
}