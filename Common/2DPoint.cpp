#include "StdAfx.h"
#include "2DPoint.h"

namespace Math {

	// Class C2DPoint //

	const C2DPoint C2DPoint::Ortx = C2DPoint(1.0, 0.0);
	const C2DPoint C2DPoint::Orty = C2DPoint(0.0, 1.0);
	const C2DPoint C2DPoint::Zero = C2DPoint(0.0, 0.0);
	//---------------------------------------------------------------------------

	C2DPoint C2DPoint::operator + (const C2DPoint& p) const
	{
		return C2DPoint(x + p.x, y + p.y);
	}

	C2DPoint C2DPoint::operator - (const C2DPoint& p) const
	{
		return C2DPoint(x - p.x, y - p.y);
	}

	C2DPoint C2DPoint::operator *  (DBL alpha) const
	{
		return C2DPoint(alpha * x, alpha * y);
	}

	C2DPoint C2DPoint::operator /  (DBL alpha) const
	{
		if (alpha < EPS) {
			//не делим на 0
			return C2DPoint(x, y);
		}

		DBL ialpha = 1 / alpha;

		return C2DPoint(ialpha * x, ialpha * y);
	}
	//---------------------------------------------------------------------------

	//Умножение на скаляр
	C2DPoint operator * (DBL alpha, const C2DPoint& p)
	{
		return C2DPoint(alpha * p.x, alpha * p.y);
	}
	//---------------------------------------------------------------------------

	C2DPoint& C2DPoint::operator = (const C2DPoint& p)
	{
		x = p.x;
		y = p.y;
		return *this;
	}

	C2DPoint& C2DPoint::operator += (const C2DPoint& p)
	{
		x += p.x;
		y += p.y;
		return *this;
	}

	C2DPoint& C2DPoint::operator -= (const C2DPoint& p)
	{
		x -= p.x;
		y -= p.y;
		return *this;
	}

	C2DPoint& C2DPoint::operator *= (DBL alpha)
	{
		x *= alpha;
		y *= alpha;
		return *this;
	}

	C2DPoint& C2DPoint::operator /= (DBL alpha)
	{
		if (alpha < EPS) {
			//не делим на 0
			return *this;
		}
		
		DBL ialpha = 1 / alpha;
		x *= ialpha;
		y *= ialpha;
		return *this;
	}
	//---------------------------------------------------------------------------

	DBL C2DPoint::operator * (const C2DPoint& p) const
	{
		return x * p.x + y * p.y;
	}
	//---------------------------------------------------------------------------

	bool C2DPoint::operator == (const C2DPoint& p) const
	{
		return (fabs(x - p.x) < EPS) && (fabs(y - p.y) < EPS);	// ==
	}

	bool C2DPoint::operator != (const C2DPoint& p) const
	{
		return (fabs(x - p.x) > EPS) || (fabs(y - p.y) > EPS);	// !=
	}
	//---------------------------------------------------------------------------

	C2DPoint::operator const DBL*() const
	{
		return coord;
	}

	//Сравнение с точкой (0; 0)
	C2DPoint::operator bool() const
	{
		return this->operator != (Zero);
	}

	//Сравнение с точкой (0; 0)
	bool C2DPoint::operator ! () const
	{
		return this->operator == (Zero);
	}
	//---------------------------------------------------------------------------


	DBL C2DPoint::Norm() const
	{
		return sqrt(x * x + y * y);
	}

	DBL C2DPoint::Norm2() const
	{
		return x * x + y * y;
	}

	DBL C2DPoint::Norm_L1() const
	{
		return fabs(x) + fabs(y);
	}
	//---------------------------------------------------------------------------

	void C2DPoint::Normalize()
	{
		DBL norm = Norm();

		if (fabs(norm) > EPS)
		{
			DBL inorm = 1 / norm;

			x *= inorm;
			y *= inorm;
		}
		// иначе - нулевой вектор
	}
	//---------------------------------------------------------------------------

	C2DPoint C2DPoint::Rotate(DBL phi, bool ccw) const
	{
		DBL s, c;
		if (!ccw) {
			phi = -phi;
		}
		//_fm_sincos(phi, &s, &c);
		s = sin(phi);
		c = cos(phi);

		return C2DPoint(x * c + y * s, -x * s + y * c);
	}
	//---------------------------------------------------------------------------

	C2DPoint C2DPoint::Min(const C2DPoint& p) const // возвращает вектор с минимальными  координатами
	{
		return C2DPoint(min(x, p.x), min(y, p.y));
	}

	C2DPoint C2DPoint::Max(const C2DPoint& p) const // возвращает вектор с максимальными координатами
	{
		return C2DPoint(max(x, p.x), max(y, p.y));
	}
	//---------------------------------------------------------------------------

	C2DPoint C2DPoint::Lerp(const C2DPoint& p, DBL t) const
	{
		return C2DPoint((p.x - x) * t + x,
			(p.y - y) * t + y);
	}
	//---------------------------------------------------------------------------

	C2DPoint& C2DPoint::FromAngle(DBL angle)
	{
		x = cos(angle);
		y = sin(angle);
		return *this;
	}

	DBL C2DPoint::ToAngle() const
	{
		if (Norm() < EPS) {
			return 0.0;
		}

		DBL angle = acos(x / Norm());

		if (y < 0) {
			angle = M_2PI - angle;
		}
		return angle;
	}

	/* Return a positive value if the points p, p1, and p2 occur in counterclockwise order;
	a negative value if they occur in clockwise order; and zero if they are collinear. */
	bool C2DPoint::Right(const C2DPoint& p1, const C2DPoint& p2) const
	{
		return ((p2.x - p1.x) * (m_y - p1.y) <= (m_x - p1.x) * (p2.y - p1.y));
	}
	//---------------------------------------------------------------------------


	//Аппроксимация
	double APPR_ex(std::vector<CPoint2D> &XY, double x) {
		
		if (XY.size() <= 0) {
			CDlgShowError cError(ID_ERROR_2DPOINT_XY_LESS_0); //CDlgShowError cError(_T("XY.size() <= 0"));
			return 0.0;
		}
		
		size_t i;

		for (i = 0; i<XY.size(); i++) {
			if (XY[i].x >= x) {
				break;
			}
		}

		if (i == 0) {
			return XY[i].y;
		}
		if (i == XY.size()) {
			return XY[i - 1].y;
		}

		return APPR(XY[i - 1].x, XY[i - 1].y, XY[i].x, XY[i].y, x);
	}

	bool CircleCenter(const C2DPoint& p1, const C2DPoint& p2, DBL r, C2DPoint& c)
	{

		C2DPoint dp = p2 - p1;

		DBL discr = 4 * r * r / dp.Norm2() - 1;

		if (discr < 0)
			return false;

		double sqrt_discr = sqrt(discr);

		if (r > 0)
		{
			c.x = p1.x - 0.5 * (-dp.x + dp.y * sqrt_discr);
			c.y = p1.y - 0.5 * (-dp.y - dp.x * sqrt_discr);
			return true;
		}

		c.x = p1.x - 0.5 * (-dp.x - dp.y * sqrt_discr);
		c.y = p1.y - 0.5 * (-dp.y + dp.x * sqrt_discr);
		return true;
		//*/
	}

	//TODO: Доделать, странное пересечение
	bool CircleCenter(const C2DPoint& pt1, const C2DPoint& pt2, const C2DPoint& pt3, C2DPoint& center)
	{
		/* Проведем через пары точек две прямые.
		Первая линия (a) пусть проходит через P1 и P2, а прямая b - через P2 и P3.
		Центр круга - находится на пересечении двух перпендикулярных прямых, проходящих через середины отрезков P1P2 и P2P3.
		Вычислим уравнения прямых, перпендикулярных a и b и проходящих через середины P1P2 и P2P3.
		Они пересекаются в центре окружности.
		*/

		bool bCenter = false;
		//TODO: Находим перпендикулярные отрезки (пересечение и есть центр окружности)
		if (!IsPerpendicular(pt1, pt2, pt3))			bCenter = CalcCircleCenter(pt1, pt2, pt3, center);
		else if (!IsPerpendicular(pt1, pt3, pt2))		bCenter = CalcCircleCenter(pt1, pt3, pt2, center);
		else if (!IsPerpendicular(pt2, pt1, pt3))		bCenter = CalcCircleCenter(pt2, pt1, pt3, center);
		else if (!IsPerpendicular(pt2, pt3, pt1))		bCenter = CalcCircleCenter(pt2, pt3, pt1, center);
		else if (!IsPerpendicular(pt3, pt2, pt1))		bCenter = CalcCircleCenter(pt3, pt2, pt1, center);
		else if (!IsPerpendicular(pt3, pt1, pt2))		bCenter = CalcCircleCenter(pt3, pt1, pt2, center);
		else { return false; }

		return bCenter;
	}

	//TODO: Доделать, странное пересечение
	bool IsPerpendicular(const C2DPoint& pt1, const  C2DPoint& pt2, const C2DPoint& pt3) {

		double yDelta_a = pt2.y - pt1.y;	//отрезок 1
		double xDelta_a = pt2.x - pt1.x;

		double yDelta_b = pt3.y - pt2.y;	//отрезок 2
		double xDelta_b = pt3.x - pt2.x;

		// проверяем, не являются ли линии вертикальными
		if (fabs(xDelta_a) <= EPS && fabs(yDelta_b) <= EPS) {
			CDlgShowError cError(ID_ERROR_2DPOINT_POINTS_XY_AXIS); //_T("The points are pependicular and parallel to x-y axis\n"));
			return false;
		}

		if (fabs(yDelta_a) <= EPS) {		//TRACE(" A line of two point are perpendicular to x-axis 1\n");
			return true;
		}
		else if (fabs(yDelta_b) <= EPS) {	//TRACE(" A line of two point are perpendicular to x-axis 2\n");
			return true;
		}
		else if (fabs(xDelta_a) <= EPS) {	//TRACE(" A line of two point are perpendicular to y-axis 1\n");
			return true;
		}
		else if (fabs(xDelta_b) <= EPS) {	//TRACE(" A line of two point are perpendicular to y-axis 2\n");
			return true;
		}

		return false;
	}

	bool CalcCircleCenter(const C2DPoint& pt1, const C2DPoint& pt2, const C2DPoint& pt3, C2DPoint& center) {

		double yDelta_a = pt2.y - pt1.y;	//отрезок 1
		double xDelta_a = pt2.x - pt1.x;

		double yDelta_b = pt3.y - pt2.y;	//отрезок 2
		double xDelta_b = pt3.x - pt2.x;

		if (fabs(xDelta_a) <= EPS && fabs(yDelta_b) <= EPS) {
			center.m_x = 0.5*(pt2.x + pt3.x);
			center.m_y = 0.5*(pt1.y + pt2.y);

			return true;
		}

		// IsPerpendicular() assure that xDelta(s) are not zero
		double aSlope = yDelta_a / xDelta_a; // 
		double bSlope = yDelta_b / xDelta_b;
		if (fabs(aSlope - bSlope) <= EPS) {	// checking whether the given points are colinear. 	
			CDlgShowError cError(ID_ERROR_2DPOINT_POINTS_COLLINEAR); //_T("The three pts are collinear\n"));
			return false;
		}

		// calc center
		center.m_x = (aSlope*bSlope*(pt1.m_y - pt3.m_y) + bSlope*(pt1.m_x + pt2.m_x) - aSlope*(pt2.m_x + pt3.m_x)) / (2 * (bSlope - aSlope));
		center.m_y = (pt1.y + pt2.y) / 2 - (center.x - (pt1.x + pt2.x) / 2) / aSlope;

		return true;
	}

	C2DPoint ClosestPoint(const C2DPoint& v, const C2DPoint& p1, const C2DPoint& p2)
	{
		C2DPoint n1(p2.x - p1.x, p2.y - p1.y), n2;

		if (fabs(n1.x) > EPS)
		{
			n2.x = -n1.y / n1.x;
			n2.y = 1;
		}
		else
		{
			n2.x = 1;
			n2.y = -n1.x / n1.y;
		}

		double t = ((p1.y - v.y) * n2.x - n2.y * (p1.x - v.x)) / (n2.y * n1.x - n1.y * n2.x);

		if (t > 1) return p2;
		if (t < 0) return p1;

		return n1 * t + p1;
	}

	C2DPoint ClosestPoint(const C2DPoint& v, const C2DPoint& p1, const C2DPoint& p2, DBL r)
	{
		if (fabs(r) < EPS) {
			return ClosestPoint(v, p1, p2);
		} //*/

		C2DPoint c;
		CircleCenter(p1, p2, r, c);
		//bool bCenter = CircleCenter(p1, p2, p3, c);

		DBL a1, a2, a3;
		a1 = ToAngle(p1 - c);
		a2 = ToAngle(p2 - c);
		a3 = ToAngle(v - c);

		r = fabs(r);
		//DBL r = (c - p1).Norm();	//радиус, положительный

		if (a1 > a2) swap(a1, a2);

		if (a3 < a1) {
			return (c + FromAngle(a1) * r);
		}
		else if (a3 > a2) {
			return (c + FromAngle(a2) * r);
		}
		else {
			return (c + FromAngle(a3) * r);
		}
	}

	DBL ToAngle(const C2DPoint& p)
	{
		if (p.Norm() < EPS) return 0.0;

		DBL angle = acos(p.x / p.Norm());

		if (p.y < 0) {
			angle = M_2PI - angle;
		}
		return angle;
	}

	C2DPoint FromAngle(DBL angle)
	{
		return C2DPoint(cos(angle), sin(angle));
	}

	inline C2DPoint i01toAB(const C2DPoint& A, const C2DPoint& B, DBL x) {
		return C2DPoint(i01toAB(DBL(A.x), DBL(B.x), x),
			i01toAB(DBL(A.y), DBL(B.y), x));
	}

	//////////////////////////////////////////////////////////
	// borkhes чудил здесь
	//////////////////////////////////////////////////////////

	//! расстояние до точки
	double C2DPoint::Len(const C2DPoint &p) const
	{
		return  sqrt((p.x - x)*(p.x - x) + (p.y - y)*(p.y - y));
	}
	//---------------------------------------------------------------------------

	bool C2DPoint::Cross(const C2DPoint &pa1, const C2DPoint &pa2, const C2DPoint &pb1, const C2DPoint &pb2) {
		if (min(pa1.x, pa2.x) > max(pb1.x, pb2.x) || max(pa1.x, pa2.x) < min(pb1.x, pb2.x) ||
			min(pa1.y, pa2.y) > max(pb1.y, pb2.y) || max(pa1.y, pa2.y) < min(pb1.y, pb2.y))return false;

		C2DPoint pa = pa2 - pa1;
		C2DPoint pb = pb2 - pb1;

		double d0 = -pa.x * pb.y + pa.y * pb.x;
		if (fabs(d0) < EPS) return false;

		double t1 = (-(pb1.x - pa1.x)*pb.y + (pb1.y - pa1.y)*pb.x) / d0;
		if (t1<0.0 || t1>1.0) return false;

		double t2 = ((pb1.y - pa1.y)*pa.x - (pb1.x - pa1.x)*pa.y) / d0;
		if (t2<0.0 || t2>1.0) return false;

		C2DPoint d1 = pa1 + pa*t1,
			d2 = pb1 + pb*t2;

		(*this) = pa1 + pa*t1;
		return true;

	}
	//---------------------------------------------------------------------------

	double C2DPoint::Splitting(const C2DPoint &p1, const C2DPoint &p2) {
		if (fabs(p2.m_x - p1.m_x) > fabs(p2.m_y - p1.m_y)) {
			return (p1.m_x - m_x) / (p1.m_x - p2.m_x);
		}

		// !=
		if (fabs(p2.m_y - p1.m_y) > EPS) {
			return (p1.m_y - m_y) / (p1.m_y - p2.m_y);
		}

		return 0.0;
	}

	//---------------------------------------------------------------------------

	bool C2DPoint::Third(const C2DPoint &p1, const C2DPoint &p2) {
		C2DPoint p0, p;
		double kff = sin(M_PI / 3);		//коэффициент длины (для высоты "треугольника")

										//Средняя точка (вспомогательная)
		p0.x = (p2.x + p1.x) / 2;
		p0.y = (p2.y + p1.y) / 2;

		if (fabs(p1.x - p2.x) >= fabs(p1.y - p2.y)) {
			p.y = p0.y + (p2.x - p1.x);
			//делим не на 0
			p.x = p0.x - (p2.y - p1.y)*(p.y - p0.y) / (p2.x - p1.x);
		}
		else {
			p.x = p0.x - (p2.y - p1.y);
			//делим не на 0
			p.y = p0.y - (p2.x - p1.x)*(p.x - p0.x) / (p2.y - p1.y);
		}
		double CurrRibLen = p1.Len(p2);	//Находим длину высоты "треугольника"
		bool bRad = Rad(p0, p, CurrRibLen*kff); //Устанавливает данную точку на прямую [p0; p] на расстоянии CurrRibLen*kff

		return true;
	}
	//---------------------------------------------------------------------------

	bool C2DPoint::Rad(const C2DPoint &p1, const C2DPoint &p2, double len) {
		double kff;	//между точками или вне их отрезка?

		if (fabs(p1.x - p2.x) < EPS && fabs(p1.y - p2.y) < EPS) {
			m_x = p1.x;
			m_y = p1.y;
			return false;
		}
		kff = len / p1.Len(p2);

		m_x = p1.x + (p2.x - p1.x)*kff;
		m_y = p1.y + (p2.y - p1.y)*kff;

		return (kff <= 1) ? true : false;
	}
	//////////////////////////////////////////////////////////////////////////////
	double C2DPoint::GetDistance(const C2DPoint &p1, const C2DPoint &p2, C2DPoint *pDPRes) const {
		double res, b, k;

		//вертикальный отрезок [p1;p2]
		if (fabs(p1.x - p2.x) < EPS) {
			pDPRes->x = p1.x;
			if ((p2.y >= p1.y && y > p2.y) ||
				(p2.y < p1.y && y < p2.y)) {
				pDPRes->y = p2.y;
				return Len(p2);
			}
			else {
				if ((p1.y >= p2.y && y > p1.y) ||
					(p1.y < p2.y && y < p1.y)) {
					pDPRes->y = p1.y;
					return Len(p1);
				}
				else {
					res = x - p1.x;
					pDPRes->y = y;
					return (res >= 0) ? res : -res;

				}
			}
		}

		//горизонтальный отрезок [p1;p2]
		if (fabs(p1.y - p2.y) < EPS) {
			pDPRes->y = p1.y;
			if ((p2.x >= p1.x && x > p2.x) ||
				(p2.x < p1.x && x < p2.x)) {
				pDPRes->x = p2.x;
				return Len(p2);
			}
			else {
				if ((p1.x >= p2.x && x > p1.x) ||
					(p1.x < p2.x && x < p1.x)) {
					pDPRes->x = p1.x;
					return Len(p1);
				}
				else {
					res = y - p1.y;
					pDPRes->x = x;
					return (res >= 0) ? res : -res;
				}
			}
		}

		//делим не на 0
		k = (p1.y - p2.y) / (p1.x - p2.x);
		b = -k*p2.x + p2.y;
		pDPRes->x = (-k*b + k*y + x) / (k*k + 1);
		pDPRes->y = (b + k*k*y + k*x) / (k*k + 1);

		if ((p2.x >= p1.x && pDPRes->x > p2.x) ||
			(p2.x < p1.x && pDPRes->x < p2.x)) {
			pDPRes->y = p2.y;
			pDPRes->x = p2.x;
			return Len(p2);
		}
		else {
			if ((p1.x >= p2.x && pDPRes->x > p1.x) ||
				(p1.x < p2.x && pDPRes->x < p1.x)) {
				pDPRes->y = p1.y;
				pDPRes->x = p1.x;
				return Len(p1);
			}
			else {
				res = (b - y + k*x) / sqrt(k*k + 1);
				return (res >= 0) ? res : -res;
			}
		}
	}
	//---------------------------------------------------------------------------

	//Близко ли находятся точки, сравниваем с EPS_3
	bool C2DPoint::IfTooClose(const C2DPoint &p1, const C2DPoint &p2) {
		if (p1.Len(p2) < EPS_3) {
			return true;
		}
		return false;
	}

	int C2DPoint::GetHash(const C2DPoint& p)
	{
		auto hash = 17;		//требуется 2 натуральных числа (выбрали 17 и 23)
		hash = hash * 23 + std::hash<double>()(p.X);
		hash = hash * 23 + std::hash<double>()(p.Y);
		return hash;
	}

/*
/////////////////
**	C2DRECT
/////////////////
*/
	C2DRect::C2DRect(DBL l, DBL r, DBL b, DBL t)
	{
		left = l;
		top = t;
		right = r;
		bottom = b;
	}

	C2DRect::C2DRect(const Math::C2DPoint& p0, const Math::C2DPoint& p1)
	{
		left = min(p0.x, p1.x);
		right = max(p0.x, p1.x);
		top = max(p0.y, p1.y);
		bottom = min(p0.y, p1.y);
	}

	void C2DRect::Fix()
	{
		Math::sortUp(left, right);
		Math::sortUp(bottom, top);
	}

	void C2DRect::AddRect(const C2DRect& rect)
	{
		left = min(left, rect.left);
		right = max(right, rect.right);
		bottom = min(bottom, rect.bottom);
		top = max(top, rect.top);
	}

	void C2DRect::AddPoint(const C2DPoint& p)
	{
		left = min(left, p.x);
		right = max(right, p.x);
		top = max(top, p.y);
		bottom = min(bottom, p.y);
	}

	void C2DRect::Union(const C2DPoint& p)
	{
		MinX = min(MinX, p.X);
		MinY = min(MinY, p.Y);
		MaxX = max(MaxX, p.X);
		MaxY = max(MaxY, p.Y);
	}

	void C2DRect::SetPoint(const C2DPoint& p)
	{
		left = right = p.x;
		top = bottom = p.y;
	}

/*
/////////////////
**	C2DLINE
/////////////////
*/

	bool C2DLine::operator == (const C2DLine& rhs) const
	{
		return S == rhs.S && E == rhs.E || S == rhs.E && E == rhs.S;	// ==
	}

	bool C2DLine::operator != (const C2DLine& rhs) const
	{
		return (S != rhs.S || S != rhs.E) && (E != rhs.E || E != rhs.S);	// !=
	}

	C2DLine & C2DLine::operator=(const C2DLine & rhs)
	{
		S = rhs.S;
		E = rhs.E;
		CalcNormC();
		return *this;
	}

	void C2DLine::CalcNormC()
	{
		m_n._0 = m_p1.y - m_p2.y;
		m_n._1 = m_p2.x - m_p1.x;
		m_c = m_p1.x * m_p2.y - m_p2.x * m_p1.y;
	}

	C2DLine::C2DLine(const C2DPoint& n, DBL c)
	{
		if (fabs(n.x) < EPS && fabs(n.y) < EPS) // нормаль слишком близка к нулю
			return;

		m_n = n, m_c = c;
		CalcPoints(true);
	}

	C2DLine::C2DLine(const C2DPoint& p1, const C2DPoint& p2)
	{
		//if (fabs(p2.x - p1.x) < EPS && fabs(p1.y - p2.y) < EPS) // точки слишком близко
		//	return;
		m_p1 = p1;
		m_p2 = p2;
		CalcNormC();
	}

	C2DLine::C2DLine(double s1, double e1, double s2, double e2) {
		S = Math::C2DPoint(s1, e1);
		E = Math::C2DPoint(s2, e2);
		CalcNormC();
	}

	bool C2DLine::Inside(const C2DPoint& p)const {
		bool bX = false;
		if (fabs(S.X - E.X) > 0.0) { //! none equal (!=)
			bX = S.X <= p.X && p.X <= E.X || E.X <= p.X && p.X <= S.X;
		}
		else {
			bX = S.Y <= p.Y && p.Y <= E.Y || E.Y <= p.Y && p.Y <= S.Y;
		}
		return (C2DPoint::CrossProduct(Vector(), p - S) < EPS) && bX;
	}

	void C2DLine::CalcPoints(bool dir) {
		if (fabs(m_n._0) < EPS)
		{
			if (m_n._1 > EPS) {
				m_p1.y = m_c / m_n._1;	// равно c / b
			}
			else {	//TODO: А вдруг у нас почти горизонтальная?
				m_p1.y = 0.0;	//TODO: Magic Number
			}
			
			m_p1.x = 0.0;
			m_p2.y = m_p1.y;
			m_p2.x = 1.0;
			return;
		}
		m_p1.x = 0.0;
		m_p1.y = GetY(m_p1.x);

		if (dir) m_p2.x = 1.0;
		else     m_p2.x = -1.0;

		m_p2.y = GetY(m_p2.x);
	}

	void C2DLine::SetPoints(const C2DPoint& p1, const C2DPoint& p2)
	{
		if (p1 == p2)
			return;

		m_p1 = p1;
		m_p2 = p2;

		CalcNormC();
	}

	void C2DLine::SetNormC(const C2DPoint& n, const DBL& c, bool dir)
	{
		if (n == C2DPoint::Zero)
			return;

		m_n = n;
		m_c = c;

		CalcPoints(dir);
	}

	DBL C2DLine::GetY(DBL x)
	{
		if (fabs(m_n._1) > EPS)
			return -(m_n._0 * x + m_c) / m_n._1; // -(Ax + C) / B

		return 0;
	}

	DBL C2DLine::GetX(DBL y)
	{
		if (fabs(m_n._0) > EPS)
			return -(m_n._1 * y + m_c) / m_n._0; // -(By + C) / A

		return 0;
	}

	int C2DLine::GetIntersection(const C2DLine& l, C2DPoint& res)
	{
		C2DPoint n1 = m_n, n2 = l.m_n;
		n1.Normalize();
		n2.Normalize();
		if (n1 == n2)
		{
			if (fabs(m_c - l.m_c) < EPS)
				return -2;
			else
				return -1;
		}
		const DBL& A1 = m_n._0;
		const DBL& A2 = l.m_n._0;
		const DBL& B1 = m_n._1;
		const DBL& B2 = l.m_n._1;
		const DBL& C1 = m_c;
		const DBL& C2 = l.m_c;

		DBL d = A1 * B2 - A2 * B1; //(A1B2-A2B1)
		res.x = (B1 * C2 - B2 * C1) / d; //(B1C2-B2C1)/(A1B2-A2B1)
		res.y = (C1 * A2 - C2 * A1) / d; //(C1A2-C2A1)/(A1B2-A2B1)
		return 0;
	}

	C2DLine C2DLine::GetPerpendicular(const C2DPoint& p) const
	{
		/*
		const DBL &A = m_n.x,
		&B = m_n.y,
		&C = m_c;

		DBL S = sqrt(A * A + B * B + C * C);

		DBL A1 = - B / S,
		B1 = A / S,
		C1 = - A1 * p.x - B1 * p.y;

		ret.SetNormC(C2DPoint(A1, B1), C1);
		*/

		// матрица поворота на 90 градусов против часовой стрелки
		// |0 -1|
		// |1  0|
		// т.е.
		// x' = -y
		// y' =  x

		return C2DLine(C2DPoint(-m_n.y, m_n.x), m_n.y * p.x - m_n.x * p.y);
	}

} // namespace Math