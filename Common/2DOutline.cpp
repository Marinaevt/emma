#include "StdAfx.h"
#include "2DOutline.h"
#include "Geometry/GeometryUtils.h"
#include "Geometry/GeometryPrimitives.h"
#include <unordered_map>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace Math;

IOIMPL(C2DNode, C2DNODE)
IOIMPL(C2DCurve, C2DCURVE)
IOIMPL(C2DCircleArc, C2DCIRCARC)
IOIMPL(C2DContour, C2DCONTOUR)
IOIMPL(C2DOutline, C2DOUTLINE)

// **************************** C2DNode ****************************

void C2DNode::DrawGL(GLParam &parameter)
{
	glBegin(GL_POINTS);
	glVertex2dv(m_point());
	glEnd();
}

Math::C2DRect C2DNode::GetBoundingBox() {
	Math::C2DRect rect;

	rect.left = m_point().x - 1.0;
	rect.right = m_point().x + 1.0;
	rect.bottom = m_point().y - 1.0;
	rect.top = m_point().y + 1.0;

	return rect;
}

bool C2DNode::GetBoundingBox(Math::C2DRect &rect) {
	rect.left = m_point().x - 1.0;
	rect.right = m_point().x + 1.0;
	rect.bottom = m_point().y - 1.0;
	rect.top = m_point().y + 1.0;
	return true;
}

/**/
// **************************** C2DCurve ****************************
/**/

C2DCurve::C2DCurve(int start, int end)
{
	m_parent = nullptr;
	m_start() = start;
	m_end() = end;

	RegisterMember(&m_start);
	RegisterMember(&m_end);
}

C2DNode* C2DCurve::GetStartNode() const
{
	if (m_parent)
		return m_parent->GetNode(m_start());

	return nullptr;
}

C2DNode* C2DCurve::GetEndNode() const
{
	if (m_parent)
		return m_parent->GetNode(m_end());

	return nullptr;
}

bool C2DCurve::SetNodes(size_t start, size_t end)
{
	if (start == end) {
		return false;
	}
	//*/

	m_start() = start;
	m_end() = end;

	return true;
}

DBL C2DCurve::CalcLength()
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return 0;

	Math::C2DPoint n0 = pn0->GetPoint(),
		n1 = pn1->GetPoint();

	return (n1 - n0).Norm();
}

// l - расстояние до точки, p - найденная точка, tau - длина вектора исходных точек\узлов
bool C2DCurve::GetPoint(DBL l, Math::C2DPoint& p, Math::C2DPoint& tau)
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1) {
		return false;
	}

	Math::C2DPoint n0 = pn0->GetPoint(),
		n1 = pn1->GetPoint();

	DBL alpha = Math::iABto01(0, CalcLength(), l);	//находим пропорции
	p = n0 * (1 - alpha) + n1 * alpha;
	tau = n1 - n0;
	tau.Normalize();	//получаем нормализованный вектор (end - start)

	return true;
}

//TODO: Должна возвращать true\false и точку
Math::C2DPoint C2DCurve::GetPointHalfLenght() {
	C2DPoint tmpF, tau;
	bool bGet = GetPoint(CalcLength() / 2, tmpF, tau);	//хоть пока CURVE == LINE, но всё может поменяться
	return tmpF;
}

// записывает точку пересечения отрезка p1p2 с кривой и возвращает:
// -1 - ошибка
//  0 - более одного пересечения
//  1 - одно пересечение (возвращает точку пересечения в result)
//  2 - нет пересечений
int C2DCurve::GetIntersection(const Math::C2DPoint& p1, const Math::C2DPoint& p2, Math::C2DPoint& result)
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return -1;

	Math::C2DPoint p3 = pn0->GetPoint(),
		p4 = pn1->GetPoint();

	// Pa = p1 + u * (p2 - p1)
	// Pb = p3 + v * (p4 - p3)

	DBL denom = (p4.y - p3.y) * (p2.x - p1.x) - (p4.x - p3.x) * (p2.y - p1.y);

	DBL u = (p4.x - p3.x) * (p1.y - p3.y) - (p4.y - p3.y) * (p1.x - p3.x),
		v = (p2.x - p1.x) * (p1.y - p3.y) - (p2.y - p1.y) * (p1.x - p3.x);

	if (fabs(denom) < EPS)
	{
		if (fabs(u) < EPS)	return 0; // совпадение прямых (отрезки могут и не пересекаться)
		return 2;
	}

	u /= denom;
	v /= denom;

	if (0 <= u && u <= 1 && 0 <= v && v <= 1)
	{
		result = p1 + u * (p2 - p1);
		return 1;
	}

	return 2;
}

bool C2DCurve::GetAllIntersections(const Math::C2DLine& l, std::vector<Math::C2DPoint>& result)
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return false;

	Math::C2DPoint p3 = pn0->GetPoint(),
		p4 = pn1->GetPoint();

	Math::C2DLine line;
	line.SetPoints(p3, p4);

	Math::C2DPoint p;
	int ret = line.GetIntersection(l, p);

	if (ret != 0)
		return false;

	Math::C2DLine pl1 = line.GetPerpendicular(p3),
		pl2 = line.GetPerpendicular(p4);

	if ((pl1.IsLeft(p) && !pl2.IsLeft(p)) || (!pl1.IsLeft(p) && pl2.IsLeft(p)))
	{
		result.push_back(p);
		return true;
	}

	return false;
}


// записывает в result точку на кривой, ближайшую к указанной и возвращает:
// -1 - ошибка
//  0 - ближайшая точка является начальным узлом кривой
//  1 - ближайшая точка является конечным узлом кривой
//  2 - ближайшая точка лежит где-то на кривой
int C2DCurve::GetClosestPoint(const Math::C2DPoint& p, Math::C2DPoint& result)
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return -1;

	Math::C2DPoint p1 = pn0->GetPoint(),
		p2 = pn1->GetPoint();

	Math::C2DPoint n1(p2.x - p1.x, p2.y - p1.y), n2;

	if (fabs(n1.x) > EPS) {
		//делим не на 0
		n2.x = -n1.y / n1.x;
		n2.y = 1;
	}
	else {
		n2.x = 1;
		if (fabs(n1.y) > EPS) {
			//делим не на 0
			n2.y = -n1.x / n1.y;
		}
		else {
			n2.y = 0.0;
		}
	}

	double t = ((p1.y - p.y) * n2.x - n2.y * (p1.x - p.x)) /
		(n2.y * n1.x - n1.y * n2.x);

	if (t > 1)
	{
		result = p2;
		return 1;
	}
	if (t < 0)
	{
		result = p1;
		return 0;
	}

	result = n1 * t + p1;

	return 2;
}


Math::C2DRect C2DCurve::GetBoundingBox()
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	Math::C2DRect ret;

	if (!pn0 || !pn1)
		return CRect2D();		//возвращаем пустой (0) прямоугольник

	Math::C2DPoint n0 = pn0->GetPoint(),
		n1 = pn1->GetPoint();

	ret.left = min(n0.x, n1.x);
	ret.right = max(n0.x, n1.x);
	ret.bottom = min(n0.y, n1.y);
	ret.top = max(n0.y, n1.y);

	return ret;
}

void C2DCurve::DrawGL(GLParam &parameter)
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return;

	glBegin(GL_LINES);
	glVertex2dv(pn0->GetPoint());
	glVertex2dv(pn1->GetPoint());
	glEnd();
}

//Puzino
//Создаёт копию отрезка и возвращает указатель на него
C2DCurve* C2DCurve::Clone() {
	C2DCurve *pLine = new C2DCurve();
	pLine->Copy(this);
	return pLine;
}

//Проверка на "равенство", т.е. совпадение основных параметров
bool C2DCurve::IsEqual(const C2DCurve* c) {
	if (c->GetType() != GetType())
		return false;

	return ((GetStart() == c->GetStart()) && (GetEnd() == c->GetEnd())) ||
		((GetStart() == c->GetEnd()) && (GetEnd() == c->GetStart()));
}

//Есть ли общий узел?
bool C2DCurve::IsConnectedTo(const C2DCurve* c) {

	if (c == nullptr) {
		CDlgShowError cError(ID_ERROR_2DOUTLINE_CURVE_NULL); // _T("C2DCurve is null"));
		return false;
	}

	return (GetStart() == c->GetStart()) ||
		(GetEnd() == c->GetEnd()) ||
		(GetStart() == c->GetEnd()) ||
		(GetEnd() == c->GetStart());
}

//Сколькими узлами соединены кривые?
int C2DCurve::HowConnectedTo(const C2DCurve* c) {

	if (c == nullptr) {
		CDlgShowError cError(ID_ERROR_2DOUTLINE_CURVE_NULL);
		return 0;
	}

	// по стандарту int(true) == 1, int(false) == 0

	return static_cast<int>(GetStart() == c->GetStart()) +
		static_cast<int>(GetEnd() == c->GetEnd()) +
		static_cast<int>(GetStart() == c->GetEnd()) +
		static_cast<int>(GetEnd() == c->GetStart());
}

/*
bool C2DCurve::Copy(IO::CInterface* ptr)
{
if (!IO::CInterface::Copy(ptr))
return false;

C2DCurve* right_ptr = dynamic_cast<C2DCurve*>(ptr);

m_start() = right_ptr->m_start();
m_end() = right_ptr->m_end();

return true;
}
*/
// **************************** C2DCircleArc ****************************

C2DCircleArc::C2DCircleArc()
{
	RegisterMember(&m_rad);
	RegisterMember(&m_cw);
	SetCW(false); ///< только для отладки
}

C2DCircleArc::C2DCircleArc(int start, int end, DBL rad)
{
	m_parent = nullptr;
	m_start() = start;
	m_end() = end;
	m_rad() = rad;
	
	RegisterMember(&m_rad);
	RegisterMember(&m_cw);
	SetCW(false); ///< только для отладки
}

bool C2DCircleArc::CalcCircleData(double& phi0, double& phi1, Math::C2DPoint& center) const
{
	C2DNode *pn0 = GetStartNode(), *pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return false;

	Math::C2DPoint n0 = pn0->GetPoint(),
		n1 = pn1->GetPoint();

	//TODO: для точного вычисления дуги необходимо использовать 3 точки
	//bool bPoint = GetPoint(CalcLength()/2,n2,tau); - плохая идея, рекурсия
	// Ищем среднюю точку на дуге


	if (!Math::CircleCenter(n0, n1, m_rad(), center))
		return false;

	phi0 = (n0 - center).ToAngle();
	phi1 = (n1 - center).ToAngle();

	if (!m_cw() && phi0 > phi1)
		phi1 += M_2PI;
	else
		if (m_cw() && phi0 < phi1)
			phi1 -= M_2PI;

	return true;
}

DBL C2DCircleArc::CalcLength()
{
	Math::C2DPoint center;
	DBL phi0, phi1;

	if (!CalcCircleData(phi0, phi1, center))
		return 0;

	return fabs((phi1 - phi0) * m_rad());
}

bool C2DCircleArc::GetPoint(DBL l, Math::C2DPoint& p, Math::C2DPoint& tau)
{
	Math::C2DPoint center;
	DBL phi0, phi1;

	if (!CalcCircleData(phi0, phi1, center))
		return false;

	DBL arcL = fabs((phi1 - phi0) * m_rad());

	DBL alpha = Math::iABto01(0, arcL, l);

	DBL phi = Math::i01toAB(phi0, phi1, alpha);

	p.FromAngle(phi);
	p *= fabs(m_rad());
	p += center;

	tau.FromAngle(phi + Math::sgn(phi1 - phi0) * M_PI2);

	return true;
}

Math::C2DPoint C2DCircleArc::GetPointHalfLenght() {
	C2DPoint tmpF, tau;
	bool bGet = GetPoint(CalcLength() / 2, tmpF, tau);
	return tmpF;
}

// проверяет, является ли данная точка p1 точкой на дуге с точностью eps
// true - да (возвращает точку пересечения в result)
//
bool C2DCircleArc::IfPoint(const Math::C2DPoint& p1, DBL eps, Math::C2DPoint& result)
{
	Math::C2DPoint tmp;
	GetClosestPoint(p1, tmp);
	if ((tmp - p1).Norm() <= eps) {
		result = tmp;
		return true;
	}
	return false;
}

// записывает точку пересечения отрезка p1p2 с дугой и возвращает:
// -1 - ошибка
//  0 - более одного пересечения
//  1 - одно пересечение (возвращает точку пересечения в result)
//  2 - нет пересечений
int C2DCircleArc::GetIntersection(const Math::C2DPoint& p1, const Math::C2DPoint& p2, Math::C2DPoint& result)
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return -1;

	Math::C2DPoint tau, center,
		n0 = pn0->GetPoint(),
		n1 = pn1->GetPoint(),
		n2;	//TODO: для точного нахождения центра дуги
			//bool bCenter = GetPoint(CalcLength()/2,n2,tau);	//для точного нахождения центра дуги


	Math::CircleCenter(n0, n1, m_rad(), center);

	DBL r2 = m_rad() * m_rad();
	Math::C2DPoint p1c = p1 - center,
		p2c = p2 - center;

	bool f1 = (p1c.Norm2() < r2),
		f2 = (p2c.Norm2() < r2);

	if (f1 && f2)
		return 2;

	/////Puzino
	DBL a1 = p1.x,
		a2 = p1.y,
		b1 = (p2 - p1).x,
		b2 = (p2 - p1).y;
	DBL r = GetRadius();
	DBL A = b1*b1 + b2*b2,
		B = 2 * (b1*(a1 - center.x) + b2*(a2 - center.y)),
		C = pow(a1 - center.x, 2) + pow(a2 - center.y, 2) - r*r;

	if (A < EPS) {
		return 2;	//иначе делим на 0
	}

	DBL diskr = B*B - 4 * A*C;
	// Нет точек пересечения
	if (diskr < 0) {
		return 2;
	}

	Math::C2DPoint tmp, result2;

	// Возможно одна точка пересечения
	if (fabs(diskr) < EPS) {
		//делим не на 0
		DBL t = -B / (2 * A);
		tmp.x = a1 + b1*t;
		tmp.y = a2 + b2*t;
		if ((IfPoint(tmp, EPS, result2) == true) && (fabs(t) <= 1.0)) {
			result = result2;
			return 1;
		}
		return 2;
	}
	// Возможно две точки пересечения
	if (diskr > 0) {
		//делим не на 0
		DBL t = (-B + sqrt(diskr)) / (2 * A);
		tmp.x = a1 + b1*t;
		tmp.y = a2 + b2*t;

		if (IfPoint(tmp, EPS, result2) == false) {
			//делим не на 0
			t = (-B - sqrt(diskr)) / (2 * A);
			if (fabs(t) <= 1.0) {
				tmp.x = a1 + b1*t;
				tmp.y = a2 + b2*t;
				if (IfPoint(tmp, EPS, result2) == false) {
					return 2;		// Нет пересечений
				}
				else { // Возможно одно пересечение
					result = result2;
					return 1;	// Одно пересечение
				}
			}
			else {
				return 2; // Нет пересечений
			}

		}
		else { // Одна точка есть, должны проверить вторую точку
			//делим не на 0
			t = (-B - sqrt(diskr)) / (2 * A);
			if (fabs(t) <= 1.0) { // Возможно второе пересечение
				tmp.x = a1 + b1*t;
				tmp.y = a2 + b2*t;

				if (IfPoint(tmp, EPS, result2) == true) {
					return 0;	// Два пересечения
				}
				else {
					result = result2;
					return 1;	// Одно пересечение
				}
			}
			else {
				result = result2;
				return 1;
			}
		}
	}

	return -1;
}

// записывает 2 точки пересечения отрезка p1p2 с дугой и возвращает:
// -1 - ошибка
//  1 - два пересечения (возвращает точку пересечения в result1 и result2)
//  2 - нет двух точек пересечений
int C2DCircleArc::GetIntersection2Points(const Math::C2DPoint& p1, const Math::C2DPoint& p2, Math::C2DPoint& result1, Math::C2DPoint& result2)
{
	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return -1;

	Math::C2DPoint tau, center,
		n0 = pn0->GetPoint(),
		n1 = pn1->GetPoint(),
		n2;	//TODO: для точного нахождения центра дуги
			//bool bCenter = GetPoint(CalcLength()/2,n2,tau);	//для точного нахождения центра дуги

	Math::CircleCenter(n0, n1, m_rad(), center);

	DBL r2 = m_rad() * m_rad();
	Math::C2DPoint p1c = p1 - center,
		p2c = p2 - center;

	bool f1 = (p1c.Norm2() < r2),
		f2 = (p2c.Norm2() < r2);

	if (f1 && f2)
		return 2;

	/////Puzino
	DBL a1 = p1.x,
		a2 = p1.y,
		b1 = (p2 - p1).x,
		b2 = (p2 - p1).y;
	DBL r = GetRadius();
	DBL A = b1*b1 + b2*b2,
		B = 2 * (b1*(a1 - center.x) + b2*(a2 - center.y)),
		C = pow(a1 - center.x, 2) + pow(a2 - center.y, 2) - r*r;

	if (fabs(A) < EPS) {
		return 2;
	}

	DBL diskr = B*B - 4 * A*C;
	Math::C2DPoint tmp;

	// Возможно две точки пересечения
	if (diskr > 0) {
		//делим не на 0
		DBL t = (-B + sqrt(diskr)) / (2 * A);
		tmp.x = a1 + b1*t;
		tmp.y = a2 + b2*t;

		if (IfPoint(tmp, EPS, result1) == false) {
			return 2; // Нет 2х пересечений		
		}
		else { // Одна точка есть, должны проверить вторую точку
			//делим не на 0
			t = (-B - sqrt(diskr)) / (2 * A);
			if (fabs(t) <= 1.0) { // Возможно второе пересечение
				tmp.x = a1 + b1*t;
				tmp.y = a2 + b2*t;

				if (IfPoint(tmp, EPS, result2) == true) {
					return 1;	// Два пересечения
				}
				else {
					return 2;	// Одно пересечение
				}
			}
			else {
				return 2;
			}
		}
	}
	else { //дискрименант == 0 или < 0
		return 2;
	}

	return -1;
}

// записывает все точки пересечения прямой l с кривой в result
bool C2DCircleArc::GetAllIntersections(const Math::C2DLine& l, std::vector<Math::C2DPoint>& result)
{
	// TODO: делает тоже самое, что GetIntersection, только возвращает все точки

	C2DNode *pn0 = GetStartNode(),
		*pn1 = GetEndNode();

	if (!pn0 || !pn1)
		return false;

	Math::C2DPoint tau, center,
		n0 = pn0->GetPoint(),
		n1 = pn1->GetPoint(),
		n2;	//TODO: для точного нахождения центра дуги
			//bool bCenter = GetPoint(CalcLength()/2,n2,tau);	//для точного нахождения центра дуги

	Math::CircleCenter(n0, n1, m_rad(), center);

	Math::C2DPoint p1, p2, tmp;
	p1 = l.GetP1();
	p2 = l.GetP2();

	DBL a1 = p1.x,
		a2 = p1.y,
		b1 = p2.x - p1.x,
		b2 = p2.y - p1.y;

	DBL r = GetRadius();

	DBL A = b1 * b1 + b2 * b2,
		B = 2 * (b1 * (a1 - center.x) + b2 * (a2 - center.y)),
		C = (a1 - center.x) * (a1 - center.x) + (a2 - center.y) * (a2 - center.y) - r * r;

	if (A < EPS) {
		return false;
	}

	DBL diskr = B * B - 4 * A * C;

	// Нет точек пересечения
	if (diskr < 0) {
		return false;
	}

	// Возможно одна точка пересечения
	if (fabs(diskr) < EPS)
	{
		//не делим на 0
		DBL t = -B / (2 * A);
		tmp.x = a1 + b1 * t;
		tmp.y = a2 + b2 * t;
		if (IfPoint(tmp, EPS, tmp) == 1) {
			result.push_back(tmp);
			return true;
		}
		return false;
	}

	// Возможно две точки пересечения
	if (diskr > 0) {
		//делим не на 0
		DBL t = (-B + sqrt(diskr)) / (2 * A);
		tmp.x = a1 + b1 * t;
		tmp.y = a2 + b2 * t;

		if (IfPoint(tmp, EPS, tmp) == true)
			result.push_back(tmp);

		//делим не на 0
		t = (-B - sqrt(diskr)) / (2 * A);
		tmp.x = a1 + b1 * t;
		tmp.y = a2 + b2 * t;

		if (IfPoint(tmp, EPS, tmp) == true)
			result.push_back(tmp);
	}

	return (result.size() != 0);
}

// записывает точку пересечения дуги с дугой c1 и возвращает:
// -1 - ошибка
//  0 - более одного пересечения
//  1 - одно пересечение (возвращает точку пересечения в result)
//  2 - нет пересечений
int C2DCircleArc::GetIntersectionArcs(C2DCircleArc& c1, Math::C2DPoint& result)
{
	C2DNode *pn01 = GetStartNode(),
		*pn02 = GetEndNode(),
		*pn11 = c1.GetStartNode(),
		*pn12 = c1.GetEndNode();

	if (!pn01 || !pn02 || !pn11 || !pn12)
		return -1;

	Math::C2DPoint n01 = pn01->GetPoint(),
		n02 = pn02->GetPoint(),
		n11 = pn11->GetPoint(),
		n12 = pn12->GetPoint();

	Math::C2DPoint o1, o2;		//Центры дуг
	DBL r1 = m_rad(),
		r2 = c1.GetRadius();

	//TODO: для точного нахождения центра дуги
	Math::C2DPoint tau, n03, n13;
	//bool bCenter1 = GetPoint(CalcLength()/2,n03,tau);
	//bool bCenter2 = c1.GetPoint(CalcLength()/2,n13,tau);

	Math::CircleCenter(n01, n02, r1, o1);
	Math::CircleCenter(n11, n12, r2, o2);

	DBL L = o1.Len(o2);			//Расстояние между центрами дуг

	if (L > (r1 + r2))
		return 2;

	Math::C2DPoint tmp, result2;
	Math::C2DPoint v = o2 - o1;
	v.Normalize();
	Math::C2DPoint vT = v.Rotate(M_PI / 2, true);
	vT.Normalize();

	if (fabs(L - (r1 + r2)) < EPS) {
		v *= r1;
		if (IfPoint(o1 + v, EPS, result2) == 1) {
			result = result2;
			return 1;
		}
		else {
			return 2;
		}
	}

	if (L < (r1 + r2)) {
		DBL l1 = (r1*r1 - r2*r2 + L*L) / (2 * L),
			h = r1*r1 - l1*l1;

		if (h<0) {
			return 2; //нету
		}
		if (fabs(h) < EPS) {
			if (IfPoint(o1 + v*l1, EPS, result2) == 1) {
				result = result2;
				return 1; //одно 
			}
			else {
				return 2; //нету
			}
		}
		if (h>0) {
			//Два возможных пересечения
			if ((IfPoint(o1 + v*l1 + vT*sqrt(h), EPS, result2) == 1) && (c1.IfPoint(o2 - v*(L - l1) + vT*sqrt(h), EPS, result2) == 1)) {
				if ((IfPoint(o1 + v*l1 - vT*sqrt(h), EPS, result2) == 1) && (c1.IfPoint(o2 - v*(L - l1) - vT*sqrt(h), EPS, result2) == 1)) {
					return 0; //2 пересечения
				}
				else {
					result = result2;
					return 1; //одно
				}
			}
			else { //одно или ноль пересечений
				if ((IfPoint(o1 + v*l1 - vT*sqrt(h), EPS, result2) == 1) && (c1.IfPoint(o2 - v*(L - l1) - vT*sqrt(h), EPS, result2) == 1)) {
					result = result2;
					return 1; //одно
				}
				else {
					return 2; //нету
				}
			}
		}
	}
	return -1;
}

// записывает 2 точки пересечения дуги с дугой c1 и возвращает:
// -1 - ошибка
//  1 - два пересечения (возвращает точки пересечения в result1 и result2)
//  2 - нет двух точек пересечений
int C2DCircleArc::GetIntersectionArcs2Points(C2DCircleArc& c1, Math::C2DPoint& result1, Math::C2DPoint& result2)
{
	C2DNode *pn01 = GetStartNode(),
		*pn02 = GetEndNode(),
		*pn11 = c1.GetStartNode(),
		*pn12 = c1.GetEndNode();

	if (!pn01 || !pn02 || !pn11 || !pn12)
		return -1;

	Math::C2DPoint n01 = pn01->GetPoint(),
		n02 = pn02->GetPoint(),
		n11 = pn11->GetPoint(),
		n12 = pn12->GetPoint();

	Math::C2DPoint o1, o2;		//Центры дуг
	DBL r1 = m_rad(),
		r2 = c1.GetRadius();

	//TODO: для точного нахождения центра дуги
	Math::C2DPoint tau, n03, n13;
	//bool bCenter1 = GetPoint(CalcLength()/2,n03,tau);
	//bool bCenter2 = c1.GetPoint(CalcLength()/2,n13,tau);

	Math::CircleCenter(n01, n02, r1, o1);
	Math::CircleCenter(n11, n12, r2, o2);

	DBL L = o1.Len(o2);			//Расстояние между центрами дуг		

	Math::C2DPoint tmp, res;
	Math::C2DPoint v = o2 - o1;
	v.Normalize();
	Math::C2DPoint vT = v.Rotate(M_PI / 2, true);
	vT.Normalize();

	if (L < (r1 + r2)) {
		DBL l1 = (r1*r1 - r2*r2 + L*L) / (2 * L),
			h = r1*r1 - l1*l1;

		if (h>0) {
			//Два возможных пересечения
			if ((IfPoint(o1 + v*l1 + vT*sqrt(h), EPS, res) == 1) && (c1.IfPoint(o2 - v*(L - l1) + vT*sqrt(h), EPS, res) == 1)) {
				result1 = res;
				if ((IfPoint(o1 + v*l1 - vT*sqrt(h), EPS, res) == 1) && (c1.IfPoint(o2 - v*(L - l1) - vT*sqrt(h), EPS, res) == 1)) {
					result2 = res;
					return 1; //2 пересечения
				}
				else {
					return 2; //одно
				}
			}
			else {
				return 2;
			}
		}
		else {
			return 2;
		}
	}
	else {
		return 2;
	}
	return 2;
}

// записывает в result точку на кривой, ближайшую к указанной и возвращает:
// -1 - ошибка
//  0 - ближайшая точка является начальным узлом кривой
//  1 - ближайшая точка является конечным узлом кривой
//  2 - ближайшая точка лежит где-то на кривой
int C2DCircleArc::GetClosestPoint(const Math::C2DPoint& p, Math::C2DPoint& result)
{
	Math::C2DPoint center;
	DBL phi0, phi1;

	if (!CalcCircleData(phi0, phi1, center))
		return -1;

	DBL phi = Math::ToAngle(p - center);

	if (!m_cw() && phi0 > phi)
		phi += M_2PI;
	else
		if (m_cw() && phi0 < phi)
			phi -= M_2PI;

	bool flag = false;

	if (m_cw())
		if (phi0 > phi && phi > phi1)
		{
			result = center + FromAngle(phi) * fabs(m_rad());
			flag = true;
		}

	if (!m_cw())
		if (phi0 < phi && phi < phi1)
		{
			result = center + FromAngle(phi) * fabs(m_rad());
			flag = true;
		}

	DBL L0 = (p - GetStartNode()->GetPoint()).Norm2(),
		L1 = (p - GetEndNode()->GetPoint()).Norm2(),
		L = (p - result).Norm2();

	if (L0 < L1)
	{
		if (flag)
			if (L < L0)
				return 2;
		result = GetStartNode()->GetPoint();
		return 0;
	}
	else
	{
		if (flag)
			if (L < L1)
				return 2;
		result = GetEndNode()->GetPoint();
		return 1;

	}

	/*
	C2DNode *pn0 = GetStartNode(),
	*pn1 = GetEndNode();

	if (!pn0 || !pn1)
	return -1;

	Math::C2DPoint p1 = pn0->GetPoint(),
	p2 = pn1->GetPoint();

	if (m_rad() == 0) {
	return C2DCurve::GetClosestPoint(p, result);
	}

	Math::C2DPoint c;
	Math::CircleCenter(p1, p2, m_rad(), c);

	DBL a1, a2, a3;
	a1 = Math::Angle(p1 - c);
	a2 = Math::Angle(p2 - c);
	a3 = Math::Angle(p - c);

	if (a1 > a2) swap(a1, a2);

	if (a3 < a1) {
	result = c + Point(a1) * fabs(m_rad());
	return 0;
	}
	else if (a3 > a2) {
	result = c + Point(a2) * fabs(m_rad());
	return 1;
	}
	else {
	result = c + Point(a3) * fabs(m_rad());
	return 2;
	}
	*/
}

Math::C2DRect C2DCircleArc::GetBoundingBox()
{
	Math::C2DPoint minp, maxp;
	Math::C2DPoint pbegin, pend, pcenter, tau;	//начальная, конечная и точка на середине дуги

	pbegin = GetStartNode()->GetPoint();
	pend = GetEndNode()->GetPoint();

	// первый прямоугольник
	minp.x = min(pbegin.x, pend.x);
	minp.y = min(pbegin.y, pend.y);
	maxp.x = max(pbegin.x, pend.x);
	maxp.y = max(pbegin.y, pend.y);

	DBL r = GetRadius();
	DBL length = CalcLength();		//длина дуги
	int N = 20;						//кол-во разбиений, больше дуга - больше число

	if (length > M_PI*r / 4) N = 45;	//дуга больше 45 градусов
	if (length > M_PI*r / 2) N = 90;	//дуга больше 90 градусов
	if (length > M_PI*r) N = 180;	//дуга больше 180 градусов
	DBL len = CalcLength() / ((DBL)N);	//единичный отрезок

	for (int i = 1; i<N; i++) {
		if (!GetPoint(len*i, pcenter, tau)) {
			return Math::C2DRect(minp, maxp);
		}

		minp.x = min(minp.x, pcenter.x);
		minp.y = min(minp.y, pcenter.y);

		maxp.x = max(maxp.x, pcenter.x);
		maxp.y = max(maxp.y, pcenter.y);
	}

	return Math::C2DRect(minp, maxp);
}

void C2DCircleArc::SwapNodes() {
	C2DCurve::SwapNodes();
	m_cw = !m_cw;
}

void C2DCircleArc::DrawGL(GLParam &parameter)
{

	Math::C2DPoint point, tau;

	DBL N = 60.0;		// кол-во частей в дуге
	DBL dLen = CalcLength();
	DBL dArcPart = dLen / N;	// хорда в 6 градусов
	//int N = static_cast<int>(dLen / LENGTHUNIT);
	
	glBegin(GL_LINE_STRIP);

	glVertex2dv(GetStartNode()->GetPoint());

	for (int i = 1; i < N; i++)
	{
		GetPoint(DBL(i) * dArcPart, point, tau);
		glVertex2dv(point);
	}

	glVertex2dv(GetEndNode()->GetPoint());

	glEnd();
}

//Puzino
//Создаёт копию дуги и возвращает указатель на него
C2DCircleArc* C2DCircleArc::Clone() {
	C2DCircleArc *pArc = new C2DCircleArc();
	pArc->Copy(this);
	return pArc;
}
/*
bool C2DCircleArc::Copy(IO::CInterface* ptr)
{
if (!IO::CInterface::Copy(ptr))
return false;

C2DCircleArc* right_ptr = dynamic_cast<C2DCircleArc*>(ptr);

m_start() = right_ptr->m_start();
m_end() = right_ptr->m_end();
m_rad() = right_ptr->m_rad();
m_cw() = right_ptr->m_cw();

return true;
}
*/

bool C2DCircleArc::CheckRadius(int acc) {

	if ((GetStartNode()->GetPoint() - GetEndNode()->GetPoint()).Norm() <= fabs(m_rad()) * 2) {
		return true;
	}
	else {
		m_rad() = floor((GetStartNode()->GetPoint() - GetEndNode()->GetPoint()).Norm() * acc * 0.5 + 1) / acc;

		return false;
	}
}

bool C2DCircleArc::IsEqual(const C2DCurve* cur) {

	//Проверяем совпадение узлов
	bool bPoints = (GetStart() == cur->GetStart()) && (GetEnd() == cur->GetEnd());

	if (cur->GetType() == C2DCIRCARC)
	{
		C2DCircleArc* right_p = (C2DCircleArc*)cur;	//dynamic_cast<>() неприменим из-за const
		bool bRad = fabs(right_p->GetRadius() - GetRadius()) < EPS;
		bool bCW = right_p->GetCW() == GetCW();
		return (bPoints && bRad && bCW);
	}

	return false;
}


// **************************** C2DContour ****************************

C2DContour::C2DContour(C2DOutline* parent)
{
	m_bcache() = false;
	m_closed() = false;
	m_parent = parent;
	RegisterMember(&m_indexes);
	RegisterMember(&m_closed);
	RegisterMember(&m_bcache);
}

/*
int C2DContour::operator [] (size_t i)
{
	if (i < m_indexes().size())
		return m_indexes()[i];

	return -1;	//TODO: не должен возвращать -1
}
//*/

void C2DContour::SetParent(C2DOutline* parent)
{
	m_bcache() = false;
	m_parent = parent;
}

int C2DContour::GetCurveIndexByIndex(size_t curve_i)
{
	size_t nCurvesNum = this->GetCurveCount();
	if (m_indexes().size() < curve_i) {
		return -1;
	}

	return m_indexes()[curve_i];

	//return -1;	//TODO: возвращать не -1, что-то другое
}

int C2DContour::GetCurveIndexByCurve(C2DCurve* pCur) {
	for (size_t i = 0; i < this->GetCurveCount(); i++)
		if (pCur == GetCurve(i))
			return i;

	return -1;
}

C2DCurve* C2DContour::GetCurve(size_t index)
{
	if (!m_parent)
		return nullptr;

	//if (!CheckBoundsLI(0, m_indexes().size(), index)) 
	if (index >= m_indexes().size()) {
		return nullptr;
	}
	return m_parent->GetCurve(m_indexes()[index]);
}

void C2DContour::Clear()
{
	m_bcache() = false;
	m_indexes().clear();
	UnClose(); //m_closed() = false;
}

bool C2DContour::AddCurve(int curve_i)
{
	if (!m_parent) {
		return false;
	}

	size_t head, tail;

	m_indexes().push_back(curve_i);

	bool bCheck = CheckContour(head, tail);	//контур остался односвязным?
	if (bCheck) {
		return true;	//кривая успешно добавлена
	}
	else {
		m_indexes().erase(m_indexes().begin() + m_indexes().size() - 1);	//кривая делает контур многосвязным - удаляем её
	}
	

	return false;

	/*
	m_closed() = false;
	int nCount = m_indexes().GetCount();
	//если кривых в контуре больше 1, то возможен замкнутый контур
	if (m_indexes().GetCount() > 1){

	int nHead,nTail;
	bool bCheck = CheckContour(nHead,nTail);
	//первая, предпоследняя и добавленная (последняя) кривые
	C2DCurve* first = m_parent->GetCurve(m_indexes()[0]);
	C2DCurve* last = m_parent->GetCurve(m_indexes()[m_indexes().GetCount() - 1]);
	C2DCurve* tba = m_parent->GetCurve(curve_i);

	//проверка на замкнутость
	//SE SE SE\ES
	if ( first->GetEndNode() == tba->GetStartNode() ){
	if ((tba->GetEndNode() == last->GetStartNode()) || (tba->GetEndNode() == last->GetEndNode()))
	m_closed() = true;
	}

	//SE ES SE\ES
	if ( first->GetEndNode() == tba->GetEndNode() ){
	if ((tba->GetStartNode() == last->GetStartNode()) || (tba->GetStartNode() == last->GetEndNode()))
	m_closed() = true;
	}

	//ES SE SE\ES
	if ( first->GetStartNode() == tba->GetStartNode() ){
	if ((tba->GetEndNode() == last->GetStartNode()) || (tba->GetEndNode() == last->GetEndNode()))
	m_closed() = true;
	}

	//ES ES SE\ES
	if ( first->GetStartNode() == tba->GetEndNode() ){
	if ((tba->GetStartNode() == last->GetStartNode()) || (tba->GetStartNode() == last->GetEndNode()))
	m_closed() = true;
	}
	}
	*/
}

// возвращает true, если:
// 1. контур пустой
// 2. контур состоит из одной кривой (nHead = nTail = 0)
// 3. контур не имеет раздвоений
bool C2DContour::CheckContour(size_t &nHead, size_t &nTail)
{
	size_t nCount = m_indexes().size();

	if (nCount == 0) {
		nHead = nTail = 0;
		return true;
	}

	if (nCount == 1)
	{
		nHead = nTail = 0;
		return true;
	}

	//Формируем список узлов в контуре
	std::vector<size_t> nodes;

	for (size_t i = 0; i < nCount; i++)
	{
		C2DCurve *pCur = GetCurve(i);
		size_t CurStart = pCur->GetStart(),
			CurEnd = pCur->GetEnd();
		size_t n[2] = { CurStart, CurEnd };
		for (int k = 0; k < 2; k++)
		{
			bool flag = false;
			for (size_t j = 0; j < nodes.size(); j++)
			{
				if (nodes[j] == n[k])
				{
					flag = true;
					break;
				}
			}
			if (!flag)
				nodes.push_back(n[k]);
		}
	}

	//int nNodesCount = nodes.GetCount();
	//node_order.InsertAt(0, 0, nodes.size());
	size_t nNodesSize = nodes.size();

	std::vector<int> node_order(nNodesSize);
	std::fill(node_order.begin(), node_order.end(), 0);	

	for (size_t i = 0; i < nNodesSize; i++) {
		for (size_t j = 0; j < nCount; j++)
		{
			if (GetCurve(j)->HasNode(nodes[i]))
				node_order[i] += 1;
		}
	}

	int flag = 0; //показывает нахождение соседей у головы и хвоста
	size_t nNodeOrderSize = node_order.size();
	for (size_t i = 0; i < nNodeOrderSize; i++)
	{
		if (node_order[i] > 2){
			return false;
		}

		if ((flag == 2) && (node_order[i] == 1)){
			return false;	//раздвоение в контуре
		}

		if ((flag == 1) && (node_order[i] == 1)){
			for (size_t j = 0; j < nCount; j++)
			{
				if (GetCurve(j)->HasNode(nodes[i])){
					nTail = j;
					break;
				}
			}

			flag = 2;
		}

		if ((flag == 0) && (node_order[i] == 1)){
			
			for (size_t j = 0; j < nCount; j++){
				if (GetCurve(j)->HasNode(nodes[i])){
					nHead = j;
					break;
				}
			}
			flag = 1;
		}
	}

	UnClose(); //контур разомкнут

			   // в случае замкнутого контура
	if (flag == 0)
	{
		nHead = 0;
		nTail = 1;
		Close(); //контур замкнут
	}

	return true;	//нормальный односвязный контур (открытый или закрытый)

					/*
					int c = m_indexes().GetCount();

					if (!c)
					return true;

					if (c == 1)
					{
					nHead = nTail = 0;
					return true;
					}

					std::vector<int> con; // количество узлов кривой (0, 1 или 2), к которым присоединены другие кривые
					con.InsertAt(0, 0, c);

					for (int i = 0; i < c; i++)
					{
					C2DCurve* pCuri = GetCurve(i);	//т.к. GetCurve от контура
					for (int j = 0; j < c; j++)
					{
					if (j == i)
					continue;

					C2DCurve* pCurj = GetCurve(j);	//т.к. GetCurve от контура

					con[i] += pCuri->HowConnectedTo(pCurj);
					}
					}

					int flag = 0;	//показывает нахождение соседей у головы и хвоста
					for (int i = 0; i < c; i++)
					{
					if (con[i] > 2) return false;
					if ((flag == 2) && (con[i] == 1))
					{
					return false;	//раздвоение в контуре
					}
					if ((flag == 1) && (con[i] == 1))
					{
					nTail = i;
					flag = 2;
					}
					if ((flag == 0) && (con[i] == 1))
					{
					nHead = i;
					flag = 1;
					}
					}

					// в случае замкнутого контура
					if (flag == 0)
					{
					nHead = 0;
					nTail = 1;
					m_closed() = true;
					}

					return true;
					*/
					/*
					if (m_closed()) {
					nCurve1 = 0;
					nCurve2 = m_indexes().GetCount() - 1;
					return true;
					}

					C2DCurve *pCur1, *pCur2;
					int t = m_indexes().GetCount() - 1;
					//первая и последняя кривые
					pCur1 = GetCurve(m_indexes()[0]);
					pCur2 = GetCurve(m_indexes()[t-1]);
					if (!(pCur1->GetStartNode()->GetPoint() == pCur2->GetStartNode()->GetPoint() ||
					pCur1->GetStartNode()->GetPoint() == pCur2->GetEndNode()->GetPoint() ||
					pCur1->GetEndNode()->GetPoint() == pCur2->GetStartNode()->GetPoint() ||
					pCur1->GetEndNode()->GetPoint() == pCur2->GetEndNode()->GetPoint())){
					nCurve1 = 0;
					nCurve2 = m_indexes().GetCount() - 1;
					return true;
					}

					//все кривые между
					for(int i=0; i < m_indexes().GetCount() - 1; i++){
					pCur1 = GetCurve(m_indexes()[i]);
					pCur2 = GetCurve(m_indexes()[i+1]);

					if (pCur1->GetStartNode()->GetPoint() == pCur2->GetStartNode()->GetPoint() ||
					pCur1->GetStartNode()->GetPoint() == pCur2->GetEndNode()->GetPoint() ||
					pCur1->GetEndNode()->GetPoint() == pCur2->GetStartNode()->GetPoint() ||
					pCur1->GetEndNode()->GetPoint() == pCur2->GetEndNode()->GetPoint()){
					continue;
					}

					//кривые не имеют общей точки
					nCurve1 = i;
					nCurve2 = i+1;
					return true;
					}
					return false;
					*/
}

bool C2DContour::OrderContour()
{
	size_t head, tail;

	if (!CheckContour(head, tail))
		return false;

	size_t c = m_indexes().size();

	if (!c)
		return true;

	std::vector<size_t> new_indexes;
	new_indexes.push_back(head);

	for (size_t i = 0; i < c; i++)
	{
		bool flag = false;
		for (size_t j = 0; j < new_indexes.size(); j++)
		{
			if (new_indexes[j] == m_indexes()[i]) {		//не было скобок
				flag = true;
				break;
			}
		}
		if (flag)
			continue;

		if (GetCurve(i)->IsConnectedTo(m_parent->GetCurve(new_indexes[new_indexes.size() - 1])))
			new_indexes.push_back(m_indexes()[i]);
	}
	m_indexes() = new_indexes;

	return true;
}

bool C2DContour::DelCurve(size_t curve_i)
{
	size_t nIndCount = m_indexes().size();

	if (nIndCount == 0) {
		CDlgShowError cError(ID_ERROR_2DOUTLINE_NO_CURVES_IN_CONTOUR); //_T("В контуре нет кривых \n"));
		return false;
	}

	if (!GetCurve(curve_i)) {
		CDlgShowError cError(ID_ERROR_2DOUTLINE_NO_SUCH_CURVE); //_T("Нет такой кривой \n"));
		return false;
	}

	size_t nFirst, nLast;
	bool bCurves = CheckContour(nFirst, nLast);	//проверка на односвязность контура
	if (!bCurves) return false;

	// если контур замкнут
	if (IsClosed())
	{
		for (size_t i = 0; i < nIndCount; i++) {
			if (m_indexes()[i] == curve_i)
			{
				m_indexes().erase(m_indexes().begin() + i);
				UnClose();
				return true;
			}
		}
		return false;
	}

	// если контур разомкнут
	if (curve_i == m_indexes()[nFirst])
	{
		m_indexes().erase(m_indexes().begin() + nFirst);
		return true;
	}

	if (curve_i == m_indexes()[nLast])
	{
		m_indexes().erase(m_indexes().begin() + nLast);
		return true;
	}

	//случай, когда нужно удалить кривую из середины незамкнутого контура
	for (size_t i = 0; i < nIndCount; i++){
		if (m_indexes()[i] == curve_i)
		{
			
			size_t nCount1 = m_parent->GetContourCount();

			//Создаём новый контур, т.к. предыдущий распадается на 2 части
			C2DContour* c = new C2DContour(m_parent);
			m_parent->AddContour(c);

			size_t nCount2 = m_parent->GetContourCount();

			//Добавляем в новый контур оставшиеся подвешенные в пространстве кривые
			bool bAdd = false;
			for (size_t j = i + 1; j < nIndCount; j++) {
				bAdd = c->AddCurve(m_indexes()[j]);
			}

			//m_indexes().RemoveAt(i, m_indexes().size() - i);
			m_indexes().erase(m_indexes().begin() + i, m_indexes().begin() + (nIndCount));
			return true;
		}
	}

	return false;
}

void C2DContour::DrawGL(GLParam &parameter)
{
	// Если есть вспомогательный кэш
	if (IsCache()) {

		//узлы контура
		glPointSize(5);
		glBegin(GL_POINTS);

		glColor3d(1, 0, 0);
		for (size_t i = 0; i<GetCache().size(); i++) {

			Math::C2DPoint DP1 = GetCache()[i];
			glVertex2dv(DP1);
		}
		glEnd();

		//кривые из контура
		glPointSize(2);
		glBegin(GL_LINE_STRIP);

		glColor3d(0, 0, 1);
		for (size_t i = 0; i<GetCache().size(); i++) {

			Math::C2DPoint DP1 = GetCache()[i];
			glVertex2dv(DP1);
		}
		glEnd();

	}

}

//! Заполнение кэша контура 
void C2DContour::FillCache()
{
	m_cache.clear();

	C2DCurve* pcurve;
	/*
	C2DCurve* pcurve = GetCurve(0);
	if (pcurve->GetType()==C2DCIRCARC && pcurve->GetEnd() != GetCurve(1)->GetStart()){
	//pcurve->SwapNodes();
	}//*/

	size_t nCount = m_indexes().size() - 1;

	for (size_t i = 0; i < m_indexes().size(); i++)
	{
		bool flag = (IsClosed() && (i == nCount));
		pcurve = GetCurve(i);

		//если отрезок
		if (pcurve->GetType() == C2DLINESEG)
		{
			DBL length = pcurve->CalcLength();

			Math::C2DPoint tmpP, tmpTau;

			pcurve->GetPoint(0, tmpP, tmpTau);
			m_cache.push_back(tmpP);
			pcurve->GetPoint(length, tmpP, tmpTau);
			m_cache.push_back(tmpP);

			continue;
		}

		//если дуга
		if (pcurve->GetType() == C2DCIRCARC) {

			C2DCircleArc* pcircle = dynamic_cast<C2DCircleArc*>(pcurve);
			if (!pcircle) return;

			// делим на отрезки и сохраняем точки
			DBL length = pcircle->CalcLength();
			int N = static_cast<int>(length / LENGTHUNIT);

			Math::C2DPoint tmpP, tmpT;
			for (int j = 0; j <= N; j++)
			{
				if (j == N) {
					if (tmpP.Len(pcircle->GetEndNode()->GetPoint()) < tmpP.Len(pcircle->GetStartNode()->GetPoint())) {
						m_cache.push_back(pcircle->GetEndNode()->GetPoint());
						continue;
					}
					else {
						m_cache.push_back(pcircle->GetStartNode()->GetPoint());
						continue;
					}
				}
				pcircle->GetPoint(j * LENGTHUNIT, tmpP, tmpT);
				m_cache.push_back(tmpP);
			}

			/*
			if (!flag){
			//m_cache.Add(((C2DCircleArc*)pcurve)->GetEndNode()->GetPoint());
			} //*/
		}
	}
	m_bcache() = true;
}

//! Находим среднюю точку для всего контура
Math::C2DPoint C2DContour::GetCenterPoint() {

	Math::C2DPoint point(0, 0);
	size_t count = m_indexes().size();
	for (size_t i = 0; i < count; i++) {
		point.x += m_parent->GetNode(m_indexes()[i])->GetPoint().x;
		point.y += m_parent->GetNode(m_indexes()[i])->GetPoint().y;
	}

	if (count>0) {
		point.x /= count;
		point.y /= count;
	}

	return point;
}

//Находит сумму длин всех кривых в контуре
DBL C2DContour::CalcLength() {

	DBL dLen = 0.0;
	size_t nCount = GetCurveCount();

	for (size_t i = 0; i < nCount; i++) {
		dLen += GetCurve(i)->CalcLength();
	}

	return dLen;
}


// **************************** C2DOutline ****************************

C2DOutline::C2DOutline()
{
	RegisterMember(&m_nodes);
	RegisterMember(&m_curves);
	RegisterMember(&m_contours);
}

void C2DOutline::SetParentToAll()
{
	for (size_t i = 0; i < m_curves().size(); i++)
		dynamic_cast<C2DCurve*>(m_curves()[i])->SetParent(this);

	for (size_t i = 0; i < m_contours().size(); i++)
		dynamic_cast<C2DContour*>(m_contours()[i])->SetParent(this);
}

int  C2DOutline::CreateNode(TYPEID type, const C2DPoint& p)
{
	C2DNode* pNode = dynamic_cast<C2DNode*>(IO::CreateObject(type));
	if (!pNode)
		return -1;
	pNode->SetPoint(p);
	return AddNode(pNode);
}

int  C2DOutline::AddNode(C2DNode* pNode)
{
	if (!pNode)
		return -1;

	// проверяем на недопустимую близость к другому узлу
	for (size_t i = 0; i < m_nodes().size(); i++)
	{
		C2DNode* pCurNode = dynamic_cast<C2DNode*>(m_nodes()[i]);
		if ((pNode->GetPoint() - pCurNode->GetPoint()).Norm() < EPS)
			return i; //TODO: [size] // если узел близок, то не добавляем новый, а возвращаем индекс существующего
	}

	m_nodes().push_back(pNode);
	return (m_nodes().size() - 1);
}

bool C2DOutline::DelNode(size_t index)
{
	if (index >= m_nodes().size())
		return false;

	// удаление узла из массива узлов
	if (m_nodes()[index])
		delete m_nodes()[index];

	m_nodes().erase(m_nodes().begin() + index);

	/*
	// перемещение всех узлов с индексом, большим чем index, на одну позицию вверх
	for (int i = index; i < m_nodes().GetCount() - 1; i++)
	m_nodes()[i] = m_nodes()[i + 1];

	// уменьшение размера массива на 1
	m_nodes().SetSize(m_nodes().GetCount() - 1);
	*/

	for (size_t i = 0; i < m_curves().size(); i++)
	{
		C2DCurve* p = dynamic_cast<C2DCurve*>(m_curves()[i]);

		// удаление кривых, содержащих этот узел из массива кривых
		// удаление индексов этих кривых из контуров
		if ((p->GetStart() == index) || (p->GetEnd() == index))
		{
			DelCurve(i);
			i -= 1;
		}
	}

	for (size_t i = 0; i < m_curves().size(); i++)
	{
		C2DCurve* p = dynamic_cast<C2DCurve*>(m_curves()[i]);

		// уменьшение номеров узлов, больших чем index, на 1
		if (p->GetStart() > index)
			p->SetNodes(p->GetStart() - 1, p->GetEnd());

		if (p->GetEnd() > index)
			p->SetNodes(p->GetStart(), p->GetEnd() - 1);
	}

	return true;
}

C2DNode* C2DOutline::GetNode(size_t index)
{
	if (index >= m_nodes().size())
		return nullptr;

	return dynamic_cast<C2DNode*>(m_nodes()[index]);
}

size_t C2DOutline::GetNodeCount() const
{
	return m_nodes().size();
}

//// Roma
void C2DOutline::DelRepNodes(double acc) {

	for (size_t i = 0; i < m_nodes().size() - 1; i++)
	{
		for (size_t j = i + 1; j < m_nodes().size(); j++)
		{
			if ((GetNode(i)->GetPoint() - GetNode(j)->GetPoint()).Norm() < acc)
			{
				for (size_t k = 0; k < m_curves().size(); k++)
				{
					C2DCurve* p = dynamic_cast<C2DCurve*>(m_curves()[k]);

					if (p->GetStart() == j) p->SetNodes(i, p->GetEnd());
					if (p->GetEnd() == j) p->SetNodes(p->GetStart(), i);

				}

				DelNode(j);
			}
		}
	}
}

void C2DOutline::DelRepNodes2(double acc)
{
	std::vector<size_t> tbd; // to be deleted

	for (size_t i = 0; i < m_nodes().size() - 1; i++)
	{
		for (size_t j = i + 1; j < m_nodes().size(); j++)
		{
			if ((GetNode(i)->GetPoint() - GetNode(j)->GetPoint()).Norm() < acc)
			{
				for (size_t k = 0; k < m_curves().size(); k++)
				{
					C2DCurve* p = dynamic_cast<C2DCurve*>(m_curves()[k]);

					if (p->GetStart() == j) p->SetNodes(i, p->GetEnd());
					if (p->GetEnd() == j) p->SetNodes(p->GetStart(), i);

				}
				tbd.push_back(j);
			}
		}
	}

	for (size_t i = 0; i < tbd.size(); i++) {
		DelNode(i);
	}

}

size_t C2DOutline::FindNumNode(C2DNode *pNode) {

	for (size_t i = 0; i < m_nodes().size(); i++) {
		if (pNode == GetNode(i)) return i;
	}

	return -1;
}

int  C2DOutline::CreateCurve(TYPEID type, int start, int end)
{
	C2DCurve* pCurve = dynamic_cast<C2DCurve*>(IO::CreateObject(type));
	if (!pCurve)
		return -1;
	pCurve->SetParent(this);
	bool bSet = pCurve->SetNodes(start, end);	//false if identical nodes ID
	//TODO:
	if (!bSet) {
		return -1;	//if nodes ID are identical
	}
	return AddCurve(pCurve);
}

ptrdiff_t C2DOutline::AddCurve(C2DCurve* pCurve)
{
	if (!pCurve)
	{
		CDlgShowError cError(ID_ERROR_2DOUTLINE_CURVE_NULL);
		return -1;
	}

	size_t start = pCurve->GetStart(),
		end = pCurve->GetEnd(),
		count = GetNodeCount();

	if(start > count || count < end)
	{
		CDlgShowError cError(ID_ERROR_2DOUTLINE_NODES_INDEX_WRONG); //_T("Неверные номера узлов: start, end, count \n"));
		return -1;
	}

	if (start == end) {
		CDlgShowError cError(ID_ERROR_2DOUTLINE_NODES_INDEX_WRONG); //_T("Неверные номера узлов: start, end, count \n"));
		return -1;
	}
	// проверка на совпадение с существующими кривыми

	/*for (int i = 0; i < m_curves().GetCount(); i++)
	{
	if (m_curves()[i]->GetType() == C2DCURVE){
	C2DCurve* pCurvei = (C2DCurve*)m_curves()[i];
	if (pCurvei->IsEqual(pCurve))
	return i; // если совпадает, то не добавлять кривую, а вернуть индекс существующей
	}
	if (m_curves()[i]->GetType() == C2DCIRCARC){
	C2DCircleArc* pCurvei = (C2DCircleArc*)m_curves()[i];
	if (pCurvei->IsEqual(pCurve))
	return i; // если совпадает, то не добавлять кривую, а вернуть индекс существующей
	}
	}
	*/
	pCurve->SetParent(this);

	m_curves().push_back(pCurve);
	//TODO: надо уточнять номер контура
	//AddCurveToContour(m_curves().GetCount() - 1,0);
	return (m_curves().size() - 1);
}

bool C2DOutline::DelCurve(size_t index)
{
	if (index >= m_curves().size())
		return false;

	//
	bool ret = DelCurveFromContours(index);

	// возвращаем память, занимаемую кривой
	if (m_curves()[index])
		delete m_curves()[index];

	//удаляем кривую из массива
	m_curves().erase(m_curves().begin() + index);

		// перемещение всех кривых с индексом, большим чем index, на одну позицию вверх
		//for (int i = index; i < m_curves().GetCount() - 1; i++)
		//	m_curves()[i] = m_curves()[i + 1];

		// уменьшение размера массива на 1
		//m_curves().SetSize(m_curves().GetCount() - 1);

	for (size_t i = 0; i < m_contours().size(); i++)
	{
		C2DContour* p = GetContour(i);  //(C2DContour*) m_contours()[i];

		// уменьшение на 1 номеров кривых (в контуре i), больших чем index
		for (size_t j = 0; j < p->GetCurveCount(); j++){
			if (p->At(j) > index){
				p->At(j) -= 1;
			}
		}
		//*/
	}//*/

	return ret;
}

C2DCurve* C2DOutline::GetCurve(size_t index)
{
	if (index >= m_curves().size()) {
		return nullptr;
	}
	return dynamic_cast<C2DCurve*>(m_curves()[index]);
}

/*
int C2DOutline::GetCurve(C2DCurve* pCurve) {
	for (int i = 0; i<m_curves().GetCount(); i++) {
		if (pCurve == dynamic_cast<C2DCurve*>(m_curves()[i])) {
			return i;
		}
	}
	return -1;
}
//*/
int C2DOutline::GetCurveIndexByCurve(C2DCurve* pCurve) {
	for (size_t i = 0; i<m_curves().size(); i++) {
		if (pCurve == dynamic_cast<C2DCurve*>(m_curves()[i])) {
			return i;
		}
	}
	return -1;
}

size_t C2DOutline::GetCurveCount() const
{
	return m_curves().size();
}

// Эта функция должна удалить "Висящие" узлы и линии, 
// найти внешний контур и перенумеровать его по часовой стрелке
bool C2DOutline::PrerareToMeshing() {

	C2DCurve *pCur;
	DBL rad;	//радиус для дуг

				//TODO: Баг - Виснет, если обход с самого начала против часовой.

				//Это просто заплатка.. предполагаем, что первая кривая пронумерована правильно.
	for (size_t i = 0; i < GetContourCount(); i++) {
		C2DContour *pCont = GetContour(i);

		//проверяем контур на висячие узлы и прочее
		size_t nHead, nTail;
		bool bIsCorrect = pCont->CheckContour(nHead, nTail);	//TODO: что-то надо делать с bIsCorrect
		//...//

		//перенумеровываем контуры по часовой стрелке
		Math::C2DPoint CenterFigPoint = pCont->GetCenterPoint();

		//перебираем кривые в контуре
		for (size_t j = 1; j < pCont->GetCurveCount(); j++) {

			Math::C2DPoint DP01 = pCont->GetCurve(j - 1)->GetEndNode()->GetPoint();  //конец 1-й кривой
			pCur = pCont->GetCurve(j);
			Math::C2DPoint DP10 = pCur->GetStartNode()->GetPoint(); //начало 2-й кривой
			Math::C2DPoint DP11 = pCur->GetEndNode()->GetPoint();  //конец 2-й кривой

																   //относительно центра смотрим, где могут быть узлы
			if (DP01.Len(DP10) > DP01.Len(DP11)) {

				//Дуги меняются после проверки
				if (pCur->GetType() == C2DCIRCARC) {
					C2DPoint pCenter_before, pCenter_after;
					C2DCircleArc *pArc = dynamic_cast<C2DCircleArc *>(pCur);
					pCenter_before = pArc->GetPointHalfLenght();

					//проверяем центры дуг до и после смены узлов

					pCur->SwapNodes();
					pCenter_after = pArc->GetPointHalfLenght();

					if (pCenter_after.Len(pCenter_before) > EPS) {
						rad = pArc->GetRadius();
						pArc->SetRadius(-rad);
					}

				}
				else {
					//иначе просто меняем местами узлы
					pCur->SwapNodes();
				}
			}
		}
	}

	return true;
}

int C2DOutline::AddNodeToCurve(int curve_i, DBL length)
{
	// добавление нового узла
	// создание 2х новых кривых из первого узла старой кривой в новый узел, и из нового узла во второй узел старой кривой
	// удаление старой кривой
	// удаление индексов этой кривой из контуров
	// проверка контуров на корректность

	C2DCurve* p_curve = GetCurve(curve_i);

	if (!p_curve)
		return -1;

	C2DNode* p_node = p_curve->GetStartNode();

	if (!p_node)
		return -1;

	C2DNode*  p_new_node = dynamic_cast<C2DNode*>  (IO::CreateObject(p_node->GetType()));
	C2DCurve* p_curve_1 = dynamic_cast<C2DCurve*> (IO::CreateObject(p_curve->GetType()));
	C2DCurve* p_curve_2 = dynamic_cast<C2DCurve*> (IO::CreateObject(p_curve->GetType()));

	if (!p_new_node || !p_curve_1 || !p_curve_2)
		return -1;

	Math::C2DPoint point, tau;
	p_curve->GetPoint(length, point, tau);
	p_new_node->SetPoint(point);

	int i_new_node = AddNode(p_new_node);

	if (i_new_node == -1)
		return -1;

	p_curve_1->SetNodes(p_curve->GetStart(), i_new_node);
	p_curve_2->SetNodes(i_new_node, p_curve->GetEnd());

	if (!DelCurve(curve_i))
		return -1;

	if (AddCurve(p_curve_1) == -1)
		return -1;

	if (AddCurve(p_curve_2) == -1)
		return -1;

	return i_new_node;
}

int  C2DOutline::CreateContour(TYPEID type)
{
	C2DContour* pContour = dynamic_cast<C2DContour*>(IO::CreateObject(type));
	if (!pContour)
		return -1;
	pContour->SetParent(this);
	return AddContour(pContour);
}

size_t C2DOutline::AddContour(C2DContour* pContour)
{
	if (!pContour)
		return -1;

	pContour->SetParent(this);

	m_contours().push_back(pContour);
	return (m_contours().size() - 1);
}

bool C2DOutline::DelContour(size_t index)
{
	if (index >= m_contours().size()) {
		return false;
	}

	//освобождаем память
	//if (m_contours()[index])
	delete m_contours()[index];

	//удаление и сдвиг номеров кривых в контуре
	m_contours().erase(m_contours().begin() + index);

	/*for (int i = index; i < m_contours().GetCount() - 1; i++)
	m_contours()[i] = m_contours()[i + 1];

	m_contours().SetSize(m_contours().GetCount() - 1); //*/

	return true;
}

C2DContour* C2DOutline::GetContour(size_t index)
{
	if (index >= m_contours().size()) {
		return nullptr;
	}
	return dynamic_cast<C2DContour*>(m_contours()[index]);
}

size_t C2DOutline::GetContourCount() const
{
	return m_contours().size();
}

bool C2DOutline::FindContours(double ribLen)
{
	//Две точки могут образовать контур (две дуги из двух точек с разными радиусами), меньше - нет
	if (GetNodeCount() < 2)
		return false;
	// заполняем пул непересекающихся отрезков
	FemLibrary::NonIntersectedCurvesPool pool;
	for (size_t p = 0; p < GetCurveCount(); p++)
	{
		C2DCurve *pCur = GetCurve(p);
		Math::C2DPoint DP = pCur->GetStartNode()->GetPoint(), DPa;

		auto lastPoint = pCur->GetStartNode()->GetPoint();

		// дуги разбиваются на несколько отрезков таким же образом, как и в CRibLine::Load
		if (pCur->GetType() != C2DLINESEG) {
			auto totalLen = pCur->CalcLength();
			for (DBL l1 = ribLen; l1 < totalLen; l1 += ribLen)
			{
				if (l1 > totalLen - ribLen)
					l1 = 0.5*(l1 + totalLen - ribLen);
				pCur->GetPoint(l1, DP, DPa);
				pool.AddLine(FemLibrary::LineSeg(lastPoint.x, lastPoint.y, DP.x, DP.y));

				lastPoint = DP;
			}
		}
		DP = pCur->GetEndNode()->GetPoint();
		pool.AddLine(FemLibrary::LineSeg(lastPoint.x, lastPoint.y, DP.x, DP.y));
	}

	// ищем контуры в массиве непересекеающихся отрезков
	auto contours = FemLibrary::GeometryUtils::FindContours(pool.Lines);

	// заново заполняем Outline 
	ClearAll();

	// хеш-таблицы уже добавленных в Outline точек и отрезков для того, чтобы не плодить дубликаты 
	std::unordered_map<FemLibrary::Vector2D, int> pointIndexes;
	std::unordered_map<FemLibrary::LineSeg, int> lineIndexes;
	for (auto contour : contours)
	{
		// создаем контур в Outline
		auto contourId = CreateContour(C2DCONTOUR);

		// сначала добавляем самую последнюю точку контура — это предыдущая точка
		auto prevPoint = contour.ContourPoints.back();
		// для каждой точки и отрезка мы сначала ищем их в соответствующей таблице
		// если примитив найден, то извлекаем из таблицы его индекс в Outline
		// если нет — то добавляем его в Outline и хеш-таблицу
		int prevInd;
		auto pointIndexesIt = pointIndexes.find(prevPoint);
		if (pointIndexesIt == pointIndexes.end())
			pointIndexes[prevPoint] = prevInd = CreateNode(C2DNODE, Math::C2DPoint(prevPoint.X, prevPoint.Y));
		else
			prevInd = pointIndexesIt->second;

		for (auto point : contour.ContourPoints)
		{
			// добавляем все точки по порядку
			int ind;
			pointIndexesIt = pointIndexes.find(point);
			if (pointIndexesIt == pointIndexes.end())
				pointIndexes[point] = ind = CreateNode(C2DNODE, Math::C2DPoint(point.X, point.Y));
			else
				ind = pointIndexesIt->second;

			// предыдущая и текущая точка формируют отрезок контура
			int lineInd;
			FemLibrary::LineSeg lineSeg(prevPoint, point);
			auto lineIndexesIt = lineIndexes.find(lineSeg);
			if (lineIndexesIt == lineIndexes.end())
				lineIndexes[lineSeg] = lineInd = CreateCurve(C2DLINESEG, prevInd, ind);
			else
				lineInd = lineIndexesIt->second;

			// добавляем индекс отрезка в контур Outline
			GetContour(contourId)->AddCurve(lineInd);

			// текущая точка становится предыдущей
			prevInd = ind;
			prevPoint = point;
		}
	}
	return true;
}

bool C2DOutline::AddCurveToContour(int curve_i, int contour_i)
{
	C2DContour* pCon = GetContour(contour_i);
	if (!pCon)
		return false;
	if (pCon->IsClosed()) {
		return false;	// если замкнут, то добавлять ничего не можем //con->UnClose();
	}

	C2DCurve* pCur = GetCurve(curve_i);
	if (!pCur)
		return false;

	return pCon->AddCurve(curve_i);
}

bool C2DOutline::DelCurveFromContour(size_t curve_i, size_t contour_i)
{
	if (contour_i >= m_contours().size()) {
		return false;
	}
	C2DContour* p = GetContour(contour_i);

	// удаление кривой из контура
	bool bDel = p->DelCurve(curve_i);
	return bDel;
}

bool C2DOutline::DelCurveFromContours(size_t curve_i)
{
	bool bRet = true;
	size_t i = GetContourCount(),
		nCount = 0;		//счётчик удалений

	while (i != 0) { //Обход справа налево принципиален из за удаления пустых контуров
		
		i--;	//Обход справа налево принципиален из за удаления пустых контуров
		
		if (!DelCurveFromContour(curve_i, i)) {
			bRet = false;
		}else {
			//удаляем пустой контур
			if (GetContour(i)->GetCurveCount() == 0) {
				DelContour(i);
			}
			nCount++;
		}
	}

	if (nCount > 0) bRet = true;	//кривая удалилась хотя бы из 1 контура
	return bRet;
}

bool C2DOutline::GetPoint(int c, DBL l, Math::C2DPoint& p, Math::C2DPoint& tau)
{
	DBL tmp = 0, cl;

	if (l < 0) {
		return false;
	}

	//C2DContour& pC = *dynamic_cast<C2DContour*>(m_contours()[c]);
	C2DContour* pC = GetContour(c);
	if (pC == nullptr) {
		return false;
	}

	for (size_t i = 0; i < pC->GetCurveCount()/* - 1*/; i++)
	{
		//int curveNum = pC->GetCurve(i);

		C2DCurve* pCurve = pC->GetCurve(i); //dynamic_cast<C2DCurve*>(m_curves()[curveNum]);

		cl = pCurve->CalcLength();
		tmp += cl;
		if (l < tmp) {
			return pCurve->GetPoint(l - (tmp - cl), p, tau);
		}
	}

	return false;
}

Math::C2DPoint C2DOutline::GetCenterPoint() {

	Math::C2DPoint point(0, 0);
	size_t count = GetNodeCount();
	for (size_t i = 0; i<count; i++) {
		point.x += GetNode(i)->GetPoint().x;
		point.y += GetNode(i)->GetPoint().y;
	}
	if (count>0) {
		point.x /= count;
		point.y /= count;
	}
	return point;
}

// записывает точку пересечения отрезка p1p2 с кривой curve_num и возвращает:
// -1 - ошибка
//  0 - более одного пересечения
//  1 - одно пересечение (возвращает точку пересечения в result)
//  2 - нет пересечений
int C2DOutline::GetIntersection(const Math::C2DPoint& p1, const Math::C2DPoint& p2, Math::C2DPoint& result, int& curve_num)
{
	for (size_t i = 0; i < GetCurveCount(); i++)
	{
		C2DCurve* p_curve = GetCurve(i);
		int ret = p_curve->GetIntersection(p1, p2, result);

		if (ret == -1)
			return -1;

		if (ret != 2)
		{
			curve_num = i;
			return ret;
		}
	}

	return 2;
}

// записывает в result точку и в curve_num номер кривой, ближайшую к указанной и возвращает:
// -1 - ошибка
//  0 - ближайшая точка является начальным узлом кривой
//  1 - ближайшая точка является конечным узлом кривой
//  2 - ближайшая точка лежит где-то на кривой
int C2DOutline::GetClosestPoint(const Math::C2DPoint& p, Math::C2DPoint& result, int& curve_num)
{
	if (!GetCurveCount())
		return -1;

	curve_num = 0;

	//в качестве начального минимума берём расстояние от курсора до (.) на 0й кривой
	C2DCurve* p_curve = GetCurve(0);
	int ret = p_curve->GetClosestPoint(p, result);

	if (ret == -1)
		return ret;

	DBL l = (p - result).Norm2();
	int ret1 = ret;
	Math::C2DPoint tmp_res;

	for (size_t i = 1; i < GetCurveCount(); i++)
	{
		p_curve = GetCurve(i);
		ret = p_curve->GetClosestPoint(p, tmp_res);

		if (ret == -1)
			return ret;

		DBL l1 = (p - tmp_res).Norm2();

		if (l1 < l)
		{
			curve_num = i;
			l = l1;
			ret1 = ret;
			result = tmp_res;
		}
	}

	return ret1;
}


int C2DOutline::GetClosestNode(const Math::C2DPoint& p, DBL& len)
{
	if (!GetNodeCount())
		return -1;

	//в качестве первой "ближайшего" узла берём 0-ю.
	int node_num = 0;
	C2DNode* p_node = GetNode(0);
	//вычисляем расстояние первого узла (НОРМА)
	len = (p - p_node->GetPoint()).Norm();

	Math::C2DPoint t_point;

	for (size_t i = 1; i < GetNodeCount(); i++)
	{
		t_point = GetNode(i)->GetPoint();

		DBL l1 = (p - t_point).Norm();

		if (l1 < len) {
			len = l1;
			node_num = i;
		}
	}

	return node_num;

}

DBL C2DOutline::GetSmallestCurveLenght()
{
	if (GetContourCount() < 1) {
		return 0.0;
	}
	if (GetContour(0)->GetCurveCount() < 1) {
		return 0.0;
	}
	//Длину первой кривой принимаем в качестве начальной минимальной
	DBL dMinLength = GetContour(0)->GetCurve(0)->CalcLength();

	//Проходим по всем контурам и всем кривым
	for (size_t i = 0; i < GetContourCount(); i++) {
		C2DContour *pCont = GetContour(i);
		for (size_t j = 0; j < pCont->GetCurveCount(); j++) {
			DBL dLen = pCont->GetCurve(j)->CalcLength();	//длина кривой
			if (dMinLength > dLen) {
				dMinLength = dLen;
			}
		}
	}
	return dMinLength;
}

int C2DOutline::GetSmallestCurveIndex()
{
	if (GetContourCount() < 1) return -1;

	//Первую кривую принимаем в качестве начальной минимальной
	DBL dMinLength = GetContour(0)->GetCurve(0)->CalcLength();
	int nMinIndex = 0;

	//Проходим по всем контурам и всем кривым
	for (size_t i = 0; i < GetContourCount(); i++) {
		C2DContour *pCont = GetContour(i);
		for (size_t j = 0; j < pCont->GetCurveCount(); j++) {
			if (dMinLength > GetContour(i)->GetCurve(j)->CalcLength()) {
				nMinIndex = GetCurveIndexByCurve(GetContour(i)->GetCurve(j));
			}
		}
	}
	return nMinIndex;
}

Math::C2DRect C2DOutline::GetBoundingBox() {
	Math::C2DRect rect;

	for (size_t i = 0; i < m_nodes().size(); i++) {
		rect.AddRect((dynamic_cast<C2DNode*>(m_nodes()[i]))->GetBoundingBox());
	}

	for (size_t i = 0; i < m_curves().size(); i++) {
		rect.AddRect((dynamic_cast<C2DCurve*>(m_curves()[i]))->GetBoundingBox());
	}

	return rect;
}

bool C2DOutline::GetBoundingBox(CRect2D &rect) {
	if (m_nodes().size() == 0 && m_curves().size() == 0) return false;

	//Начальное состояние прямоугольника - вокруг первого узла
	//Иначе сравниваем с прямоугольником (0,0,0,0).
	rect = (dynamic_cast<C2DNode*>(m_nodes()[0]))->GetBoundingBox();

	//Добавляем оставшееся
	for (size_t i = 0; i < m_nodes().size(); i++) {
		rect.AddRect((dynamic_cast<C2DNode*>(m_nodes()[i]))->GetBoundingBox());
	}

	for (size_t i = 0; i < m_curves().size(); i++) {
		rect.AddRect((dynamic_cast<C2DCurve*>(m_curves()[i]))->GetBoundingBox());
	}

	return true;
}

void C2DOutline::DrawGL(GLParam &parameter)
{

	//glColor3d(0,0,0);
	for (size_t i = 0; i < m_curves().size(); i++) {
		GetCurve(i)->DrawGL(parameter);
	}
		

	// затем узлы
	glPointSize(5);
	for (size_t i = 0; i < m_nodes().size(); i++) {
		GetNode(i)->DrawGL(parameter);
	}
	glPointSize(1);

	// затем контура
	for (size_t i = 0; i < m_contours().size(); i++) {
		GetContour(i)->DrawGL(parameter);
	}
		

}

DBL C2DOutline::GetSumLengthOfAllCurves() {

	DBL dSum = 0.0;
	for (size_t i = 0; i < GetCurveCount(); i++) {
		dSum += GetCurve(i)->CalcLength();
	}
	return dSum;
}

CString C2DOutline::SaveMetaToString(size_t precisizon)
{
	CString sMetaToSave;
	CString sMetaFormatTmp;
	//Проходим по всем узлам
	for (size_t i = 0; i < GetNodeCount(); i++) {
		double dx = GetNode(i)->GetPoint().x;
		double dy = GetNode(i)->GetPoint().y;
		//проблема с унарными и бинарными операторами при загрузке meta файла...
		//setting precision 12
		sMetaFormatTmp.Format(_T("point p%d ( %."+AllToString(precisizon)+"f, %." + AllToString(precisizon) + "f );\n"), i, dx, dy);		
		
		sMetaToSave += sMetaFormatTmp;
	}
	//Enter
	sMetaToSave += _T("\n");
	CString sLineOrArc = _T("");
	CString sRadius = _T("");

	//Проходим по всем кривым (отрезкам и дугам)
	size_t nOutlineCurveCount = GetCurveCount();
	for (size_t i = 0; i < nOutlineCurveCount; i++) {
		int nPointStart = GetCurve(i)->GetStart();
		int nPointEnd = GetCurve(i)->GetEnd();

		//смотрим тип кривой
		if (GetCurve(i)->GetType() == C2DCIRCARC) {
			DBL dRad = ((C2DCircleArc*)(GetCurve(i)))->GetRadius();
			sLineOrArc = _T("arc A");
			//sRadius = _T(", ") + AllToString(dRad + EPS);
			sRadius.Format(_T(", %." + AllToString(precisizon) + "f"), dRad);		//setting precision 12
		}
		else {
			sLineOrArc = _T("line L");
			sRadius.Empty();
		}
		sMetaToSave += sLineOrArc + AllToString(i) + _T(" (p") + AllToString(nPointStart) + _T(", p") + AllToString(nPointEnd) + sRadius + _T(");\n");
	}
	//Enter
	sMetaToSave += _T("\n");
	CString sZapt = _T(", ");
	CString sL_Or_A = _T("L");

	//Проходим по всем контурам (смотрим отрезки и дуги)
	size_t nContourCount = GetContourCount();
	for (size_t i = 0; i < nContourCount; i++) {
		sMetaToSave += _T("contour Loop") + AllToString(i) + _T(" (");	//начало строки контура
		size_t nCurveCount = GetContour(i)->GetCurveCount();
		for (size_t j = 0; j < nCurveCount; j++) {

			int nCurve = GetContour(i)->GetCurveIndexByIndex(j);
			if (j == nCurveCount - 1) {
				sZapt.Empty();
			}
			else {
				sZapt = _T(", ");
			}

			//смотрим тип кривой
			if (GetCurve(j)->GetType() == C2DCIRCARC) {
				sL_Or_A = _T("A");
			}else{
				sL_Or_A = _T("L");
			}
			sMetaToSave += sL_Or_A + AllToString(nCurve) + sZapt;
		}
		sMetaToSave += _T(");\n");	//конец строки контура
	}

	return sMetaToSave;
}

bool C2DOutline::Load(CStorage& file)
{
	if (!IO::CClass::Load(file))
		return false;

	SetParentToAll();

	return true;
}

bool C2DOutline::Copy(IO::CInterface* ptr)
{
	if (!IO::CClass::Copy(ptr))
		return false;

	SetParentToAll();

	return true;
}


void C2DOutline::WriteToLog()
{
	DLOG(CString(_T("C2DOutline<Type:")) + AllToString(GetType()) + CString(_T(">")), log_info);

	for (size_t i = 0; i < GetNodeCount(); i++) {
		DLOG(CString(_T("node<")) + AllToString(GetNode(i)->GetType()) + CString(_T(">[")) + AllToString(i) + CString(_T("]: ")) + AllToString(GetNode(i)->GetPoint().x) + CString(_T(", ")) + AllToString(GetNode(i)->GetPoint().y), log_info);
	}

	for (size_t i = 0; i < GetCurveCount(); i++) {
		DLOG(CString(_T("curve<")) + AllToString(GetCurve(i)->GetType()) + CString(_T(">[")) + AllToString(i) + CString(_T("]: ")) + AllToString(GetCurve(i)->GetStart()) + CString(_T(", ")) + AllToString(GetCurve(i)->GetEnd()), log_info);
	}

	//int c = GetContourCount();
	for (size_t i = 0; i < GetContourCount(); i++)
	{
		CString tmp;
		tmp = CString(_T("contour<")) + AllToString(GetContour(i)->GetType()) + CString(_T(">[")) + AllToString(i) + CString(_T("]: "));

		for (size_t j = 0; j < GetContour(i)->GetCurveCount(); j++)
			tmp = tmp + AllToString(GetContour(i)->At(j)) + CString(_T(", "));

		if (GetContour(i)->IsClosed()) {
			tmp = tmp + CString(_T("closed"));
		}
		else {
			tmp = tmp + CString(_T("opened"));
		}

		tmp = tmp + CString(_T(", total:")) + AllToString(GetContour(i)->GetCurveCount()) + CString(_T("."));

		DLOG(tmp, log_info);
	}
}

void C2DOutline::ClearAll()
{
	m_nodes.DeleteObjects();
	m_nodes().clear();
	m_curves.DeleteObjects();
	m_curves().clear();
	m_contours.DeleteObjects();
	m_contours().clear();
}

bool C2DOutline::GetAllIntersections(size_t c, const Math::C2DLine& line, std::vector<Math::C2DPoint>& points)
{
	if (GetContourCount() < c) {
		return 0;
	}

	C2DContour* pCo = GetContour(c);
	std::vector<Math::C2DPoint> curpoints;

	for (size_t i = 0; i < pCo->GetCurveCount(); i++)
	{
		C2DCurve* pCu = pCo->GetCurve(i);
		curpoints.clear();
		bool ret = pCu->GetAllIntersections(line, curpoints);
		if (ret)
		{
			for (size_t j = 0; j < curpoints.size(); j++)
				points.push_back(curpoints[j]);
		}
	}

	return (points.size() != 0);
}

// проверяет, находится ли точка p внутри контура c (если он замкнут)
// возвращает -2 если ошибка, -1 если точка снаружи, 0 - если на контуре, 1 - если внутри
int C2DOutline::IsInsideContour(const Math::C2DPoint& p, int c)
{
	C2DContour* pcontour = GetContour(c);
	if (!pcontour)
		return -2;

	//size_t nHead, nTail;
	//bool bCheck = pcontour->CheckContour(nHead, nTail);	//for debug only

	if (!pcontour->IsClosed()) // не замкнутый контур => выход
		return -2;

	//перезаполняем кэш контура
	pcontour->FillCache();

	UNINT nTop = 0, nBottom = 0;

	Math::C2DPoint DP1 = pcontour->GetCache()[0];

	bool bLeft = DP1.x < p.x ? true : false;
	size_t r = pcontour->GetCache().size();

	while (r--)
	{
		if (p == DP1) return 0;

		Math::C2DPoint DP2 = DP1;
		DP1 = pcontour->GetCache()[r];

		if ((DP1.x < p.x && !bLeft) || (DP1.x >= p.x && bLeft))
		{
			bLeft = !bLeft;
			double dVal = (DP1.x - DP2.x)*((DP1.y - DP2.y)*(p.x - DP2.x) - (DP1.x - DP2.x)*(p.y - DP2.y));
			if (dVal > 0) {
				nTop++;
			}
			else if (dVal < 0) {
				nBottom++;
			}
			else return 0;
		}
		else if (fabs(DP1.x - p.x) < EPS && fabs(DP2.x - p.x) < EPS && (p.y - DP1.y)*(p.y - DP2.y) <= 0) {
			return 0;
		}
	}

	return 2 * (nTop % 2) - 1;



	//проходим по контуру кроме последней точки
	/*for(int i=0; i<r; i++)
	{
	DP1 = pcontour->GetCache()[i];
	DP2 = pcontour->GetCache()[(i+1 < r ? i+1 : i+1-r)];

	if( (DP1.y >= p.y && DP2.y <= p.y) || (DP1.y <= p.y && DP2.y >= p.y) ){

	double A=0,B=0,C=0;
	A=DP1.y-DP2.y;
	B=DP2.x-DP1.x;
	C=DP1.x*DP2.y-DP2.x*DP1.y;

	if( (-C-B*p.y)/A > p.x ) N++;
	}

	if (DP1.x == p.x && DP2.x == p.x){
	N++;
	if ((p.y - DP1.y) * (p.y - DP2.y) <= 0)
	return 0;
	}
	if (DP1.y == p.y && DP2.y == p.y){
	N++;
	if ((p.x - DP1.x) * (p.x - DP2.x) <= 0)
	return 0;
	}
	}

	if(N%2 == 0) return -1;	//не в контуре
	else return 1;			//в контуре
	//*/

	/*
	while(r--)
	{
	if(p == DP1)
	return 0;

	Math::C2DPoint DP2 = DP1;
	DP1 = pcontour->GetCache()[r];

	if((DP1.m_x < p.x && !bLeft) || (DP1.m_x >= p.x && bLeft)){

	bLeft = !bLeft;
	//double y0 = (DP1.m_y-DP2.m_y)/(DP1.m_x-DP2.m_x)*(x-DP2.m_x) + DP2.m_y;
	//double dVal = (DP1.m_y-DP2.m_y)*(x-DP2.m_x)*(DP1.m_x-DP2.m_x) > (DP1.m_x-DP2.m_x)*(DP1.m_x-DP2.m_x)*(y - DP2.m_y);
	double dVal = (DP1.m_y - DP2.m_y) * (p.x - DP2.m_x) - (DP1.m_x - DP2.m_x) * (p.y - DP2.m_y);

	if(dVal > 0) {
	nTop++;
	}else if(dVal < 0){
	nBottom++;
	}else return 0;
	}
	else if (DP1.m_x == p.x && DP2.m_x == p.x && (p.y - DP1.m_y) * (p.y - DP2.m_y) <= 0)
	return 0;
	}
	return 2 * (nTop%2) - 1;
	//*/

}

/*
int CFrontier::HaveInside(double x, double y, int nLoop)
{// возвращает -1 если точка снаружи, 0 - если на контуре, 1 - если внутри

if(!m_Loops[nLoop].m_bLoop) return -1;
int nTop=0, nBottom=0;
Math::C2DPoint DP1 = GetPoint(nLoop,0), DP(x,y);
bool bLeft = DP1.m_x < x ? true : false;
int r = m_Loops[nLoop].m_nSize;
while(r--)
{
if(DP==DP1)return 0;
Math::C2DPoint DP2 = DP1;
DP1 = GetPoint(nLoop, r);
if((DP1.m_x < x && !bLeft) || (DP1.m_x >= x && bLeft))
{
bLeft = !bLeft;
//double y0 = (DP1.m_y-DP2.m_y)/(DP1.m_x-DP2.m_x)*(x-DP2.m_x) + DP2.m_y;
//double dVal = (DP1.m_y-DP2.m_y)*(x-DP2.m_x)*(DP1.m_x-DP2.m_x) > (DP1.m_x-DP2.m_x)*(DP1.m_x-DP2.m_x)*(y - DP2.m_y);
double dVal = (DP1.m_y-DP2.m_y)*(x-DP2.m_x) - (DP1.m_x-DP2.m_x)*(y - DP2.m_y);
if(dVal > 0) nTop++;
else if(dVal < 0) nBottom++;
else return 0;
} else if(DP1.m_x == x && DP2.m_x == x && (y - DP1.m_y)*(y - DP2.m_y) <=0 )	{
return 0;
}
}
return 2 * (nTop%2) - 1;
}
*/

// проверяет, находится ли точка p внутри хотя бы одного из контуров
int C2DOutline::IsInside(const Math::C2DPoint& p)
{
	size_t nContours = m_contours().size();
	for (size_t i = 0; i < nContours; i++)
	{
		int ret = IsInsideContour(p, i);
		if ((ret != -1) && (ret != -2))
			return ret;
	}

	return -1;
}

C2DOutline::~C2DOutline()
{
	//for (int i = 0; i < m_nodes().GetCount(); i++)
	//	if (m_nodes()[i])
	//		delete m_nodes()[i];
	m_nodes.DeleteObjects();

	//for (int i = 0; i < m_curves().GetCount(); i++)
	//	if (m_curves()[i])
	//		delete m_curves()[i];
	m_curves.DeleteObjects();

	//for (int i = 0; i < m_contours().GetCount(); i++)
	//	if (m_contours()[i])
	//		delete m_contours()[i];
	m_contours.DeleteObjects();

}

// Фаска, кодил Юра
/* возвращает 2 вектора, выходящих из угловой точки */
bool C2DFacet_common::Facet_ugl(C2DPoint *vect1, C2DPoint *vect2, C2DPoint p1b, C2DPoint p1e, C2DPoint p2b, C2DPoint p2e)
{
	if (p1b == p2b && p2b == m_ugl) {
		*vect1 = p1e - m_ugl;
		*vect2 = p2e - m_ugl;
		return true;
	}
	//Если ugl совпадает с началом 1ого отрезка и концом 2ого
	if (p1b == p2e && p2e == m_ugl) {
		*vect1 = p1e - m_ugl;
		*vect2 = p2b - m_ugl;
		return true;
	}

	//Если ugl совпадает с концом 1ого и началом 2ого отрезка
	if (p1e == p2b && p2b == m_ugl) {
		*vect1 = p1b - m_ugl;
		*vect2 = p2e - m_ugl;
		return true;
	}//Если ugl совпадает с концами 1ого и 2ого отрезков
	if (p1e == p2e && p2e == m_ugl) {
		*vect1 = p1b - m_ugl;
		*vect2 = p2b - m_ugl;
		return true;
	}
	return false;
}

/* возвращает повёрнутый на 90 и смасштабированный вектор */
Math::C2DPoint C2DFacet_common::Facet_scale_rotate(C2DPoint pbegin, C2DPoint p1, C2DPoint p2, C2DPoint p3) {

	pbegin.Normalize();
	pbegin *= fabs(m_dRad);

	if (p1.Right(p2, p3)) {
		pbegin = pbegin.Rotate(M_PI / 2, true);	//90 градусов по часовой (y,-x)
	}
	else {
		pbegin = pbegin.Rotate(M_PI / 2, false);	//90 градусов против часовой (-y,x)
	}
	return pbegin;
}

/* cоздаёт дугу и возвращает её номер */
int C2DFacet_common::Facet_create_arc(C2DPoint p1, C2DPoint p2, C2DPoint p3, int n11, int n21, DBL rad) {

	int c;
	if (p1.Right(p2, p3)) {
		c = CreateCurve(C2DCIRCARC, n11, n21);
	}
	else {
		c = CreateCurve(C2DCIRCARC, n21, n11);
	}

	if (c != -1) {
		C2DCircleArc* pArc = dynamic_cast<C2DCircleArc *>(GetCurve(c));
		pArc->SetRadius(fabs(rad));
		pArc->SetCW(false);
	}
	else {
		return c;
	}

	return c;
}

// Выстраивает узлы в нужном порядке для вспомогательного контура
void C2DFacet_common::SetTempFacetPointsOrder(C2DPoint s1, C2DPoint e1, C2DPoint s2, C2DPoint e2) {

	if (s1 == s2 && s2 == m_ugl) {
		m_nNode1 = CreateNode(C2DNODE, e1);
		m_nNode2 = CreateNode(C2DNODE, s1);
		m_nNode3 = CreateNode(C2DNODE, e2);
	}
	if (s1 == e2 && e2 == m_ugl) {
		m_nNode1 = CreateNode(C2DNODE, e1);
		m_nNode2 = CreateNode(C2DNODE, s1);
		m_nNode3 = CreateNode(C2DNODE, s2);
	}
	if (e1 == s2 && s2 == m_ugl) {
		m_nNode1 = CreateNode(C2DNODE, s1);
		m_nNode2 = CreateNode(C2DNODE, e1);
		m_nNode3 = CreateNode(C2DNODE, e2);
	}
	if (e1 == e2 && e2 == m_ugl) {
		m_nNode1 = CreateNode(C2DNODE, s1);
		m_nNode2 = CreateNode(C2DNODE, e1);
		m_nNode3 = CreateNode(C2DNODE, s2);
	}
	//случай когда точек две
	if (s1 == s2 && e1 == e2 && s1 == m_ugl) {
		m_nNode1 = CreateNode(C2DNODE, s1);
		m_nNode2 = CreateNode(C2DNODE, e1);
		m_nNode3 = CreateNode(C2DNODE, s2);
	}
	if (s1 == e2 && e1 == s2 && s1 == m_ugl) {
		m_nNode1 = CreateNode(C2DNODE, s1);
		m_nNode2 = CreateNode(C2DNODE, e1);
		m_nNode3 = CreateNode(C2DNODE, e2);
	}
}

/*void C2DFacet_common::CorrectArcs(short t1, short t2, Math::C2DPoint s1, Math::C2DPoint e1, Math::C2DPoint s2, Math::C2DPoint e2){

if ( t1 == C2DCIRCARC && ((s1==s2) || (s1==e2)) ){
((C2DCircleArc *) GetCurve(0))->SwapNodes();
}
if ( t2 == C2DCIRCARC && ((s1!=s2) && (s2!=e1)) ){
((C2DCircleArc *) GetCurve(1))->SwapNodes();
}

//Если отрезок-дуга замкнуты (полумесяц)
short type;
if ((t1 == C2DCIRCARC && t2 == C2DCURVE) || (t1 == C2DCURVE && t2 == C2DCIRCARC)){
if (t1 == C2DCIRCARC){
type = 0;
}else{
type = 1;
}

if ( ((s2==s1) && (e1==e2)) || ((e2==s1) && (e1==s2)) ){
((C2DCircleArc *) GetCurve(type))->SwapNodes();
}
}

//Если дуга-дуга замкнуты (месяц или овал)
if (t1==C2DCIRCARC && t2==C2DCIRCARC){
if ( (s2==s1) && (e1==e2) ){
((C2DCircleArc *) GetCurve(0))->SwapNodes();
((C2DCircleArc *) GetCurve(1))->SwapNodes();
}
// Второй случай
if ( (e2==s1) && (e1==s2) ){
((C2DCircleArc *) GetCurve(0))->SwapNodes();
}
}
}//*/

// ********2_LINES***************
bool C2DFacet_2lines::Create() {

	//Создание векторов перемещения
	Math::C2DPoint a, b; //вектора для параллельного перемещения
	a = Facet_scale_rotate(m_pVect2, m_ugl + m_pVect1, m_ugl, m_ugl + m_pVect2);
	b = Facet_scale_rotate(m_pVect1, m_ugl + m_pVect2, m_ugl, m_ugl + m_pVect1);

	//C2DOutline Outline;
	//Создание вспомогательного отрезка, параллельного исходному
	int n1 = CreateNode(C2DNODE, m_ugl + b),		//номера вспомогательных узлов
		n2 = CreateNode(C2DNODE, m_ugl + b + m_pVect1),	//номер вспомогательных кривых
		c = CreateCurve(C2DCURVE, n1, n2);

	Math::C2DPoint tper, per1, per2; //точки пересечения

	if (GetCurve(c)->GetIntersection(m_ugl + a, m_ugl + a + m_pVect2, tper) != 1) { return false; };

	//Находим точки на отрезках, в которых будет начало дуги
	per1 = tper - a;
	per2 = tper - b;

	//Удаляем исходные кривые
	DelCurve(c);
	DelCurve(1);
	DelCurve(0);
	//Создаём фаску
	if (!CreateFacet(per1, per2)) { return false; };
	return true;
};

/* AddFacet_2lines создаёт фаску для отрезка-отрезка */
bool C2DFacet_2lines::CreateFacet(C2DPoint per1, C2DPoint per2) {

	int n11 = CreateNode(C2DNODE, per1),
		n12 = CreateNode(C2DNODE, m_ugl + m_pVect2),
		n21 = CreateNode(C2DNODE, per2),
		n22 = CreateNode(C2DNODE, m_ugl + m_pVect1);
	if ((n11 == -1) || (n12 == -1) || (n21 == -1) || (n22 == -1)) return false;

	//Если проекция опустилась на концы векторов, то строим ТОЛЬКО дугу
	if (per1 == (m_ugl + m_pVect2) && per2 == (m_ugl + m_pVect1)) {
		int c3 = Facet_create_arc(m_ugl + m_pVect1, m_ugl, m_ugl + m_pVect2, n11, n21, m_dRad);
		if (c3 == -1) { return false; };
		return true;
	}

	//Создание кривых
	int c1 = CreateCurve(C2DCURVE, n11, n12),
		c2 = CreateCurve(C2DCURVE, n21, n22),
		c3 = Facet_create_arc(m_ugl + m_pVect1, m_ugl, m_ugl + m_pVect2, n11, n21, m_dRad);

	if ((c1 == -1) || (c2 == -1) || (c3 == -1)) { return false; };
	return true;
}
void C2DFacet_2lines::CreateData() {

	CreateCurve(C2DCURVE, m_nNode1, m_nNode2);
	CreateCurve(C2DCURVE, m_nNode2, m_nNode3);
}

// ********LINE_ARC***************
// Определяет точки для случая отрезок-дуга.
bool C2DFacet_line_arc::Create() {

	C2DCurve* v1 = GetCurve(m_Cur);
	C2DCircleArc* v2 = dynamic_cast<C2DCircleArc*>(GetCurve(m_Arc));

	Math::C2DPoint o1;			//центр исходной дуги
	DBL r = v2->GetRadius();	//радиус исходной дуги
	Math::C2DPoint p3, tau;		//TODO: точка на середине дуги (центр дуги по 3м точкам)
								//bool bCenter = v2->GetPoint(v2->CalcLength()/2,p3,tau);

	if (!CircleCenter(v2->GetStartNode()->GetPoint(), v2->GetEndNode()->GetPoint(), r, o1)) return false; //Получаем центр дуги
	Math::C2DPoint l2, o2; //середина отрезка, дуги и касательная
	if (!v2->GetPoint(v2->CalcLength() / 2, o2, tau)) return false;

	if (m_Arc == 0) {
		tau = m_pVect2;
		m_pVect2 = m_pVect1;
		m_pVect1 = tau;
	}

	if (!v1->GetPoint(v1->CalcLength() / 2, l2, tau)) { return false; }
	if (((m_ugl + m_pVect1) == (m_ugl + m_pVect2)) && (o2.Len(l2) < m_dRad*2.0)) {
		return false; //если расстояние в замкнутой области между отрезком и дугой мало
	}

	//Создание векторов копирования отрезка
	Math::C2DPoint a; //вектора для параллельного перемещения отрезка
	a = Facet_scale_rotate(m_pVect1, o2, m_ugl, m_ugl + m_pVect1);

	//Точки 2х вспомогательных отрезков
	m_P1up = m_ugl + a;
	m_P2up = m_ugl + a + m_pVect1,
		m_P1down = m_ugl - a,
		m_P2down = m_ugl - a + m_pVect1;

	//Вспомогательные вектора для дуги
	Math::C2DPoint	tmpvect3 = m_ugl - o1,
		tmpvect4 = (m_ugl + m_pVect2) - o1;
	tmpvect3.Normalize();
	tmpvect4.Normalize();

	//Узлы новых дуг (расширение\сужение)
	int narc1up = CreateNode(C2DNODE, o1 + tmpvect3*fabs(m_dRad + r)),
		narc2up = CreateNode(C2DNODE, o1 + tmpvect4*fabs(m_dRad + r)),
		narc1down = CreateNode(C2DNODE, o1 + tmpvect3*fabs(m_dRad - r)),
		narc2down = CreateNode(C2DNODE, o1 + tmpvect4*fabs(m_dRad - r));

	//Создание вспомогательных дуг
	int carc_up = Facet_create_arc(o2, m_ugl, m_ugl + m_pVect2, narc1up, narc2up, fabs(m_dRad + r));
	int carc_down = Facet_create_arc(o2, m_ugl, m_ugl + m_pVect2, narc1down, narc2down, fabs(m_dRad - r));

	//Нахождение центра фаски
	Math::C2DPoint result;
	if ((m_ugl + m_pVect1) != (m_ugl + m_pVect2)) {
		if (!Facet_line_arc_point(carc_up, carc_down, &result)) {
			/*DelCurve(carc_down); DelCurve(carc_up);
			DelNode(narc2down); DelNode(narc1down); DelNode(narc2up); DelNode(narc1up); */
			return false;
		}
	}
	else { // возможна замкнутая область
		if (!Facet_line_arc_point(carc_up, carc_down, &result)) {
			//если не одна точка пересечения
			if (!Facet_line_arc_2points(carc_up, carc_down, &result)) {
				/*DelCurve(carc_down); DelCurve(carc_up);
				DelNode(narc2down); DelNode(narc1down); DelNode(narc2up); DelNode(narc1up); */
				return false;
			}
		}
	}
	//Точки пересечения с исходными кривыми
	Math::C2DPoint per1, per2;
	v1->GetClosestPoint(result, per1);
	v2->GetClosestPoint(result, per2);

	if (per2 == m_ugl || per1 == m_ugl) {
		DelCurve(carc_down); DelCurve(carc_up);
		DelNode(narc2down); DelNode(narc1down); DelNode(narc2up); DelNode(narc1up);
		return false;
	}
	//Удаляем вспомогательные отрезки, создаём 2 новых отрезка и дугу
	DelCurve(carc_down); DelCurve(carc_up);
	DelNode(narc2down); DelNode(narc1down); DelNode(narc2up); DelNode(narc1up);

	//Создаём фаску
	if (!CreateFacet(per1, per2, o2)) { return false; };
	//Удаляем исходные кривые
	DelCurve(1); DelCurve(0);
	return true;
}
/* Facet_line_arc_point возвращает центр фаски для отрезка-дуги */
bool C2DFacet_line_arc::Facet_line_arc_point(int up, int down, C2DPoint* result) {

	Math::C2DPoint mas[4], res;
	short count = 0, i;
	C2DCircleArc *pCurve_up = dynamic_cast<C2DCircleArc *>(GetCurve(up)),
		*pCurve_down = dynamic_cast<C2DCircleArc *>(GetCurve(down));

	if (pCurve_up->GetIntersection(m_P1up, m_P2up, res) == 1) {
		mas[count] = res;
		count++;
	};
	if (pCurve_up->GetIntersection(m_P1down, m_P2down, res) == 1) {
		mas[count] = res;
		count++;
	};
	if (pCurve_down->GetIntersection(m_P1up, m_P2up, res) == 1) {
		mas[count] = res;
		count++;
	};
	if (pCurve_down->GetIntersection(m_P1down, m_P2down, res) == 1) {
		mas[count] = res;
		count++;
	};

	if (count == 0) { return false; };

	DBL len = mas[0].Len(m_ugl + m_pVect1);
	res = mas[0];
	for (i = 1; i<count; i++) {
		if (len > mas[i].Len(m_ugl + m_pVect1)) {
			len = mas[i].Len(m_ugl + m_pVect1);
			res = mas[i];
		}
	}
	*result = res;
	return true;
}
/* Facet_line_arc_2points возвращает центр фаски для отрезка-дуги при замыкании */
bool C2DFacet_line_arc::Facet_line_arc_2points(int up, int down, C2DPoint* result) {

	Math::C2DPoint mas[4][2], res1, res2;
	short count1 = 0, count2 = 0,
		i, j;
	C2DCircleArc *pCurve_up = dynamic_cast<C2DCircleArc *>(GetCurve(up)),
		*pCurve_down = dynamic_cast<C2DCircleArc *>(GetCurve(down));

	if (pCurve_up->GetIntersection2Points(m_P1up, m_P2up, res1, res2) == 1) {
		mas[count1][count2] = res1;
		count2++;
		mas[count1][count2] = res2;
		count1++;
	};
	count2 = 0;
	if (pCurve_up->GetIntersection2Points(m_P1down, m_P2down, res1, res2) == 1) {
		mas[count1][count2] = res1;
		count2++;
		mas[count1][count2] = res2;
		count1++;
	};
	count2 = 0;
	if (pCurve_down->GetIntersection2Points(m_P1up, m_P2up, res1, res2) == 1) {
		mas[count1][count2] = res1;
		count2++;
		mas[count1][count2] = res2;
		count1++;
	};
	count2 = 0;
	if (pCurve_down->GetIntersection2Points(m_P1down, m_P2down, res1, res2) == 1) {
		mas[count1][count2] = res1;
		count2++;
		mas[count1][count2] = res2;
		count1++;
	};

	if (count1 == 0) { return false; };

	DBL len = mas[0][0].Len(m_ugl);
	res1 = mas[0][0];
	for (i = 0; i<count1; i++) {
		for (j = 0; j<2; j++) {
			if (len > mas[i][j].Len(m_ugl)) {
				len = mas[i][j].Len(m_ugl);
				res1 = mas[i][j];
			}
		}
	}
	*result = res1;
	return true;
}
/* CreateFacet создаёт фаску для отрезка-дуги */
bool C2DFacet_line_arc::CreateFacet(C2DPoint per1, C2DPoint per2, C2DPoint o2) {

	C2DCircleArc* v2 = dynamic_cast<C2DCircleArc*>(GetCurve(m_Arc));
	int n11 = CreateNode(C2DNODE, per1),
		n12 = CreateNode(C2DNODE, m_ugl + m_pVect1),
		n21 = CreateNode(C2DNODE, per2),
		n22 = CreateNode(C2DNODE, m_ugl + m_pVect2);
	if ((n11 == -1) || (n12 == -1) || (n21 == -1) || (n22 == -1)) { return false; }

	//Создание кривых
	int c1, c2, c3;
	c1 = CreateCurve(C2DCURVE, n11, n12);
	c2 = Facet_create_arc(o2, m_ugl, m_ugl + m_pVect1, n21, n22, v2->GetRadius());

	//Проверка на правильную выгнутость исходной дуги
	C2DCircleArc* vc3 = dynamic_cast<C2DCircleArc*>(GetCurve(c2));
	Math::C2DPoint o3, o4, tmp;
	if (!vc3->GetPoint((vc3->CalcLength()) / 2, o3, tmp)) {
		return false;
	}
	if (v2->GetClosestPoint(o3, o4) == -1) { return false; }
	if (o4.Len(o3) >= EPS) {
		vc3->SwapNodes();
		vc3->SetCW(!vc3->GetCW());
	}

	c3 = Facet_create_arc(o2, m_ugl, m_ugl + m_pVect1, n11, n21, m_dRad);

	if ((c1 == -1) || (c2 == -1) || (c3 == -1)) { return false; };
	return true;
}
/* CreateData создаёт кривые */
void C2DFacet_line_arc::CreateData(short t1, short t2, bool cw, DBL r) {

	int n1, n2;

	n1 = CreateCurve(t1, m_nNode1, m_nNode2);
	n2 = CreateCurve(t2, m_nNode2, m_nNode3);

	if (t1 == C2DCURVE) {
		m_Cur = 0;
		m_Arc = 1;
	}
	else {
		m_Cur = 1;
		m_Arc = 0;
	};

	C2DCircleArc *pArc = dynamic_cast<C2DCircleArc *>(GetCurve(m_Arc));
	pArc->SetCW(cw);
	pArc->SetRadius(r);

}
/* CorrectArcs корректирует дугу */
void C2DFacet_line_arc::CorrectArcs(short type1, short type2, int nCurve1, int nCurve2, C2DOutline* out) {
	Math::C2DPoint tmpF, tmpO, tau;
	C2DCircleArc *pOut;

	if (type1 == C2DCURVE) {
		//тип второй кривой - дуга
		pOut = dynamic_cast<C2DCircleArc *>(out->GetCurve(nCurve2));
		pOut->GetPoint(pOut->CalcLength() / 2, tmpO, tau);
	}
	else {
		//тип первой кривой - дуга
		pOut = dynamic_cast<C2DCircleArc *>(out->GetCurve(nCurve1));
		pOut->GetPoint(pOut->CalcLength() / 2, tmpO, tau);
	}

	pOut = dynamic_cast<C2DCircleArc *>(GetCurve(m_Arc));
	pOut->GetPoint(pOut->CalcLength() / 2, tmpF, tau);

	//Сверяем построенную дугу с исходной
	if (tmpF != tmpO) {
		pOut->SetCW(!pOut->GetCW());
	}
}

// ********2_ARCS***************
// Определяет точки для случая дуга-дуга.
bool C2DFacet_2arcs::Create() {

	C2DCircleArc* v1 = dynamic_cast<C2DCircleArc*>(GetCurve(0));
	C2DCircleArc* v2 = dynamic_cast<C2DCircleArc*>(GetCurve(1));
	Math::C2DPoint l1, l2, tau; //середины дуг и касательная

	if (!v1->GetPoint((v1->CalcLength()) / 2, l1, tau)) { return false; };
	if (!v2->GetPoint((v2->CalcLength()) / 2, l2, tau)) { return false; };

	if (((m_ugl + m_pVect1) == (m_ugl + m_pVect2)) && (l1.Len(l2) < m_dRad*2.0)) {
		return false; //если расстояние в замкнутой области между дугами мало
	}

	Math::C2DPoint o1, o2;			//центры исходных дуг
	DBL r1 = v1->GetRadius(),
		r2 = v2->GetRadius();	//радиусы исходных дуг

	Math::C2DPoint p03, p13;		//TODO: точка на середине дуги (центр дуги по 3м точкам)
									//bool bCenter1 = v1->GetPoint(v1->CalcLength()/2,p03,tau);
									//bool bCenter2 = v2->GetPoint(v2->CalcLength()/2,p13,tau);

									//Получаем центры дуг
	if (!CircleCenter(v1->GetStartNode()->GetPoint(), v1->GetEndNode()->GetPoint(), r1, o1)) { return false; };
	if (!CircleCenter(v2->GetStartNode()->GetPoint(), v2->GetEndNode()->GetPoint(), r2, o2)) { return false; };

	//Вспомогательные вектора для дуг
	Math::C2DPoint	tmpvect11 = m_ugl - o1,
		tmpvect12 = (m_ugl + m_pVect1) - o1,
		tmpvect21 = m_ugl - o2,
		tmpvect22 = (m_ugl + m_pVect2) - o2;
	tmpvect11.Normalize();
	tmpvect12.Normalize();
	tmpvect21.Normalize();
	tmpvect22.Normalize();

	//Узлы новых дуг (расширение\сужение)
	int narc11up = CreateNode(C2DNODE, o1 + tmpvect11*fabs(m_dRad + r1)),
		narc12up = CreateNode(C2DNODE, o1 + tmpvect12*fabs(m_dRad + r1)),
		narc11down = CreateNode(C2DNODE, o1 + tmpvect11*fabs(m_dRad - r1)),
		narc12down = CreateNode(C2DNODE, o1 + tmpvect12*fabs(m_dRad - r1)),

		narc21up = CreateNode(C2DNODE, o2 + tmpvect21*fabs(m_dRad + r2)),
		narc22up = CreateNode(C2DNODE, o2 + tmpvect22*fabs(m_dRad + r2)),
		narc21down = CreateNode(C2DNODE, o2 + tmpvect21*fabs(m_dRad - r2)),
		narc22down = CreateNode(C2DNODE, o2 + tmpvect22*fabs(m_dRad - r2));

	//Создание вспомогательных дуг
	m_nArc1_up = Facet_create_arc(l1, m_ugl, m_ugl + m_pVect1, narc11up, narc12up, fabs(m_dRad + r1));
	m_nArc1_down = Facet_create_arc(l1, m_ugl, m_ugl + m_pVect1, narc11down, narc12down, fabs(m_dRad - r1));
	m_nArc2_up = Facet_create_arc(l2, m_ugl, m_ugl + m_pVect2, narc21up, narc22up, fabs(m_dRad + r2));
	m_nArc2_down = Facet_create_arc(l2, m_ugl, m_ugl + m_pVect2, narc21down, narc22down, fabs(m_dRad - r2));

	//Нахождение центра фаски
	Math::C2DPoint result;
	//Если дуги не образуют кольцо, то пересечений 1 шт.
	if ((m_ugl + m_pVect1) != (m_ugl + m_pVect2)) {
		if (!Facet_2arcs_point(l1, l2, &result)) {
			/*DelCurve(m_nArc2_down); DelCurve(m_nArc2_up); DelCurve(m_nArc1_down); DelCurve(m_nArc1_up);
			DelNode(narc22down); DelNode(narc21down); DelNode(narc22up); DelNode(narc21up);
			DelNode(narc12down); DelNode(narc11down); DelNode(narc12up); DelNode(narc11up); */
			return false;
		};
	}
	else {
		//Дуги замкнулись, если 1 пересечение
		if (!Facet_2arcs_point(l1, l2, &result)) {
			// если 2 пересечения
			if (!Facet_2arcs_2points(l1, l2, &result)) {
				/*DelCurve(m_nArc2_down); DelCurve(m_nArc2_up); DelCurve(m_nArc1_down); DelCurve(m_nArc1_up);
				DelNode(narc22down); DelNode(narc21down); DelNode(narc22up); DelNode(narc21up);
				DelNode(narc12down); DelNode(narc11down); DelNode(narc12up); DelNode(narc11up); */
				return false;
			};
		};
	};
	//Точки пересечения с исходными кривыми
	Math::C2DPoint per1, per2;
	v1->GetClosestPoint(result, per1);
	v2->GetClosestPoint(result, per2);

	//Удаляем вспомогательные отрезки и узлы
	DelCurve(m_nArc2_down); DelCurve(m_nArc2_up); DelCurve(m_nArc1_down); DelCurve(m_nArc1_up);
	DelNode(narc22down); DelNode(narc21down); DelNode(narc22up); DelNode(narc21up);
	DelNode(narc12down); DelNode(narc11down); DelNode(narc12up); DelNode(narc11up);

	//Удаляем исходные кривые
	DelCurve(1); DelCurve(0);
	if (!CreateFacet(per1, per2, l1, l2, r1, r2)) { return false; };
	return true;
}
/* CreateFacet cоздаёт фаску для дуги-дуги */
bool C2DFacet_2arcs::CreateFacet(C2DPoint per1, C2DPoint per2, C2DPoint l1, C2DPoint l2, DBL r1, DBL r2) {

	int n11 = CreateNode(C2DNODE, per1),
		n12 = CreateNode(C2DNODE, m_ugl + m_pVect1),
		n21 = CreateNode(C2DNODE, per2),
		n22 = CreateNode(C2DNODE, m_ugl + m_pVect2);
	if ((n11 == -1) || (n12 == -1) || (n21 == -1) || (n22 == -1)) return false;

	//Создание кривых
	int c1, c2, c3, c4, c5;

	if (l1.Right(m_ugl, (m_ugl + m_pVect1))) {
		c1 = CreateCurve(C2DCIRCARC, n11, n12);
	}
	else {
		c1 = CreateCurve(C2DCIRCARC, n12, n11);
	}

	if (l2.Right(m_ugl, (m_ugl + m_pVect2))) {
		c2 = CreateCurve(C2DCIRCARC, n21, n22);
	}
	else {
		c2 = CreateCurve(C2DCIRCARC, n22, n21);
	}

	//Фаска
	if (l1.Right(m_ugl, (m_ugl + m_pVect1)) && (!l2.Right(m_ugl, (m_ugl + m_pVect2)))) {
		c3 = CreateCurve(C2DCIRCARC, n11, n21);
	}
	else {
		c3 = CreateCurve(C2DCIRCARC, n21, n11);
	}

	//Если дуги замкнуты
	if (m_ugl + m_pVect1 == m_ugl + m_pVect2) {

		DelCurve(c3);

		c4 = CreateCurve(C2DCIRCARC, n11, n21);
		C2DCircleArc *pCirc4 = dynamic_cast<C2DCircleArc *>(GetCurve(c4));
		pCirc4->SetRadius(fabs(m_dRad));
		pCirc4->SetCW(false);

		c5 = CreateCurve(C2DCIRCARC, n21, n11);
		C2DCircleArc *pCirc5 = dynamic_cast<C2DCircleArc *>(GetCurve(c5));
		pCirc5->SetRadius(fabs(m_dRad));
		pCirc5->SetCW(false);

		Math::C2DPoint pl1, pl2, tau; //середины дуг и касательная
		if (!pCirc4->GetPoint(pCirc4->CalcLength() / 2, pl1, tau)) { return false; };
		if (!pCirc5->GetPoint(pCirc5->CalcLength() / 2, pl2, tau)) { return false; };

		DelCurve(c5);
		DelCurve(c4);

		if (pl1.Len(m_ugl) < pl2.Len(m_ugl)) {
			c3 = CreateCurve(C2DCIRCARC, n11, n21);
		}
		else {
			c3 = CreateCurve(C2DCIRCARC, n21, n11);
		}
	}

	if ((c1 == -1) || (c2 == -1) || (c3 == -1)) return false;

	//Параметры для дуг
	dynamic_cast<C2DCircleArc *>(GetCurve(c1))->SetRadius(fabs(r1));
	dynamic_cast<C2DCircleArc *>(GetCurve(c1))->SetCW(false);
	dynamic_cast<C2DCircleArc *>(GetCurve(c2))->SetRadius(fabs(r2));
	dynamic_cast<C2DCircleArc *>(GetCurve(c2))->SetCW(false);
	dynamic_cast<C2DCircleArc *>(GetCurve(c3))->SetRadius(fabs(m_dRad));
	dynamic_cast<C2DCircleArc *>(GetCurve(c3))->SetCW(false);
	return true;
}
/* Facet_2arcs_point возвращает центр фаски для дуги-дуги */
bool C2DFacet_2arcs::Facet_2arcs_point(C2DPoint l1, C2DPoint l2, C2DPoint* result) {

	Math::C2DPoint mas[4], res;
	short count = 0, i;
	C2DCircleArc *c1_up = dynamic_cast<C2DCircleArc *>(GetCurve(m_nArc1_up)),
		*c1_down = dynamic_cast<C2DCircleArc *>(GetCurve(m_nArc1_down)),
		*c2_up = dynamic_cast<C2DCircleArc *>(GetCurve(m_nArc2_up)),
		*c2_down = dynamic_cast<C2DCircleArc *>(GetCurve(m_nArc2_down));

	if (c1_up->GetIntersectionArcs(*c2_up, res) == 1) {
		mas[count] = res;
		count++;
	};
	if (c1_up->GetIntersectionArcs(*c2_down, res) == 1) {
		mas[count] = res;
		count++;
	};
	if (c1_down->GetIntersectionArcs(*c2_up, res) == 1) {
		mas[count] = res;
		count++;
	};
	if (c1_down->GetIntersectionArcs(*c2_down, res) == 1) {
		mas[count] = res;
		count++;
	};

	if (count == 0) { return false; };

	Math::C2DPoint tmp = (l1 + l2) / 2.0;
	DBL len = mas[0].Len(tmp);
	res = mas[0];
	for (i = 1; i<count; i++) {
		if (len > mas[i].Len(tmp)) {
			res = mas[i];
		}
	}
	*result = res;
	return true;
};
/* Facet_2arcs_2points возвращает центр фаски для дуги-дуги при замыкании */
bool C2DFacet_2arcs::Facet_2arcs_2points(Math::C2DPoint l1, Math::C2DPoint l2, C2DPoint* result) {

	Math::C2DPoint mas[4][2], res1, res2;
	short count1 = 0, count2 = 0,
		i, j;
	C2DCircleArc *c1_up = dynamic_cast<C2DCircleArc *>(GetCurve(m_nArc1_up)),
		*c1_down = dynamic_cast<C2DCircleArc *>(GetCurve(m_nArc1_down)),
		*c2_up = dynamic_cast<C2DCircleArc *>(GetCurve(m_nArc2_up)),
		*c2_down = dynamic_cast<C2DCircleArc *>(GetCurve(m_nArc2_down));

	if (c1_up->GetIntersectionArcs2Points(*c2_up, res1, res2) == 1) {
		mas[count1][count2] = res1;
		count2++;
		mas[count1][count2] = res2;
		count1++;
	};
	count2 = 0;
	if (c1_up->GetIntersectionArcs2Points(*c2_down, res1, res2) == 1) {
		mas[count1][count2] = res1;
		count2++;
		mas[count1][count2] = res2;
		count1++;
	};
	count2 = 0;
	if (c1_down->GetIntersectionArcs2Points(*c2_up, res1, res2) == 1) {
		mas[count1][count2] = res1;
		count2++;
		mas[count1][count2] = res2;
		count1++;
	};
	count2 = 0;
	if (c1_down->GetIntersectionArcs2Points(*c2_down, res1, res2) == 1) {
		mas[count1][count2] = res1;
		count2++;
		mas[count1][count2] = res2;
		count1++;
	};

	if (count1 == 0) { return false; };

	DBL len = mas[0][0].Len(m_ugl);
	res1 = mas[0][0];
	for (i = 0; i<count1; i++) {
		for (j = 0; j<2; j++) {
			if (len > mas[i][j].Len(m_ugl)) {
				len = mas[i][j].Len(m_ugl);
				res1 = mas[i][j];
			}
		}
	}
	*result = res1;
	return true;
}
/* CreateData заполняет C2DFacet_2arcs */
void C2DFacet_2arcs::CreateData(bool cw1, DBL r1, bool cw2, DBL r2) {

	int n1, n2;
	n1 = CreateCurve(C2DCIRCARC, m_nNode1, m_nNode2);
	n2 = CreateCurve(C2DCIRCARC, m_nNode2, m_nNode3);

	//n1 & n2 == (1 OR 0)
	C2DCircleArc *pCur0 = dynamic_cast<C2DCircleArc *>(GetCurve(n1));
	pCur0->SetCW(cw1);
	pCur0->SetRadius(r1);

	C2DCircleArc *pCur1 = dynamic_cast<C2DCircleArc *>(GetCurve(n2));
	pCur1->SetCW(cw2);
	pCur1->SetRadius(r2);
}
/* CorrectArcs корректирует дуги */
void C2DFacet_2arcs::CorrectArcs(int nCurve1, int nCurve2, C2DOutline* out) {

	Math::C2DPoint tmpF1, tmpF2, tmpO1, tmpO2, tau;
	C2DCircleArc *pOut1 = dynamic_cast<C2DCircleArc *>(out->GetCurve(nCurve1));
	C2DCircleArc *pOut2 = dynamic_cast<C2DCircleArc *>(out->GetCurve(nCurve2));

	//Середины созданных дуг
	pOut1->GetPoint(pOut1->CalcLength() / 2, tmpO1, tau);
	pOut2->GetPoint(pOut2->CalcLength() / 2, tmpO2, tau);

	//Середины исходных дуг
	C2DCircleArc *pCur1 = dynamic_cast<C2DCircleArc *>(GetCurve(0));
	C2DCircleArc *pCur2 = dynamic_cast<C2DCircleArc *>(GetCurve(1));

	pCur1->GetPoint(pCur1->CalcLength() / 2, tmpF1, tau);
	pCur2->GetPoint(pCur2->CalcLength() / 2, tmpF2, tau);

	if (tmpF1 != tmpO1) {
		pCur1->SetCW(!pCur1->GetCW());
	}
	if (tmpF2 != tmpO2) {
		pCur2->SetCW(!pCur2->GetCW());
	}
}
//Конец кодоводства Юры