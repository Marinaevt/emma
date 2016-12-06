#include "StdAfx.h"
#include "3DPoint.h"


namespace Math {

	/**/
	// Class C3DPoint //
	/**/

	const C3DPoint C3DPoint::Ortx = C3DPoint(1.0, 0.0, 0.0);
	const C3DPoint C3DPoint::Orty = C3DPoint(0.0, 1.0, 0.0);
	const C3DPoint C3DPoint::Ortz = C3DPoint(0.0, 0.0, 1.0);
	const C3DPoint C3DPoint::Zero = C3DPoint(0.0, 0.0, 0.0);
	//---------------------------------------------------------------------------

	C3DPoint& C3DPoint::operator = (const C3DPoint& p)
	{
		m_x = p.m_x;
		m_y = p.m_y;
		m_z = p.m_z;
		return *this;
	}

	//---------------------------------------------------------------------------
	C3DPoint C3DPoint::operator + (const C3DPoint& p) const
	{
		return C3DPoint(x + p.x, y + p.y, z + p.z);
	}
	//---------------------------------------------------------------------------

	C3DPoint C3DPoint::operator - (const C3DPoint& p) const
	{
		return C3DPoint(x - p.x, y - p.y, z - p.z);
	}
	//---------------------------------------------------------------------------

	C3DPoint C3DPoint::operator *  (DBL alpha) const
	{
		return C3DPoint(alpha * x, alpha * y, alpha * z);
	}
	//---------------------------------------------------------------------------

	C3DPoint C3DPoint::operator /  (DBL alpha) const
	{
		DBL ialpha = 1 / alpha;

		return C3DPoint(ialpha * x, ialpha * y, ialpha * z);
	}
	//---------------------------------------------------------------------------

	C3DPoint operator * (DBL alpha, const C3DPoint& p)
	{
		return C3DPoint(alpha * p.x, alpha * p.y, alpha * p.z);
	}
	//---------------------------------------------------------------------------

	C3DPoint& C3DPoint::operator += (const C3DPoint& p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}
	//---------------------------------------------------------------------------

	C3DPoint& C3DPoint::operator -= (const C3DPoint& p)
	{
		x -= p.x;
		y -= p.y;
		z -= p.z;
		return *this;
	}
	//---------------------------------------------------------------------------

	C3DPoint& C3DPoint::operator *= (DBL alpha)
	{
		x *= alpha;
		y *= alpha;
		z *= alpha;
		return *this;
	}
	//---------------------------------------------------------------------------

	C3DPoint& C3DPoint::operator /= (DBL alpha)
	{
		DBL ialpha = 1 / alpha;
		x *= ialpha;
		y *= ialpha;
		z *= ialpha;
		return *this;
	}
	//---------------------------------------------------------------------------

	DBL C3DPoint::operator * (const C3DPoint& p) const // скалярное произведение
	{
		return x * p.x + y * p.y + z * p.z;
	}
	//---------------------------------------------------------------------------

	bool C3DPoint::operator == (const C3DPoint& p) const
	{
		return (fabs(x - p.x) < EPS) && (fabs(y - p.y) < EPS) && (fabs(z - p.z) < EPS);
	}
	//---------------------------------------------------------------------------

	bool C3DPoint::operator != (const C3DPoint& p) const
	{
		return (fabs(x - p.x) > EPS) || (fabs(y - p.y) > EPS) || (fabs(z - p.z) > EPS);
	}
	//---------------------------------------------------------------------------

	C3DPoint::operator const DBL*() const
	{
		return c;
	}
	//---------------------------------------------------------------------------

	C3DPoint::operator bool() const
	{
		return this->operator != (Zero);
	}
	//---------------------------------------------------------------------------

	bool C3DPoint::operator ! () const
	{
		return this->operator == (Zero);
	}
	//---------------------------------------------------------------------------

	DBL C3DPoint::Norm() const  // норма
	{
		return sqrt(x * x + y * y + z * z);
	}
	//---------------------------------------------------------------------------

	DBL C3DPoint::Norm2() const // норма ^ 2
	{
		return x * x + y * y + z * z;
	}
	//---------------------------------------------------------------------------

	void C3DPoint::Normalize()
	{
		DBL norm_ = Norm();	// всегда >= 0

		if (norm_ > EPS)
		{
			DBL inorm = 1 / norm_;

			x *= inorm;
			y *= inorm;
			z *= inorm;
		}
	}
	//---------------------------------------------------------------------------

	C3DPoint C3DPoint::Rotate(DBL phi, bool ccw) const
	{
		/*
		DBL s, c;
		if (!ccw) {
			phi = - phi;
		}
		//_fm_sincos(phi, &s, &c);
		s = sin(phi);
		c = cos(phi);

		return C3DPoint(x * c + y * s, -x * s + y * c);
		*/
		return C3DPoint(0, 0, 0);
	}
	//---------------------------------------------------------------------------

	C3DPoint C3DPoint::Min(const C3DPoint& p) const // возвращает вектор с минимальными  координатами
	{
		return C3DPoint(min(x, p.x), min(y, p.y), min(z, p.z));
	}
	//---------------------------------------------------------------------------

	C3DPoint C3DPoint::Max(const C3DPoint& p) const // возвращает вектор с максимальными координатами
	{
		return C3DPoint(max(x, p.x), max(y, p.y), max(z, p.z));
	}
	//---------------------------------------------------------------------------

	C3DPoint C3DPoint::Lerp(const C3DPoint& p, DBL t) const
	{
		return C3DPoint((p.x - x) * t + x,
			(p.y - y) * t + y,
			(p.z - z) * t + z);
	}
	//---------------------------------------------------------------------------
	C3DPoint& C3DPoint::FromAngle(DBL angle)
	{
		x = cos(angle);
		y = sin(angle);
		return *this;
	}
	//---------------------------------------------------------------------------


	C3DPoint C3DPoint::Trans(C3DPoint Log, C3DPoint Vector, double Angle, C3DPoint Centr) const
	{
		C3DPoint TrRes;
		double l = Vector.x;
		double m = Vector.y;
		double n = Vector.z;
		double asd = l*l + m*m + n*n;
		double asd_sq = sqrt(asd);
		double sin_ug = sin(Angle);
		double cos_ug = cos(Angle);
		C3DPoint CX;

		CX = Log - Centr;

		TrRes.x = Centr.x + ((l*l + (m*m + n*n)*cos_ug)*CX.x + (l*m - l*m*cos_ug - n*asd_sq*sin_ug)*CX.y + (l*n - l*n*cos_ug + m*asd_sq*sin_ug)*CX.z) / asd;
		TrRes.y = Centr.y + ((m*m + (l*l + n*n)*cos_ug)*CX.y + (n*m - n*m*cos_ug - l*asd_sq*sin_ug)*CX.z + (l*m - l*m*cos_ug + n*asd_sq*sin_ug)*CX.x) / asd;
		TrRes.z = Centr.z + ((n*n + (m*m + l*l)*cos_ug)*CX.z + (l*n - l*n*cos_ug - m*asd_sq*sin_ug)*CX.x + (m*n - m*n*cos_ug + l*asd_sq*sin_ug)*CX.y) / asd;

		return TrRes;
	}


	DBL C3DPoint::ToAngle() const
	{
		DBL angle = acos(x / Norm());

		if (y < 0)
			angle = M_2PI - angle;

		return angle;
	}
	//---------------------------------------------------------------------------

	/**/
	// Class C3DRect //
	/**/
	C3DRect::C3DRect(double min_X, double maxX, double min_Y, double maxY, double min_Z, double maxZ) {
		m_minX = min_X;
		m_maxX = maxX;
		m_minY = min_Y;
		m_maxY = maxY;
		m_minZ = min_Z;
		m_maxZ = maxZ;
	}
		
}