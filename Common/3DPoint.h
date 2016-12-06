#pragma once
#include "AllTypes.h"
#include "2DPoint.h"	//для конструктора 3D точки
#include "SimpleMath.h"	//константы вида EPS, M_PI ...


namespace Math {
	
	class C3DPoint
	{
	public:
		union
		{
			DBL c[3];
			struct { DBL  x, y, z; };
			struct { DBL m_x, m_y, m_z; };
			struct { DBL _0, _1, _2; };
		};

		static const C3DPoint Ortx;
		static const C3DPoint Orty;
		static const C3DPoint Ortz;
		static const C3DPoint Zero;

		C3DPoint() {}
		C3DPoint(DBL x_, DBL y_, DBL z_) : x(x_), y(y_), z(z_) {}
		C3DPoint(const C3DPoint& p) : x(p.x), y(p.y), z(p.z) {}
		C3DPoint(const DBL *p) : x(p[0]), y(p[1]), z(p[2]) {}	//TODO: danger
		C3DPoint(DBL x_, Math::C2DPoint p) : x(x_), y(p.y), z(p.x) {}
		C3DPoint(Math::C2DPoint p, DBL z_) : x(p.x), y(p.y), z(z_) {}

		C3DPoint operator + (const C3DPoint& p) const;
		C3DPoint operator - (const C3DPoint& p) const;
		C3DPoint operator *  (DBL alpha) const;
		C3DPoint operator /  (DBL alpha) const;

		friend C3DPoint operator * (DBL alpha, const C3DPoint& p);
		C3DPoint& operator = (const C3DPoint& p);

		C3DPoint& operator += (const C3DPoint& p);
		C3DPoint& operator -= (const C3DPoint& p);
		C3DPoint& operator *= (DBL alpha);
		C3DPoint& operator /= (DBL alpha);

		DBL operator * (const C3DPoint& p) const; // скалярное произведение

		bool operator == (const C3DPoint& p) const;
		bool operator != (const C3DPoint& p) const;

		operator const DBL*() const;
		operator bool() const;
		bool operator ! () const;

		DBL Norm() const;  // норма
		DBL Norm2() const; // норма ^ 2

		void Normalize();

		C3DPoint Rotate(DBL phi, bool ccw = true) const;
		C3DPoint Trans(C3DPoint Log, C3DPoint Vector, double Angle, C3DPoint Centr) const;


		C3DPoint Min(const C3DPoint& p) const; // возвращает вектор с минимальными  координатами
		C3DPoint Max(const C3DPoint& p) const; // возвращает вектор с максимальными координатами

		C3DPoint Lerp(const C3DPoint& p, DBL t) const;

		C3DPoint& FromAngle(DBL angle);
		DBL ToAngle() const;

		~C3DPoint() {}
	};

	//Класс параллелепипеда в 3D
	//TODO: почти всё
	class C3DRect {
	public:
		C3DRect(double m_X = 0, double mX = 0, double m_Y = 0, double mY = 0, double m_Z = 0, double mZ = 0);
		~C3DRect(void) {};
	public:
		double m_minX;
		double m_maxX;
		double m_minY;
		double m_maxY;
		double m_minZ;
		double m_maxZ;
	public:
		double dX() { return fabs(m_maxX - m_minX); }
		double dY() { return fabs(m_maxY - m_minY); }
		double dZ() { return fabs(m_maxZ - m_minZ); }
	};
}

//Сокращённые названия
typedef Math::C3DPoint CPoint3D;
typedef Math::C3DRect CRect3D;
