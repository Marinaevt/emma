#pragma once
#include "AllTypes.h"
#include "SimpleMath.h"
#include <vector>
#include "DlgShowError.h"
#include "ResourceCommon.h"

/*
Содержит классы:
C2DPoint - 2D точка
C2DRect - 2D прямоугольник
C2DLine - 2D линия (Ax + By + C = 0)
*/

namespace Math {

	/**!
	Класс C2DPoint реализует функционал двумерной точки, имеющей две координаты: (X ; Y).
	*/
	class C2DPoint
	{
	public:
		union //без названия
		{
			//именовать структуры нежелательно (много переделывать)
			DBL coord[2];
			struct { DBL x, y; };
			struct { DBL X, Y; };
			struct { DBL _0, _1; };
			struct { DBL m_x, m_y; };
			struct { DBL m_X, m_Y; };
		};

		static const C2DPoint Ortx;
		static const C2DPoint Orty;
		static const C2DPoint Zero;

	public:
		C2DPoint() : x(0), y(0) {}
		C2DPoint(DBL x_, DBL y_) : x(x_), y(y_) {}
		C2DPoint(const C2DPoint& p) : x(p.x), y(p.y) {}
		C2DPoint(const DBL *p) : x(p[0]), y(p[1]) {}

		C2DPoint operator + (const C2DPoint& p) const;
		C2DPoint operator - (const C2DPoint& p) const;
		C2DPoint operator *  (DBL alpha) const;
		C2DPoint operator /  (DBL alpha) const;

		friend C2DPoint operator * (DBL alpha, const C2DPoint& p);

		C2DPoint& operator = (const C2DPoint& p);
		C2DPoint& operator += (const C2DPoint& p);
		C2DPoint& operator -= (const C2DPoint& p);
		C2DPoint& operator *= (DBL alpha);
		C2DPoint& operator /= (DBL alpha);

		// скалярное произведение
		DBL operator * (const C2DPoint& p) const;

		bool operator == (const C2DPoint& p) const;
		bool operator != (const C2DPoint& p) const;

		operator const DBL*() const;
		operator bool() const;
		bool operator ! () const;

		//Норма (евклидова)
		DBL Norm() const;
		// Норма ^ 2 (евклидова)
		DBL Norm2() const;
		// Норма (L1, манхэттенское расстояние) - сумма координат 
		DBL Norm_L1() const;

		// Нормализация вектора
		void Normalize();

		// Поворот вектора на угол PHI по или против часовой стрелки
		C2DPoint Rotate(DBL phi, bool ccw = true) const;

		// Возвращает вектор с минимальными координатами (из обоих координат)
		C2DPoint Min(const C2DPoint& p) const;
		// Возвращает вектор с максимальными координатами (из обоих координат)
		C2DPoint Max(const C2DPoint& p) const;

		// Линейная интерполяция вектора по t (начало - в основной точке)
		C2DPoint Lerp(const C2DPoint& p, DBL t) const;

		// Задаёт координаты вида (cos(x); sin(x))
		C2DPoint& FromAngle(DBL angle);
		// Возвращает угол по (x; y) и arccos() диапазон [0; 2*PI]
		DBL ToAngle() const;

		bool Right(const C2DPoint& p1, const C2DPoint& p2) const;

		//расстояние до точки
		double Len(const C2DPoint &p) const;
		double LenX(const C2DPoint &p) const;
		//Расстояние до отрезка [p1;p2] и ближайшая точка на нём
		double GetDistance(const C2DPoint &DP1, const C2DPoint &DP2, C2DPoint *pDPRes) const;
		//Расщепление
		double Splitting(const C2DPoint &p1, const C2DPoint &p2);
		//объект устанавливается в "третью" точку (достраивая равносторонний треугольник)
		bool Third(const C2DPoint &p1, const C2DPoint &p2);
		
		//объект устанавливается на прямую, соединяющую p1 и p2 на расстоянии len от p1 в направлении p2
		bool Rad(const C2DPoint &p1, const C2DPoint &p2, double len);
		//пересечение отрезков
		bool Cross(const C2DPoint &pa1, const C2DPoint &pa2, const C2DPoint &pb1, const C2DPoint &pb2);
		//слишком близко
		bool IfTooClose(const C2DPoint &p1, const C2DPoint &p2);
		///////////////////////////////////////////////////////////////////////////////////

		//Проверка на нулевую точку
		bool IsZero() const { return fabs(X) < EPS && fabs(Y) < EPS; };
		static double CrossProduct(const C2DPoint& p1, const C2DPoint& p2) { return p1.X * p2.Y - p1.Y * p2.X; }
		//Вычисление хэша узла (наследие .Net и C# для быстрого определения элемента коллекции)
		static int GetHash(const C2DPoint& p);

		~C2DPoint() {}
	};

	// Для векторов

	//Аппроксимация 
	double APPR_ex(std::vector<C2DPoint> &XY, double x);

	//Интервал [0;1] в отрезок [A;B]
	inline C2DPoint i01toAB(const C2DPoint& A, const C2DPoint& B, DBL x);

	// Нахождение координат центра круга по 2м точкам и радиусу
	bool CircleCenter(const C2DPoint& p1, const C2DPoint& p2, DBL r, C2DPoint& c);
	// Нахождение координат центра круга по 3м точкам (НЕ ТЕСТИРОВАЛОСЬ)
	bool CircleCenter(const C2DPoint& p1, const C2DPoint& p2, const C2DPoint& p3, C2DPoint& c);
	bool IsPerpendicular(const C2DPoint& pt1, const  C2DPoint& pt2, const C2DPoint& pt3);
	bool CalcCircleCenter(const C2DPoint& pt1, const C2DPoint& pt2, const C2DPoint& pt3, C2DPoint& center);

	//Ближайшая точка на отрезке p1p2 к точке V (на отрезке)
	C2DPoint ClosestPoint(const C2DPoint& v, const C2DPoint& p1, const C2DPoint& p2);
	//Ближайшая точка на отрезке p1p2 к точке V (на дуге)
	C2DPoint ClosestPoint(const C2DPoint& v, const C2DPoint& p1, const C2DPoint& p2, DBL r);

	DBL ToAngle(const C2DPoint& p);
	C2DPoint FromAngle(DBL angle);

	/**!
	Класс C2DRect реализует функционал прямоугольника.
	Он же CRect2D.
	*/
	class C2DRect
	{
	public:
		C2DRect(DBL l = 0.0, DBL r = 0.0, DBL b = 0.0, DBL t = 0.0);

		C2DRect(const Math::C2DPoint& p0, const Math::C2DPoint& p1);

		C2DRect(const std::initializer_list<Math::C2DPoint>& points) {
			for (auto point : points) Union(point); 
		}

		~C2DRect(void) {};

	public:

		union //без названия
		{
			//именовать структуры нежелательно (много переделывать)
			struct { DBL left, top, right, bottom; };
			struct { DBL MinX, MinY, MaxX, MaxY; };
		};

	public:
		DBL Width() const { return right - left; }
		DBL Heigth() const { return top - bottom; }
		
		//BoundingBox(const std::vector<CPoint2D>& points) { for (auto point : points) Union(point); }
		bool IntersectsWith(const C2DRect& bbox) const { 
			return !(MaxX < bbox.MinX || MinX > bbox.MaxX || MaxY < bbox.MinY || MinY > bbox.MaxY); 
		}

		void Fix();

		//Левая верхняя точка
		C2DPoint GetLTPoint() const
		{
			return C2DPoint(left, top);
		}

		//Правая нижняя точка
		C2DPoint GetRBPoint() const
		{
			return C2DPoint(right, bottom);
		}

		// находит объемлющий прямоугольник
		void AddRect(const C2DRect& rect);

		// Добавляет точку и расширяет прямоугольник
		void AddPoint(const C2DPoint& p);

		// Добавляет точку и расширяет прямоугольник
		void Union(const C2DPoint& p);

		//Превращается в точку
		void SetPoint(const C2DPoint& p);

	};

	/**!
	Класс C2DLine реализует функционал двумерной прямой, имеющей направление.
	Уравнение прямой: Ax + By + C = 0.
	Параметры прямой хранятся в двух видах:
	1. Нормаль m_n и число m_c. (A = m_n.x, B = m_n.y, C = m_c);
	2. Точки m_p1 и m_p2, через которые проходит прямая.
	*/
	class C2DLine
	{
	protected:
		C2DPoint m_n; // нормаль или коэффициенты a, b уравнения прямой
		DBL m_c; // коэффициент c уравнения прямой		

		void CalcNormC();

		void CalcPoints(bool dir = true);

	public:

		//координаты точек отрезка
		union {
			struct { C2DPoint m_p1, m_p2; }; // две точки, через которые проходит прямая
			struct { C2DPoint S, E; };		//start, end
		};

		C2DLine(): m_n(0, 1), m_c(0)
		{
			CalcPoints(true);
		}
		
		//конструктор по нормали
		C2DLine(const C2DPoint& n, DBL c);
		//конструктор по 2м точкам
		C2DLine(const C2DPoint& p1, const C2DPoint& p2);
		//конструктор по 4м координатам
		C2DLine(double s1, double e1, double s2, double e2);
		//конструктор по вектору
		C2DLine(const C2DLine& p) : S(p.S), E(p.E) { CalcNormC(); }

		~C2DLine(void) {};

		C2DPoint Vector() const { return E - S; }
		double Length() const { return Vector().Norm(); }
		// точка лежит внутри отрезка
		bool Inside(const C2DPoint& p) const;
		
		C2DRect GetBoundingBox() const { return C2DRect({ S, E }); }
		
		bool operator == (const C2DLine& p) const;
		bool operator != (const C2DLine& p) const;
		C2DLine& operator = (const C2DLine& p);

		// хеш для работы std::unordered_set
		static int GetHash(const C2DLine& p)
		{
			return C2DPoint::GetHash(p.S) ^ C2DPoint::GetHash(p.E);
		}

		///! инициализирует прямую по двум точкам
		void SetPoints(const C2DPoint& p1, const C2DPoint& p2);

		///! инициализирует прямую по нормали, числу и направлению (направление нужно, чтобы определить порядок точек m_p1 и m_p2)
		void SetNormC(const C2DPoint& n, const DBL& c, bool dir = true);

		const C2DPoint& GetP1() const
		{
			return m_p1;
		}

		const C2DPoint& GetP2() const
		{
			return m_p2;
		}

		const C2DPoint& GetN() const
		{
			return m_n;
		}

		const DBL& GetC()
		{
			return m_c;
		}

		DBL GetY(DBL x);

		DBL GetX(DBL y);

		//! возвращает точку пересечения с прямой l, если таковая существует
		//! -2 - прямые совпадают, -1 - прямые параллельны, 0 - прямые пересекаются
		int GetIntersection(const C2DLine& l, C2DPoint& res);

		//! проверяет, находится ли точка p слева от прямой (слева - это там, куда смотрит нормаль)
		bool IsLeft(const C2DPoint& p)
		{
			return (p.x * m_n._0 + p.y * m_n._1 + m_c) > 0;
		}

		//! проверяет, находится ли точка p на прямой с точностью EPS
		bool IsOnLine(const C2DPoint& p)
		{
			return fabs(m_n.x * p.x + m_n.y * p.y + m_c) < EPS;
		}

		///! возвращает прямую, перпендикулярную данной в точке p
		C2DLine GetPerpendicular(const C2DPoint& p) const;

	};

}; // namespace Math

class  ContourGeometry
{
public: std::vector<Math::C2DPoint> ContourPoints;
		// std::vector<ContourGeometry*> Holes; пока не используется
};

   //Сокращённые названия
typedef Math::C2DPoint CPoint2D;

typedef Math::C2DRect CRect2D;

typedef Math::C2DLine CLine2D;