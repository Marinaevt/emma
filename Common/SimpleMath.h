#pragma once
#include "AllTypes.h"
#include "math.h"

//Погрешность измерений  10^(-15)
#define EPS	     0.000000000000001
//const double EPS = 0.000000000000001;	//говорят const <type> лучше для разных компиляторов
//Погрешность измерений  10^(-10)
#define EPS_10	 0.0000000001
//Погрешность измерений  10^(-3)
#define EPS_3	 0.001

//Двадцать (20) знаков после запятой

#define M_PI     3.14159265358979323846
#define M_2PI    6.28318530717958647692
#define M_PI2    1.57079632679489661923
#define M_1_2PI  0.15915494309189533577
#define M_1_180  0.00555555555555555556
#define M_PI_180 0.01745329251994329577
#define M_180_PI 57.2957795130823208768
#define M_E		 2.71828182845904523536

#define LENGTHUNIT 0.1 //Число чуть меньше 1, для кривых, контуров и геометрии

namespace Math {

// Интерполяционные процедуры
// Если x задается не в пределах [A, B], то
// эти процедуры также являются линейной
// экстраполяцией.

// [A, B] -> [0, 1], если B > A и A <= x <= B
// [A, B] -> [-1, 0], если B < A и B<= x <= A
// B != A !!!
inline DBL iABto01(DBL A, DBL B, DBL x)
{
	if (fabs(B - A) > EPS) {
		return (x - A) / (B - A);
	}
	
	return 0.0;
}

// [0, 1] -> [A, B]
inline DBL i01toAB(DBL A, DBL B, DBL x)
{
	if (x >= 0.0 && x <= 1.0) {
		return A + (B - A) * x;
	}

	return 0.0;
}

// if > 0 -> +1, if < 0 -> -1, if 0 -> 0
inline DBL sgn(DBL x)
{
    return ((x) > 0.0 ? 1 : ((x) < 0.0 ? -1 : 0));
}

//Гиперболический Арксинус arcsinh()
inline DBL arcsinh(double x) {
	return log(x + sqrt(x*x + 1));
}

double APPR(double z1, double w1, double z2, double w2, double zn);


/*
template <typename T>
inline T min(const T& a, const T& b)
{
	if ((a) < (b)) return (a);
	else return (b);
}

template <class T> inline T max(const T& a, const T& b)
{
	if ((a) > (b)) return (a);
	else return (b);
}
*/

template <typename T>
void swap(T& a, T& b)
{
    T tmp = a;
    a = b;
    b = tmp;
}

template <typename T>
void sortUp(T& a, T& b)
{
    if (a > b) swap(a, b);
}

template <typename T>
void sortDown(T& a, T& b)
{
    if (a < b) swap(a, b);
}

template <typename T>
bool CheckBounds(T a, T b, T x)
{
	if ((x > a) && (x < b)) return true;
	return false;
}

template <typename T>
bool CheckBounds(T a, T b, T x1, T x2)
{
	return (CheckBounds(a, b, x1) && CheckBounds(a, b, x2));
}

template <typename T>
bool CheckBounds(T a, T b, T x1, T x2, T x3)
{
	return (CheckBounds(a, b, x1, x2) && CheckBounds(a, b, x3));
}

template <typename T>
bool CheckBoundsInc(T a, T b, T x)
{
	if ((x < a) || (x > b)) return false;
	return true;
}

template <typename T>
bool CheckBoundsInc(T a, T b, T x1, T x2)
{
	return (CheckBoundsInc(a, b, x1) && CheckBoundsInc(a, b, x2));
}

template <typename T>
bool CheckBoundsInc(T a, T b, T x1, T x2, T x3)
{
	return (CheckBoundsInc(a, b, x1, x2) && CheckBoundsInc(a, b, x3));
}

template <typename T>
bool CheckBoundsLI(T a, T b, T x)
{
	if ((x >= a) && (x < b)) return true;
	return false;
}

template <typename T>
bool CheckBoundsLI(T a, T b, T x1, T x2)
{
	return (CheckBoundsLI(a, b, x1) && CheckBoundsLI(a, b, x2));
}

template <typename T>
bool CheckBoundsLI(T a, T b, T x1, T x2, T x3)
{
	return (CheckBoundsLI(a, b, x1, x2) && CheckBoundsLI(a, b, x3));
}

template <typename T>
bool CheckBoundsRI(T a, T b, T x)
{
	if ((x > a) && (x <= b)) return true;
	return false;
}

template <typename T>
bool CheckBoundsRI(T a, T b, T x1, T x2)
{
	return (CheckBoundsRI(a, b, x1) && CheckBoundsRI(a, b, x2));
}

template <typename T>
bool CheckBoundsRI(T a, T b, T x1, T x2, T x3)
{
	return (CheckBoundsRI(a, b, x1, x2) && CheckBoundsRI(a, b, x3));
}

// -----------------------------------------------------------------------

template<int M>
int UMod(int x)
{
	if (x <  0) return x + M;
	if (x >= M) return x - M;

	return x;
}

template<int M>
int UXMod(int x)
{
	while (x < 0) x += M;
	while (x >= M) x -= M;

	return x;
}

} // namespace Math