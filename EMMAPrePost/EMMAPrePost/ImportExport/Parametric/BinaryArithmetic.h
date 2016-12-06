#pragma once

#include <cmath>

#include "Operands.h"
#include "Parameter.h"
#include "Point.h"

namespace MM {

// ---- Операции с двумя числами ------------------------------------------------------------------

inline CNumber operator + (const CNumber & n1, const CNumber & n2) { return n1.Value() + n2.Value(); }
inline CNumber operator - (const CNumber & n1, const CNumber & n2) { return n1.Value() - n2.Value(); }
inline CNumber operator * (const CNumber & n1, const CNumber & n2) { return n1.Value() * n2.Value(); }
inline CNumber operator / (const CNumber & n1, const CNumber & n2) { return n1.Value() / n2.Value(); }

inline CNumber pow (const CNumber & n1, const CNumber & n2) { return ::pow(n1.Value(), n2.Value()); }

// ---- Операции с двумя параметрами --------------------------------------------------------------

inline CNumber operator + (const CParameter & p1, const CParameter & p2) { return p1.Value() + p2.Value(); }
inline CNumber operator - (const CParameter & p1, const CParameter & p2) { return p1.Value() - p2.Value(); }
inline CNumber operator * (const CParameter & p1, const CParameter & p2) { return p1.Value() * p2.Value(); }
inline CNumber operator / (const CParameter & p1, const CParameter & p2) { return p1.Value() / p2.Value(); }

inline CNumber pow (const CParameter & p1, const CParameter & p2) { return pow(p1.Value(), p2.Value()); }

// ---- Операции с двумя точками ------------------------------------------------------------------

inline CPoint operator + (const CPoint & p1, const CPoint & p2) { return CPoint(p1.X() + p2.X(), p1.Y() + p2.Y()); }
inline CPoint operator - (const CPoint & p1, const CPoint & p2) { return CPoint(p1.X() - p2.X(), p1.Y() - p2.Y()); }

// ---- Операции с параметром и числом ------------------------------------------------------------

inline CNumber operator + (const CParameter & parameter, const CNumber & number) { return parameter.Value() + number; }
inline CNumber operator - (const CParameter & parameter, const CNumber & number) { return parameter.Value() - number; }
inline CNumber operator * (const CParameter & parameter, const CNumber & number) { return parameter.Value() * number; }
inline CNumber operator / (const CParameter & parameter, const CNumber & number) { return parameter.Value() / number; }

inline CNumber pow (const CParameter & parameter, const CNumber & number) { return pow(parameter.Value(), number); }

inline CNumber operator + (const CNumber & number, const CParameter & parameter) { return number + parameter.Value(); }
inline CNumber operator - (const CNumber & number, const CParameter & parameter) { return number - parameter.Value(); }
inline CNumber operator * (const CNumber & number, const CParameter & parameter) { return number * parameter.Value(); }
inline CNumber operator / (const CNumber & number, const CParameter & parameter) { return number / parameter.Value(); }

inline CNumber pow (const CNumber & number, const CParameter & parameter) { return pow(number, parameter.Value()); }

// ---- Операции с точкой и числом ----------------------------------------------------------------

inline CPoint operator * (const CPoint & point, const CNumber & number) { return CPoint(point.X() * number, point.Y() * number); }
inline CPoint operator / (const CPoint & point, const CNumber & number) { return CPoint(point.X() / number, point.Y() / number); }

inline CPoint operator * (const CNumber & number, const CPoint & point) { return point * number; }

} // namespace MM