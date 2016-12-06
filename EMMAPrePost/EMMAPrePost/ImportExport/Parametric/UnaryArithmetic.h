#pragma once

#include <cmath>

#include "Operands.h"
#include "Parameter.h"

namespace MM {

// ---- Операции с числом -------------------------------------------------------------------------

inline CNumber operator - (const CNumber & number) { return -number.Value(); }

inline CNumber sin (const CNumber & number) { return ::sin(number.Value()); }
inline CNumber cos (const CNumber & number) { return ::cos(number.Value()); }
inline CNumber log (const CNumber & number) { return ::log(number.Value()); }

// ---- Операции с параметром ---------------------------------------------------------------------

inline CNumber operator - (const CParameter & parameter) { return -parameter.Value(); }

inline CNumber sin (const CParameter & parameter) { return sin(parameter.Value()); }
inline CNumber cos (const CParameter & parameter) { return cos(parameter.Value()); }
inline CNumber log (const CParameter & parameter) { return log(parameter.Value()); }

// ---- Операции с точкой -------------------------------------------------------------------------

//inline CPoint operator - (const CPoint & point) { return CPoint(-point.X(), -point.Y()); }

} // namespace MM