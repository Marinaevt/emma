#pragma once

#include <cmath>

#include "Element.h"
#include "UnaryOperators.h"
#include "Multimethods.h"

namespace MM {

//!		Базовый класс для операторов-функций.
class CFunction : public CPrefixOperator {
public:
	CFunction () {}
	CFunction (const CUnitSptr & operand) : CPrefixOperator(operand) {}

private:
	void Init () { m_priority = CPriority::Function; }
};

//!		Синус.
class CSine : public CFunction {
public:
	CSine () { Init(); }
	CSine (const CUnitSptr & operand) : CFunction(operand) { Init(); }
	
	CUnitSptr Result () const { return Sine(m_operand->Result()); }
	
private:
	void Init () { m_definition = "sin"; }
};

//!		Косинус.
class CCosine : public CFunction {
public:
	CCosine () { Init(); }
	CCosine (const CUnitSptr & operand) : CFunction(operand) { Init(); }
	
	CUnitSptr Result () const { return Cosine(m_operand->Result()); }
	
private:
	void Init () { m_definition = "cos"; }
};

//!		Логарифм по основанию e.
class CLogarithm : public CFunction {
public:
	CLogarithm () { Init(); }
	CLogarithm (const CUnitSptr & operand) : CFunction(operand) { Init(); }
	
	CUnitSptr Result () const { return Logarithm(m_operand->Result()); }
	
private:
	void Init () { m_definition = "ln"; }
};

} // namespace MM