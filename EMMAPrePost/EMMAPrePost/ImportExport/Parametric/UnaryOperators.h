#pragma once

#include "Expression.h"
#include "Reader.h"
#include "Multimethods.h"

namespace MM {

//!		Унарный оператор.
/*!
		Базовый класс для всех операторов, имеющих дело с одним операндом.
*/
class CUnaryOperator : public CExpression {
public:
	CUnaryOperator () {}
	CUnaryOperator (const CUnitSptr & operand) : m_operand(operand) {}

protected:
	CUnitSptr m_operand;
};

// ------------------------------------------------------------------------------------------------

//!		Префиксный оператор.
/*!
		Базовый класс для всех унарных операторов, стоищих слева от операнда.
*/
class CPrefixOperator : public CUnaryOperator {
public:
	CPrefixOperator () { Init(); }
	CPrefixOperator (const CUnitSptr & unit) : CUnaryOperator(unit) { Init(); }
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const;
	
private:
	void Init () { m_priority = CPriority::PrefixOperator; }
};

class CNegation : public CPrefixOperator {
public:
	CNegation () { Init(); }
	CNegation (const CUnitSptr & operand) : CPrefixOperator(operand) { Init(); }
	
	CUnitSptr Result() const;
	
private:
	void Init () { m_definition = "-"; }
};

class CUnaryPlus : public CPrefixOperator {
public:
	CUnaryPlus () { Init(); }
	CUnaryPlus (const CUnitSptr & operand) : CPrefixOperator(operand) { Init(); }
	
	CUnitSptr Result() const;
	
private:
	void Init () { m_definition = "+"; }
};

// ------------------------------------------------------------------------------------------------

//!		Постфиксный оператор.
/*!
		Базовый класс для всех унарных операторов, стоищих справа от операнда.
*/
class CPostfixOperator : public CUnaryOperator {
public:
	CPostfixOperator () { Init(); }
	CPostfixOperator (const CUnitSptr & operand) : CUnaryOperator(operand) { Init(); }
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const;
	
private:
	void Init () { m_priority = CPriority::PostfixOperator; }
};


} // namespace MM