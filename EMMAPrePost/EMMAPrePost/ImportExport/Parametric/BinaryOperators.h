#pragma once

#include "Expression.h"
#include "Reader.h"
#include "Multimethods.h"
#include "Operands.h"

namespace MM {

//!		Бинарный оператор.
/*!
		Базовый класс для всех операторов, имеющих дело с двумя операндами.
*/
class CBinaryOperator : public CExpression {
public:
	CBinaryOperator () {}
	CBinaryOperator (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) :
		m_leftOperand(leftOperand), m_rightOperand(rightOperand) {}
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const;

protected:
	CUnitSptr m_leftOperand;
	CUnitSptr m_rightOperand;
};

// ------------------------------------------------------------------------------------------------

//!		Аддитивный оператор.
/*!
		Базовый класс для всех аддитивных операторов, таких как сложение и вычитание.
*/
class CAdditiveOperator : public CBinaryOperator {
public:
	CAdditiveOperator () { Init(); }
	CAdditiveOperator (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) :
		CBinaryOperator(leftOperand, rightOperand) { Init(); }
		
private:
	void Init () { m_priority = CPriority::AdditiveOperator; }
};

//!		Сумма.
class CSum : public CAdditiveOperator {
public:
	CSum () { Init(); }
	CSum (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) : CAdditiveOperator(leftOperand, rightOperand) { Init(); }

	CUnitSptr Result() const;

private:
	void Init() { m_definition = "+"; }
};

//!		Разность.
class CDifference : public CAdditiveOperator {
public:
	CDifference () { Init(); }
	CDifference (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) : CAdditiveOperator(leftOperand, rightOperand) { Init(); }

	CUnitSptr Result() const;

private:
	void Init() { m_definition = "-"; }
};


// ------------------------------------------------------------------------------------------------

//!		Мультипликативный оператор.
/*!
		Базовый класс для всех мультипликативных операторов, таких как умножение и деление.
*/
class CMultiplicativeOperator : public CBinaryOperator {
public:
	CMultiplicativeOperator () { Init(); }
	CMultiplicativeOperator (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) :
		CBinaryOperator(leftOperand, rightOperand) { Init(); }
		
private:
	void Init () { m_priority = CPriority::MultiplicativeOperator; }
};

//!		Произведение.
class CProduct : public CMultiplicativeOperator {
public:
	CProduct () { Init(); }
	CProduct (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) :
		CMultiplicativeOperator(leftOperand, rightOperand) { Init(); }

	CUnitSptr Result () const { return m_leftOperand->Result() * m_rightOperand->Result(); }

private:
	void Init() { m_definition = "*"; }
};

//!		Частное.
class CQuotient : public CMultiplicativeOperator {
public:
	CQuotient () { Init(); }
	CQuotient (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) :
		CMultiplicativeOperator(leftOperand, rightOperand) { Init(); }

	CUnitSptr Result () const { return m_leftOperand->Result() / m_rightOperand->Result(); }

private:
	void Init() { m_definition = "/"; }
};

// ------------------------------------------------------------------------------------------------

//!		Степень.
class CPower : public CBinaryOperator {
public:
	CPower () { Init(); }
	CPower (const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) :
		CBinaryOperator(leftOperand, rightOperand) { Init(); }
	
	CUnitSptr Result () const { return Power(m_leftOperand->Result(), m_rightOperand->Result()); }
	
private:
	void Init ();
};

} // namespace MM