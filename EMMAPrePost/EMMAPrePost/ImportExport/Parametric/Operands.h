#pragma once

#include <iostream>

#include "xtl/xtl_string.h"

#include "Expression.h"
#include "Reader.h"

namespace MM {

//!		Число.
/*!
		Простейший операнд.
*/
class CNumber : public CExpression {
public:
	CNumber () { 
		m_value = 0.0; 
		Init(); 
	}
	CNumber (double value) : m_value(value) { Init(); }
	
	// ---- Доступ к переменным -------------------------------------------------------------------
	
	double Value () const { return m_value; }
	
	CUnitSptr Result () const { return new CNumber(m_value); }
	
	// ---- Индикаторы ----------------------------------------------------------------------------
	
	bool DefinitionIs (const std::string & string) const { return xtd::contains_all_of(m_definition, string); }
	
	// ---- Ввод-вывод ----------------------------------------------------------------------------
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const { stream << m_value; }

private:
	void Init ();

	double m_value;
};

//!		Симольная константа.
/*!
		Имеет имя и фиксированное значение.
		Например, числа пи и е.
*/
class CCharacterLiteral : public CExpression {
public:	
	// ---- Доступ к переменным -------------------------------------------------------------------
	
	CUnitSptr Result () const { return new CNumber(m_value); }
	
	// ---- Ввод-вывод ----------------------------------------------------------------------------
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const { stream << m_definition; }
	
protected:
	double m_value;
};

//!		Число пи.
class CPi : public CCharacterLiteral {
public:
	CPi () { Init(); }
	
private:
	void Init ();
};

//!		Число е.
class CExponent : public CCharacterLiteral {
public:
	CExponent () { Init(); }
	
private:
	void Init ();
};

} // namespace MM