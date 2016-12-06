#pragma once

#include <list>
#include "xtl/xtl_log.h"
#include "Compare.h"
#include "Expression.h"
#include "Reader.h"

namespace MM {

//!		Мультиоперандное выражение.
/*!
		Является базовым классом для выражений, которые использую для вычислений несколько операндов.
	Синтаксис таких выражений сходен с синтаксисом функций нескольких переменных:
		"название_функции(операнд, операнд, ...)"
*/
class CMultioperandExpression : public CExpression {
public:
	CMultioperandExpression () {}
	
	CUnitSptr Result () const;
	
	size_t OperandsCount () const { return m_operands.size(); }
	
	std::list<CUnitSptr>::const_iterator OperandsBegin () const { return m_operands.begin(); }
	std::list<CUnitSptr>::const_iterator   OperandsEnd () const { return   m_operands.end(); }
	
protected:
	std::list<CUnitSptr> m_operands;
};

//!		Открывающая скобка.
class COpeningBracket : public CMultioperandExpression {
public:
	COpeningBracket () { Init(); }
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const;

private:
	void Init ();
};

//!		Закрывающая скобка.
class CClosingBracket : public CExpression {
public:
	CClosingBracket () { Init(); }

	CUnitSptr Result () const { return nullptr; }

	bool Load (CReader & reader);
	void Save (std::ostream & stream) const { xtd::warning::report("Закрывающей скобке нечего сохранять. Это служебный символ."); }

private:
	void Init ();
};

//!		Разделитель.
/*
		Является дополнением к скобкам, делая возможной загрузку нескольких операндов.
*/
class CComma : public CExpression {
public:
	CComma () { Init(); }
	
	CUnitSptr Result () const { return nullptr; }
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const { xtd::warning::report("Запятой нечего сохранять. Это служебный символ."); }
	
private:
	void Init ();
};

// ------------------------------------------------------------------------------------------------

////!		Комментарий в стиле C++.
//class CComment : public CExpression {
//public:
//	CComment () { Init(); }
//	
//	CUnitSptr Result () const { return nullptr; }
//	
//	bool Load (CReader & reader);
//	void Save (std::ostream & stream) const { xtd::warning::report("Комментарию нечего сохранять. Это служебный символ."); }
//
//private:
//	void Init () { m_definition = "//"; }
//};

//!		Комментарий в стиле C.
class CCommentStart : public CExpression {
public:
	CCommentStart () { Init(); }
	
	CUnitSptr Result () const { return nullptr; }
	
	bool Load (CReader & reader);
	void Save (std::ostream & stream) const { xtd::warning::report("Комментарию нечего сохранять. Это служебный символ."); }

private:
	void Init () { m_definition = "/*"; }
};

class CCommentEnd : public CExpression {
public:
	CCommentEnd () { Init(); }
	
	CUnitSptr Result () const { return nullptr; }
	
	bool Load (CReader & reader) { return xtd::error::report("Найден символ окончания комментария, хотя начала не было."); }
	void Save (std::ostream & stream) const { xtd::warning::report("Комментарию нечего сохранять. Это служебный символ."); }

private:
	void Init () { m_definition = "*/"; }
};

//!		Пробелы.
class CSpace : public CExpression {
public:
	CSpace () { Init(); }
	
	CUnitSptr Result () const { return nullptr; }
	
	bool DefinitionIs (const std::string & string) const { return xtd::contains_all_of(m_definition, string); }
	
	bool Load (CReader & reader) { 
		return xtd::success::report(); 
	}
	void Save (std::ostream & stream) const { xtd::warning::report("Пробелу нечего сохранять. Это служебный символ."); }

private:
	void Init () { m_definition = " \t\n\v\f\r"; }
};

} // namespace MM