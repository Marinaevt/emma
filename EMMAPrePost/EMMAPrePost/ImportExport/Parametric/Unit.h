#pragma once

#include <iostream>
#include <string>

#include "xtl/xtl_memory.h"
#include "xtl/xtl_log.h"

namespace MM {

class CReader;

//!		Базовый класс в структуре метамодели.
/*!
		Предок всех геометрических элементов и арифметических выражений, то есть классов программы,
	которые обрабатываются интерпретатором и выстраивают дерево структуры метамодели.
*/
class CUnit {
public:
	CUnit() {
		m_definition.clear();
	}
	virtual ~CUnit () {}	

	// ---- Доступ к переменным -------------------------------------------------------------------
	
	const std::string & Definition () const { return m_definition; }
	
	virtual xtd::smart_ptr<CUnit> Result () const = 0;
	
	// ---- Индикаторы ----------------------------------------------------------------------------
	
	virtual bool DefinitionIs (const std::string & string) const { return m_definition == string; }
	
	// ---- Ввод-вывод ----------------------------------------------------------------------------
	
	virtual bool Load (CReader & reader) = 0;
	virtual void Save (std::ostream & stream) const = 0;
	virtual void SaveAsExpression (std::ostream & stream) const = 0;

private:
	//! Инициалицация.
	/*!
		Эта функция обязана быть определена в каждом неабстрактном потомке и вызываться в каждом его конструкторе.
		В ней указываются уникальный идентификатор каждого конкретного элемента.
	*/
	virtual void Init () = 0;
	
protected:
	std::string m_definition; //! Уникальный идентификатор элемента.
};

typedef xtd::smart_ptr<CUnit> CUnitSptr;

} // namespace MM