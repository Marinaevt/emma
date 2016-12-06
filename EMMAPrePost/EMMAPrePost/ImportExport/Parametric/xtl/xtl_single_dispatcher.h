#ifndef _XTL_SINGLE_DISPATCHER_H
#define _XTL_SINGLE_DISPATCHER_H

#include <typeinfo>
#include <map>

#include "xtl_memory.h"
#include "xtl_functional.h"

namespace xtd {

template <
	typename base_type,
	template <
		typename argument_type,
		typename result_type
	> class unary_predicate
>
class single_dispatcher {
private:
	// ---- Определения ---------------------------------------------------------------------------
	
	typedef smart_ptr<base_type> base_type_sptr;
	typedef unary_function<base_type_sptr> base_type_unary_function;
	typedef smart_ptr<base_type_unary_function> base_type_unary_function_sptr;

	//!		Исполнитель операции "unary_predicate" с наследником класса "base_type".
	/*!
			Является промежуточным звеном между диспетчером и соответствующей предикату "unary_predicate" операцией,
		которая может быть как методом класса "child_type", так и независимой функцией этой переменной (класс "child_type"
		должен быть наследником класса "base_type"). Результат этой операции используется для создания объекта класса "result_type".
		\return Указатель на вновь созданный объект класса "result_type", который автоматически преобразуется к умному указателю
		на класс "base_type".
	*/
	template <
		typename child_type,
		typename result_type
	>
	class unary_executor : public base_type_unary_function {
	public:
		base_type_sptr operator () (const base_type_sptr & operand) const {
			unary_predicate<child_type, result_type> predicate;
			
			return new result_type(predicate(dynamic_cast<child_type &>(*operand)));
		}
	};

public:
	// ---- Модификация ---------------------------------------------------------------------------
	
	//!		Регистрация унарного функтора.
	/*!
			Записывает в таблицу новый идентификатор, полученный от класса "type", и соответствующий ему
		функтор "unary_executor", который при вызове выполнит над объектом типа "type" операцию "unary_predicate",
		создаст из результата новый экземпляр класса "result_type" и вернёт его в виде указателя на класс "base_type".
		Если операция "unary_predicate" над типом "type" не определена, то программа не скомпилируется.
	*/
	template <typename type, typename result_type>
	void put () {
		m_dispatch_table.insert(
			std::make_pair(
				typeid(type).name(),
				new unary_executor<type, result_type>()
			)
		); 
	}
	
	//!		Упрощённый шаблон регистрации функтора.
	template <typename type>
	void put () { put<type, type>(); }
	
	// ---- Доступ --------------------------------------------------------------------------------
	
	//!		Диспетчеризация.
	/*!
			Ищет в таблице идентификатор объекта "operand" и если находит, то выполняет
		соответствующую (заранее зарегистрированную) операцию над этим объектом и возвращает результат.
		В противном случае возвращает нулевой указатель ("nullptr").
	*/
	base_type_sptr dispatch (const base_type_sptr & operand) const {
		typename dispatch_table_type::const_iterator i;
		
		i = m_dispatch_table.find(typeid(*operand).name());
		
		if (i != m_dispatch_table.end()) return (*i->second)(operand);
		else return nullptr;
	}
	
private:
	typedef std::map<std::string, base_type_unary_function_sptr> dispatch_table_type;
	dispatch_table_type m_dispatch_table;
};

}

#endif