#ifndef _XTL_DOUBLE_DISPATCHER_H
#define _XTL_DOUBLE_DISPATCHER_H

#include <typeinfo>
#include <map>

#include "xtl_memory.h"
#include "xtl_functional.h"

namespace xtd {

template <
	typename base_type,
	template <
		typename first_argument_type,
		typename second_argument_type,
		typename result_type
	> class binary_predicate
>
class double_dispatcher {
private:
	// ---- Определения ---------------------------------------------------------------------------
	
	typedef smart_ptr<base_type> base_type_sptr;
	typedef binary_function<base_type_sptr> base_type_binary_function;
	typedef smart_ptr<base_type_binary_function> base_type_binary_function_sptr;

	//!		Исполнитель операции "binary_predicate" над двумя наследниками класса "base_type".
	/*!
			Является промежуточным звеном между диспетчером и соответствующей предикату "binary_predicate" операцией,
		которая может быть как методом класса "one_child_type", принимающим на вход экземпляр класса "another_child_type",
		так и независимой функцией этих двух переменных (классы "one_child_type" и "another_child_type" должны быть
		наследниками класса "base_type"). Результат этой операции используется для создания объекта класса "result_type".
		\return Указатель на вновь созданный объект класса "result_type", который автоматически преобразуется к умному указателю
		на класс "base_type".
	*/
	template <
		typename one_child_type,
		typename another_child_type,
		typename result_type
	>
	class binary_executor : public base_type_binary_function {
	public:
		base_type_sptr operator () (const base_type_sptr & left_operand, const base_type_sptr & right_operand) const {
			binary_predicate<one_child_type, another_child_type, result_type> predicate;
			
			return new result_type(
				predicate(
					dynamic_cast<one_child_type &>(*left_operand),
					dynamic_cast<another_child_type &>(*right_operand)
				)
			);
		}
	};

public:	
	// ---- Модификация ---------------------------------------------------------------------------
	
	//!		Регистрация бинарного функтора.
	/*!
			Записывает в таблицу новую пару идентификаторов, полученных от классов "first_type" и "second_type",
		и соответствующий этой паре функтор "binary_executor", который при вызове выполнит над объектами типов
		"first_type" и "second_type" операцию "binary_predicate", создаст из результата новый экземпляр класса "result_type"
		и вернёт его в виде указателя на класс "base_type".
			Если операция "binary_predicate" над типами "first_type" и "second_type" не определена, то программа не скомпилируется.
	*/
	template <typename first_type, typename second_type, typename result_type>
	void put () {	
		m_dispatch_table.insert(
			std::make_pair(
				key_type(
					typeid(first_type).name(),
					typeid(second_type).name()
				),
				new binary_executor<first_type, second_type, result_type>()
			)
		);
	}
	
	//!		Упрощённый шаблон регистрации бинарного функтора.
	template <typename type>
	void put () { put<type, type, type>(); }
	
	// ---- Доступ --------------------------------------------------------------------------------
	
	//!		Диспетчеризация.
	/*!
			Ищет в таблице пару идентификаторов объектов "left_operand" и "right_operand", и если находит,
		то выполняет соответствующую (заранее зарегистрированную) операцию над этими объектами и возвращает результат.
		В противном случае возвращает нулевой указатель ("nullptr").
	*/
	base_type_sptr dispatch (const base_type_sptr & left_operand, const base_type_sptr & right_operand) const {
		typename dispatch_table_type::const_iterator i;
		
		i = m_dispatch_table.find(
			key_type(
				typeid(*left_operand).name(),
				typeid(*right_operand).name()
			)
		);
		
		if (i != m_dispatch_table.end()) return (*i->second)(left_operand, right_operand);
		else return nullptr;
	}
	
private:
	typedef std::pair<std::string, std::string> key_type;
	typedef std::map<key_type, base_type_binary_function_sptr> dispatch_table_type;
	dispatch_table_type m_dispatch_table;
};

}

#endif