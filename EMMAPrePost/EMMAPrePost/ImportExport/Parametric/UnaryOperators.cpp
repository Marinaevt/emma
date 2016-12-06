#include "StdAfx.h"
#include "UnaryOperators.h"

#include <iostream>

/*
	Унарные операторы имеют бОльший приоритет по сравнению с бинарными!
	Это связано с загрузкой meta файлов (координаты + и -)
	Иначе при координате -1 бинарный оператор не найдёт левого числа и свалится с ошибкой!
	Знак минус считается унарным.
*/

namespace MM {

bool CPrefixOperator::Load (CReader & reader) {
	// Если непосредственно перед префиксным оператором стоит число или закрывающая скобка,
	// то оно сообщает о невозможности загрузки. В будущем это можно будет переработать таким образом,
	// чтобы вместо сообщения об ошибке туда автоматически вставлялся оператор умножения.
//	if (!reader.HasReadyUnit()) {
		reader.AddWaitingUnit(this);
		
		if (reader.Read() == xtd::good) {
			reader.ExtractWaitingUnit();
		
			m_operand = reader.ExtractReadyUnit();
			reader.AddReadyUnit(this);
			
			return xtd::success::report();
		}
		else {
			return xtd::error::report("Нет операнда для оператора \"" + m_definition + "\".");
		}
	//}
	//else {
//		return xtd::error::report("Оператор \"" + m_definition + "\" не может быть загружен.");
//	}

	// В данный момент сообщение о невозможности загрузки отключено для правильной работы операторов
	// внутри скобок.
}

void CPrefixOperator::Save (std::ostream & stream) const {
	stream << m_definition;
	m_operand->SaveAsExpression(stream);
}

CUnitSptr CNegation::Result() const {
	CUnitSptr t = -m_operand->Result();
	return t;
}

CUnitSptr CUnaryPlus::Result() const {
	CUnitSptr t = m_operand->Result();
	return t;
}

bool CPostfixOperator::Load (CReader & reader) {
	if (reader.HasReadyUnit()) {
		m_operand = reader.ExtractReadyUnit();
		reader.AddReadyUnit(this);
		
		return xtd::success::report();
	}
	
	//иначе
	return xtd::error::report("Нет операнда для оператора \"" + m_definition + "\".");
	
}

void CPostfixOperator::Save (std::ostream & stream) const {
	m_operand->Save(stream);
	stream << m_definition;
}

} // namespace MM