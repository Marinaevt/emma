#include "StdAfx.h"
#include "BinaryOperators.h"

#include <iostream>

/*
	Унарные операторы имеют бОльший приоритет по сравнению с бинарными!
	Это связано с загрузкой meta файлов (координаты + и -)
	Иначе при координате -1 бинарный оператор не найдёт левого числа и свалится с ошибкой!
	Знак минус считается унарным.
*/

namespace MM {

bool CBinaryOperator::Load (CReader & reader) {
	if (reader.HasReadyUnit()) {
		if (!reader.HasWaitingUnit() || dynamic_cast<CExpression &>(*reader.GetWaitingUnit()).Priority() < m_priority) {
			m_leftOperand = reader.ExtractReadyUnit();
			reader.AddWaitingUnit(this);
			
			if (reader.Read() == xtd::good) {
				reader.ExtractWaitingUnit();
				
				m_rightOperand = reader.ExtractReadyUnit();
				reader.AddReadyUnit(this);
				
				return xtd::success::report();
			} else {
				reader.AddReadyUnit(m_leftOperand);
			
				return xtd::error::report("Нет правого операнда для оператора \"" + m_definition + "\".");
			}
		} else {
			reader.AddSuspendedUnit(this);
				
			return xtd::no_good; // Сигнал о приостановке чтения, но не об ошибке.
		}
	}
	else {
		return xtd::error::report("Нет левого операнда для оператора \"" + m_definition + "\".");
		//если нет левого оператора, то он == 0
	}
}

void CBinaryOperator::Save (std::ostream & stream) const {
	m_leftOperand->SaveAsExpression(stream);
	stream << " " << m_definition << " ";
	m_rightOperand->SaveAsExpression(stream);
}

CUnitSptr CDifference::Result() const {
	CUnitSptr t = m_leftOperand->Result() - m_rightOperand->Result();
	return t;
}

CUnitSptr CSum::Result() const {
	CUnitSptr t = m_leftOperand->Result() + m_rightOperand->Result();
	return t;
}

void CPower::Init () {
	m_definition = "^";
	m_priority = CPriority::Power;
}

} // namespace MM