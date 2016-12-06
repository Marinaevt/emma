#include "StdAfx.h"
#include "SpecialExpressions.h"

namespace MM {

CUnitSptr CMultioperandExpression::Result () const {
	if (m_operands.size() == 1) {
		return (*m_operands.begin())->Result();
	}
	else {
		xtd::error::report("Неправильное использование мультиоперандного выражения.");
		// Для того, чтобы использовать несколько переменных, нужно явно преобразовать указатель на юнит
		// к указателю на мультиоперандное выражение (dynamic_cast<multioperand_expression &>(*указатель_на_выражение))
		// и обращаться непосредственно к массиву операндов.
	
		return nullptr;
	}
}

bool COpeningBracket::Load (CReader & reader) {
	CUnitSptr notToOwn;

	if (reader.HasReadyUnit()) {
		notToOwn = reader.GetReadyUnit();
	}
	
	reader.AddWaitingUnit(this);
	// Скобки притворяются оператором с низким приоритетом для того, чтобы искуственно завысить приоритеты
	// операторов внутри, и чтобы те не приостанавливали чтение, увидев ждущий оператор с более высоким приоритетом.
	// Эта хитрость возможна, т.к. рано или поздно (если, конечно, выражение составлено правильно)
	// должна встретиться закрывающая скобка, которая позволит правильно завершить чтение и восстановить приоритет.
	m_priority = CPriority::Lowest;
	
	if (reader.Read() == xtd::good) {		
		while (reader.HasReadyUnit() && reader.GetReadyUnit() != notToOwn) {
			m_operands.push_back(reader.ExtractReadyUnit());
		}
		
		reader.ExtractWaitingUnit();
		m_priority = CPriority::Brackets;
		reader.AddReadyUnit(this);
		
		return xtd::success::report();
	}
	
	//иначе ошибка
	return xtd::error::report("Нет операнда для оператора \"" + m_definition + "\".");
}

bool CClosingBracket::Load (CReader & reader) {
	if (reader.HasWaitingUnit()) {
		if (reader.GetWaitingUnit()->DefinitionIs("(") == false) {
			reader.AddSuspendedUnit(this);
		}
		
		return xtd::no_good;	// равносильно return true, без установки флага\статуса ошибки
	}
	else {
		return xtd::error::report("Оператор \"" + m_definition + "\" не может быть загружен.");
	}
}

void COpeningBracket::Save (std::ostream & stream) const {
	stream << "(";
	
	std::list<CUnitSptr>::const_iterator i = m_operands.begin();
	CUnitSptr pUnitSptr = *i;
	pUnitSptr->Save(stream);
	
	for (++i; i != m_operands.end(); ++i) {
		stream << ", ";
		pUnitSptr = *i;
		pUnitSptr->Save(stream);
	}

	stream << ")";
}

void COpeningBracket::Init () {
	m_definition = "(";
	m_priority = CPriority::Brackets;
}

void CClosingBracket::Init () {
	m_definition = ")";
	m_priority = CPriority::Brackets;
}

bool CComma::Load (CReader & reader) {
	if (reader.HasWaitingUnit()) {
		if (reader.GetWaitingUnit()->DefinitionIs("(")) {
			return xtd::success::report();
		}
		else {
			reader.AddSuspendedUnit(this);
		
			return xtd::no_good;	// равносильно return true, без установки флага\статуса ошибки
		}
	}
	else {
		return xtd::error::report("Разделитель «запятая» должен быть использован внутри скобок.");
	}
}

void CComma::Init() {
	m_definition = ",";
	m_priority = CPriority::Lowest;
}

bool CCommentStart::Load (CReader & reader) {
	while (reader.HasToken()) {
		if (FStringIsDefinitionOf(CCommentEnd())(reader.GetToken())) {
			break;
		}
	}

	return xtd::success::report();
}

} // namespace MM