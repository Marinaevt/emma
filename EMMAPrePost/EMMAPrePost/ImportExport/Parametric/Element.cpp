#include "StdAfx.h"
#include "xtl/xtl_log.h"

#include "Element.h"
#include "SpecialExpressions.h"
#include "Compare.h"

namespace MM {

bool CElement::Load (CReader & reader, parametersContainer & parameters) {
	bool bLoadName = LoadName(reader);	//false == successful Load (o_O), true == unsuccessful (bad)
	if (bLoadName == xtd::no_good) {
		return xtd::error::report("Неверное имя");
	}
	
	if (PrepareForLoading(reader) == xtd::good) {
		CUnitSptr notToOwn = reader.HasReadyUnit() ? reader.GetReadyUnit() : nullptr;
	
		if (reader.Read() == xtd::good) {
			parametersContainer units;
		
			while (reader.HasReadyUnit() && (reader.GetReadyUnit() != notToOwn)) {
				units.push_back(reader.ExtractReadyUnit());
			}
			
			std::reverse(units.begin(), units.end());
			std::swap(units, parameters);
		}
		else {
			xtd::warning::report("Параметры элемента \"" + m_definition + "\" не указаны, поэтому будут установлены значения по-умолчанию.");
			return xtd::no_good;
		}
	}
	else {
		return xtd::error::report("Нет параметров для элемента \"" + m_definition + "\".");
	}
	
	return xtd::success::report();
}

bool CElement::LoadName (CReader & reader) {
	
	std::string sPeekToken = reader.PeekToken();
	//если НЕ открывающая скобка, то имя
	if (!FUnitDefinitionIs(sPeekToken)(COpeningBracket())) {
		m_name = reader.GetToken();
		return xtd::success::report();
	}
	
	return xtd::error::report("Couldn't load a name");
}

bool CElement::PrepareForLoading (CReader & reader) const {
	COpeningBracket * openingBracket = new COpeningBracket();
	openingBracket->SetPriority(CExpression::CPriority::Lowest);
	
	if (FUnitDefinitionIs(reader.PeekToken())(*openingBracket)) {
		reader.GetToken();
		reader.AddWaitingUnit(openingBracket);
		
		return xtd::success::report();
	}
	else {
		return xtd::error::report("Параметры каждого элемента должны быть указаны в скобках.");
	}
}

bool CElement::FinishLoading (CReader & reader) {
	reader.ExtractWaitingUnit(); // Скобка вытаскивается из стека ожидающих элементов.
	reader.AddReadyUnit(this);
	if (!NameIs("")) {
		reader.AddReadyElement(this);
	}
	else {
		return xtd::error::report("NameIs is Empty");
	}
	
	return xtd::success::report();
}

bool CElement::UndoLoading (CReader & reader, parametersContainer & parameters) {
//	reader.ExtractReadyUnit(); ?
//	reader.ExtractReadyElement(); ?

//	можно записать считанные ранее параметры в ридера.
	
	return xtd::error::report("Элементу \"" + m_definition + "\" передано неправильное число параметров.");
}

int CElement::ExportToOutline(C2DOutline * outline) const {
	xtd::warning::report("Элемент \"" + m_definition + "\"не поддерживает экспорт в C2DOutline.");

	return -1;
}

void CElement::SaveBeginning (std::ostream & stream) const {
	stream << m_definition;

	if (!NameIs("")) stream << " " << m_name;
	
	stream << " (";
}

void CElement::SaveDelimeter (std::ostream & stream) const {
	stream << ", ";
}

void CElement::SaveEnding (std::ostream & stream) const {
	stream << ")";
}
	
void CElement::SaveAsExpression (std::ostream & stream) const {
	if (!NameIs("")) stream << m_name;
	else Save(stream);
}

} // namespace MM