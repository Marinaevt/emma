#include "StdAfx.h"
#include <list>

#include "SpecialElements.h"
#include "Compare.h"

namespace MM {

bool CVariable::Load (CReader & reader) {
	parametersContainer parameters;

	if (CElement::Load(reader, parameters) == xtd::good) {
		if (parameters.size() == 1) {
			m_value = parameters[0];
		}
		else if (parameters.size() != 0) {
			return UndoLoading(reader, parameters);
		}
	}
	else {
		return xtd::error::report("Элемент \"" + m_definition + "\" описан неправильно.");
	}
	
	return FinishLoading(reader);
}

void CVariable::Save (std::ostream & stream) const {
	SaveBeginning(stream);
	m_value->SaveAsExpression(stream);
	SaveEnding(stream);
}

void CVariable::Init () {
	m_definition = "variable";
	m_name.clear();
}

bool CVariableLoader::Load (CReader & reader) {
	typedef xtd::stack<CUnitSptr>::container container_type;
	const container_type & elements = reader.GetReadyElements();
	
	std::string sReaderToken;
	for (container_type::const_iterator i = elements.begin(); i != elements.end(); ++i) {
		
		CUnitSptr pUnitSptr = *i;
		sReaderToken = reader.CurrentToken();

		if (FElementNameIs(sReaderToken)(dynamic_cast<const CElement &>(*pUnitSptr))) {
			CVariable * variable = new CVariable(pUnitSptr);
			variable->SetName(sReaderToken);
			
			reader.AddReadyUnit(variable);
			
			return xtd::success::report();
		}
	}
	
	sReaderToken = reader.CurrentToken();
	return xtd::error::report("Элемент с именем \"" + sReaderToken + "\" не найден.");
}

bool CElementEnd::Load (CReader & reader) {
	if (reader.HasReadyUnit()) {
		reader.ExtractReadyUnit();
		return xtd::success::report();
	}

	return xtd::error::report("Перед символом завершения должен стоять хотя бы один элемент.");
}

void CVariableLoader::Init () {
	m_definition = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
}

} // namespace MM