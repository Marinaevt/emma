#include "StdAfx.h"
#include "Parameter.h"

namespace MM {

bool CParameter::Load (CReader & reader) {
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

void CParameter::Save (std::ostream & stream) const {
	SaveBeginning(stream);
	m_value->SaveAsExpression(stream);
	SaveEnding(stream);
}

} // namespace MM