#include "StdAfx.h"
#include <cmath>

#include "xtl/xtl_log.h"

#include "Operands.h"
#include "../../../../Common/2DPoint.h"

namespace MM {

bool CNumber::Load (CReader & reader) {
	m_value = atof(reader.CurrentToken().c_str());
	reader.AddReadyUnit(this);
	
	return xtd::success::report();
}

void CNumber::Init () {
	m_definition = ".0123456789";
	m_priority = CPriority::Number;
}

bool CCharacterLiteral::Load (CReader & reader) {
	reader.AddReadyUnit(this);
	
	return xtd::success::report();
}

void CPi::Init () {
	m_definition = "pi";
	m_priority = CPriority::Number;
	m_value = M_PI;
}

void CExponent::Init () {
	m_definition = "e";
	m_priority = CPriority::Number;
	m_value = M_E;
}

} // namespace MM