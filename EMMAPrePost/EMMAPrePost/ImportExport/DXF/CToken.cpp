#include "StdAfx.h"
#include "CToken.h"

namespace DXF {

bool CToken::Is (GroupCode code, const char * command) const {
	if (m_code == code && m_command == command) return true;
	else return false;
}

} // namespace DXF