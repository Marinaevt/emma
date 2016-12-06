#include "StdAfx.h"
#include "CReader.h"

namespace DXF {

bool CReader::SeekForEntities () {
	GetToken();
	
	while (m_nextToken.IsSectionStart() == true) {
		GetToken();
		
		if (m_nextToken.IsEntitiesStart() == false) {
			SkipUntil(END_OF_SECTION);
			GetToken();
		} else {
			GetToken();
			return true;
		}
	}
	
	return false;
}

CToken CReader::GetToken () {
	CToken token = m_nextToken;

	m_nextToken.m_code = GetGroupCode();
	m_nextToken.m_command = GetString();

	return token;
}

GroupCode CReader::GetGroupCode () {
	GroupCode code;
	
	*this >> (int &)code;
	
	return code;
}

std::string CReader::GetString () {
	std::string string;
	char c;
	
	glance();
	
	while (!eof() && get(c)) {
		if (c == '\n' || c == '\r') break;
		else string.push_back(c);
	}
	
	return string;
}

int CReader::glance () {
	int c;
	
	while (isspace(c = peek())) get();
	
	return c;
}

void CReader::SkipUntil (const char * command) {
	std::string string = GetString();
	
	while (string != command) string = GetString();
}

void CReader::SkipEntity () {
	while (NextToken().m_code != ELEMENT_START) GetToken();
}


}