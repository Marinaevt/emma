#pragma once

#include <fstream>

#include "Constants.h"
#include "CToken.h"

namespace DXF {

//! Класс ввода-вывода DXF-файлов.
/*!
	Взаимодействует непосредственно с данными на жёстком диске, позволяя читать и записывать DXF-файлы.
*/
class CReader : public std::ifstream {
public:
	CReader (const char * filename, ios_base::openmode mode = ios_base::in) { open(filename, mode); }
	
	bool SeekForEntities ();
	
	GroupCode GetGroupCode ();
	std::string GetString ();
	CToken GetToken ();
	
	// "Подсматривает" следующий непробельный (не "white space") символ, не извлекая его из потока.
	int glance ();
	
	void SkipUntil (const char * command);
	void SkipEntity ();
	
	const CToken & NextToken () const { return m_nextToken; }
	
private:
	CToken m_nextToken;
};

}