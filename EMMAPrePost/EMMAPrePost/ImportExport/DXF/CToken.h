#pragma once

#include <string>
#include "Constants.h"

namespace DXF {

//!	DXF-команда.
/*!
	Каждая команда состоит из двух частей: код группы, определяющий формат последующей команды и
	собственно команда.
	Возможные значения групповых кодов и команд описаны в заголовочном файле Constants.h.
*/
struct CToken {
	CToken () {
		m_code = COMMENT;
		m_command.clear();
	}
	~CToken () {}
	
	// -----------
	// Индикаторы.
	// -----------
	
	bool Is (GroupCode code, const char * command) const;

	bool IsSectionStart () const { return Is(ELEMENT_START, SECTION); }
		bool IsHeaderStart () const { return Is(ELEMENT_NAME, HEADER); }
		bool IsTablesStart () const { return Is(ELEMENT_NAME, TABLES); }
			bool IsTable () const { return Is(ELEMENT_START, TABLE); }
				bool IsLayerTable () const { return Is(ELEMENT_NAME, LAYER); }
				bool IsLTypeTable () const { return Is(ELEMENT_NAME, LTYPE); }
				bool IsStyleTable () const { return Is(ELEMENT_NAME, STYLE); }
				bool IsViewTable () const { return Is(ELEMENT_NAME, VIEW); }
				bool IsVPortTable () const { return Is(ELEMENT_NAME, VPORT); }
				bool IsUCSTable () const { return Is(ELEMENT_NAME, UCS); }
			bool IsEndOfTable () const { return Is(ELEMENT_START, END_OF_TABLE); }
		bool IsBlocksStart () const { return Is(ELEMENT_NAME, BLOCKS); }
		bool IsEntitiesStart () const { return Is(ELEMENT_NAME, ENTITIES); }
			bool IsLine () const { return Is(ELEMENT_START, LINE); }
			bool Is3DLine () const { return Is(ELEMENT_START, LINE3D); }
			bool IsPoint () const { return Is(ELEMENT_START, POINT); }
			bool IsCircle () const { return Is(ELEMENT_START, CIRCLE); }
			bool IsArc () const { return Is(ELEMENT_START, ARC); }
			bool IsTrace () const { return Is(ELEMENT_START, myTRACE); }
			bool IsSolid () const { return Is(ELEMENT_START, SOLID); }
			bool IsText () const { return Is(ELEMENT_START, TEXT); }
			bool IsShape () const { return Is(ELEMENT_START, SHAPE); }
			bool IsInsert () const { return Is(ELEMENT_START, INSERT); }
			bool IsAttrib () const { return Is(ELEMENT_START, ATTRIB); }
			bool IsPolyline () const { return Is(ELEMENT_START, POLYLINE); }
			bool IsVertex () const { return Is(ELEMENT_START, VERTEX); }
			bool IsSequenceEnd () const { return Is(ELEMENT_START, SEQEND); }
			bool Is3DFace () const { return Is(ELEMENT_START, FACE3D); }
			bool IsDimension () const { return Is(ELEMENT_START, DIMENSION); }
	bool IsSectionEnd () const { return Is(ELEMENT_START, END_OF_SECTION); }
	bool IsEndOfFile () const { return Is(ELEMENT_START, END_OF_FILE); }
	
	// -----------
	// Переменные.
	// -----------
	
	GroupCode m_code;
	std::string m_command;
};

}