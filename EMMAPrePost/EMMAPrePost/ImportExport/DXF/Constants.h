/*!
	\file Constants.h
	\brief Стандартные элементы файлов формата DXF.
	
	Файл описывает числовые и строковые константы, встречающиеся в файлах формата DXF.
*/

#pragma once

#include <iostream>

//! Пространство имён для элементов, связанных с файлами типа DXF.
/*!
	Включает классы, перечисления, константы и прочие штуки, имеющие непосредственное отношение к
	DXF-файлам.
*/
namespace DXF {

//! Групповые коды.
/*!
	Числа, с которых начинаются все команды, и которые определяют тип этих команд.
*/
enum GroupCode {
	ELEMENT_START = 0,			//!< Начало элемента.
	TEXT_VALUE = 1,				//!< Текстовое значение.
	ELEMENT_NAME = 2,			//!< Имя элемента 
	OTHER_TEXT_VALUES = 3,		//!< 3-5 Прочие текстовые значения или имена.
	LINE_TYPE = 6,				//!< Тип линии.
	TYPE_NAME = 7,				//!< Имя шрифта.
	LAYER_NAME = 8,				//!< Имя слоя.
	VARIABLE_NAME = 9,			//!< Имя переменной (используется только в заголовке).
	FIRST_X_COORDINATE = 10,	//!< Первая x-координата.
	SECOND_X_COORDINATE = 11,	//!< Вторая x-координата.
	LAST_X_COORDINATE = 18,		//!< 12-18 — прочие x-координаты.
	FIRST_Y_COORDINATE = 20,	//!< Первая y-координата.
	SECOND_Y_COORDINATE = 21,	//!< Вторая y-координата.
	LAST_Y_COORDINATE = 28,		//!< 22-28 — прочие x-координаты.
	FIRST_Z_COORDINATE = 30,	//!< Первая z-координата.
	SECOND_Z_COORDINATE = 31,	//!< Вторая z-координата.
	LAST_Z_COORDINATE = 37,		//!< 32-37 — прочие x-координаты.
	OBJECT_LEVEL = 38,			//!< Уровень объекта.
	// 39 — уточнить
	RADIUS = 40,				//!< Радиус окружности.
	START_WIDTH = 40,			//!< Ширина линии в начале сегмента.
	END_WIDTH = 41,				//!< Ширина линии в конце сегмента.
	CURVATURE = 42,				//!< Кривизна сегмента полилинии.
								//   43...48 разные параметры размеров.
	REPEAT = 49,				//!< Для организации таблиц с повторяющимися значениями (флаги 70-78 показывают число повторений 49-й группы).
	TANGENT_DIRECTION = 50,		//!< Касательное направление для сегмента полилинии.
	ANGLE = 50,					//!< Угловые данные.
	START_ANGLE = 50,			//!< Начальный угол дуги окружности.
	END_ANGLE = 51,				//!< Конечный угол дуги окружности.
	COLOR = 62,					//!< Номер цвета.
	ENTITIES_FLAG = 66,			//!< Флаг: 1 — если далее идут графические примитивы, 0 — если других элементов нет.
	FLAG70 = 70,
	FLAG71 = 71,
	FLAG72 = 72,
	FLAG73 = 73,
	COMMENT = 999				//!< Комментарий.
};

inline bool isXCoordinate (GroupCode code) {
	return (code >= FIRST_X_COORDINATE && code <= LAST_X_COORDINATE);
}

inline bool isYCoordinate (GroupCode code) {
	return (code >= FIRST_Y_COORDINATE && code <= LAST_Y_COORDINATE);
}

inline bool isZCoordinate (GroupCode code) {
	return (code >= FIRST_Z_COORDINATE && code <= LAST_Z_COORDINATE);
}

// ------------------------------------------------------------------------------------------------
// Команды.
// ------------------------------------------------------------------------------------------------

const char SECTION[] = "SECTION";
	const char HEADER[] = "HEADER";
	const char TABLES[] = "TABLES";
		const char TABLE[] = "TABLE";
			const char LAYER[] = "LAYER";
			const char LTYPE[] = "LTYPE";
			const char STYLE[] = "STYLE";
			const char VIEW[] = "VIEW";
			const char VPORT[] = "VPORT";
			const char UCS[] = "UCS";
		const char END_OF_TABLE[] = "ENDTAB";
	const char BLOCKS[] = "BLOCKS";
	const char ENTITIES[] = "ENTITIES";
		const char LINE[] = "LINE";
		const char LINE3D[] = "3DLINE";
		const char POINT[] = "POINT";
		const char CIRCLE[] = "CIRCLE";
		const char ARC[] = "ARC";
		const char myTRACE[] = "TRACE";
		const char SOLID[] = "SOLID";
		const char TEXT[] = "TEXT";
		const char SHAPE[] = "SHAPE";
		const char INSERT[] = "INSERT";
		const char ATTRIB[] = "ATTRIB";
		const char POLYLINE[] = "POLYLINE";
		const char VERTEX[] = "VERTEX";
		const char SEQEND[] = "SEQEND";
		const char FACE3D[] = "3DFACE";
		const char DIMENSION[] = "DIMENSION";
const char END_OF_SECTION[] = "ENDSEC";

const char END_OF_FILE[] = "EOF";

}