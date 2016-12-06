#pragma once
#include "AllTypes.h"
#include "2DPoint.h"

//! Класс для работы с файлами. Наследник CStdioFile, который наследник CFile
//! Пуст, т.к. будем дописывать свои методы в "ближайшем" будущем (сохранения\загрузки и прочее).
class CStorage : public CStdioFile
{
	//	DECLARE_DYNAMIC(CStorage);

	public:
		CStorage() {}
		CStorage(FILE *pOpenStream);
		~CStorage() {}
};

//Сохраняем содержимое с wide_char символами (UNICODE)
bool SaveWideChar(CString path, CString string);

bool FindString(CStorage &file, const char* str);
bool FindString(CStorage &File, const CString &str);

bool read_dtoch(CStorage &file, double *d, char ch);
bool read_dtoch(CStorage &File, double *d, CString ch);

double read_d(CStorage *pfile, UNSHORT *pLastSymb);
bool read_array(CStorage &File, char *str, std::vector<Math::C2DPoint> &Array);

bool ReadValue(CStorage &file, const char *name, double *pVal);
bool ReadValue(CStorage &file, const CString &str, double *pVal);
