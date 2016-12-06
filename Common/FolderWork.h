#pragma once
#include "afxtempl.h"
#include "afx.h"
//#include "cstringt.h"
#include "AllTypes.h"
#include "Tchar.h"
#include <vector>

/*! 
    \brief Получить список файлов папки.
    \param folder Папка, файлы которой просматриваются.
    \param filter Маска, по которой фильтруются имена файлов (по-умолчанию, "*.*").
	\param fullpath true = (в списке полные пути файлов), false = (только имена) (по-умолчанию, true).
*/
std::vector<CString> GetFolderListing(const CString& folder, const CString& filter, bool fullpath = true);
/*! 
    \brief Генерирует имя файла вида "000001" по заданному номеру, заданной длины.
    \param filenumber Номер файла.
    \param namelength Длина имени файла от 1 до 16 (по-умолчанию, 16).
*/
CString GenFileName(UNINT filenumber, UNINT namelength = 16);

/*!
	\brief проверяет, существует ли директория по данному пути szPath
	\param szPath Путь до папки
*/
bool DirectoryExists(CString szPath);
