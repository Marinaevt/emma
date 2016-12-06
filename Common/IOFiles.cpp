#include "stdafx.h"
#include "IOFiles.h"

CStorage::CStorage(FILE *pOpenStream) : CStdioFile(pOpenStream) {
	//If you pass a FILE* stream to the CStdioFile constructor, then it will not auto-delete the file in the destructor!
	//To prevent this:
	m_bCloseOnDelete = 1;
}

bool SaveWideChar(CString path, CString string){

	// MFC only can handle ANSI-textfiles with CStdioFile, because it uses the default behaviour of the CRT.
	// But you can open the file - stream by yourself and then you can use the new encoding functionality in the CRT to specify the correct encoding.
	// Open the file with the specified encoding
	FILE *fStream;

	//"r" Opens for reading.If the file does not exist or cannot be found, the fopen_s call fails.
	//"w" Opens an empty file for writing.If the file exists, its contents are destroyed.
	//"a" Opens for writing at the end of the file(appending) without removing the EOF marker before writing new data to the file.Creates the file if it doesn't exist.
	//"r+" Opens for both reading and writing.(The file must exist.)
	//"w+" Opens an empty file for both reading and writing.If the file exists, its contents are destroyed.
	//"a+" Opens for reading and appending.The appending operation includes the removal of the EOF marker before new data is written to the file and the EOF marker is restored after writing is complete.Creates the file if it doesn't exist.
	// ccs == encoding from MS table
	
	//"wt" Open for Write in Text mode

	int e = _wfopen_s(&fStream, path, _T("wt,ccs=UNICODE"));
	if (e != 0) {
		return false; // failed.
	}
	
	CStorage f(fStream);  // open the file from this stream
	f.WriteString(string);
	f.Close();
	return true;
}

///////////////////////////////////////////////////////////////////////////////////
//GLobal
///////////////////////////////////////////////////////////////////////////////////
bool FindString(CStorage &File, const char* str) {
	File.Seek(0, CStorage::begin);
	char symb;	//wchar_t symb;
	int i = 0;
	while (str[i] != '\0') {

		//if(!File.Read((void*)&symb, sizeof(wchar_t))){
		if (!File.Read(&symb, 1)) {
			return false;
		}
		if (symb == str[i]) {
			i++;
		}
		else {
			i = 0;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////
bool FindString(CStorage &File, const CString &str) {
	File.Seek(0, CStorage::begin);
	
	CString sLine;
	while (File.ReadString(sLine))
	{
		if (sLine.Find(str) != -1) {
			return true;	//первое вхождение
		};
	}

	return false;
}

//////////////////////////////////////////////////////////////////////
bool read_dtoch(CStorage &File, double *d, char ch) {
	// считывает из File число d, при обнаружении символа ch возвращает false
	char symb = ' ';
	bool flag = false;
	int i = 0;
	char string[1000];

	do {
		//if (!) { break;	}

		if (symb == ch) break;

		if (symb == ',')symb = '.';
		if ((symb >= '0' && symb <= '9') || symb == '.' || symb == 'E' || symb == '-' || symb == '+') {
			flag = true;
		}
		else if (flag) {
			break;
		}

		if (flag) {
			string[i] = symb;
			i++;
			if (i >= 999) break;
		}
	} while (File.Read(&symb, 1)); //symb != EOF);

	string[i] = '\0';
	if (!flag) {
		return false;
	}
	*d = atof(string);
	return true;
}

//////////////////////////////////////////////////////////////////////
bool read_dtoch(CStorage &File, double *d, CString ch) {
	// считывает из File число d, при обнаружении символа ch возвращает false
	CString symb = _T(" ");
	bool flag = false;
	int i = 0;
	CString string;
	CString sLine;

	do {
		//if (!File.Read(&symb, sizeof(symb))) {
		//	break;
		//}

		int n = File.ReadString(sLine);
		symb = sLine[i];
		
		if (symb == ch) break;

		if (symb == ",") symb = ".";
		if ((symb >= _T("0") && symb <= _T("9")) || symb == _T(".") || symb == _T("E") || symb == _T("-") || symb == _T("+")) {
			flag = true;
		}
		else if (flag) {
			break;
		}

		if (flag) {
			string += symb;
			i++;
			if (i >= 999) break;	//слишком большое число
		}
	} while (symb != EOF);

	string += _T("\0");
	if (!flag) {
		return false;
	}
	*d = _ttof(string);	//конвертируем в double, for Unicode builds _wtof()
	return true;
}


//////////////////////////////////////////////////////////////////////
double read_d(CStorage *pFile, UNSHORT *pLastSymb) {
	char symb = ' ';
	bool flag = false;
	int i = 0;
	char string[1000];
	double ret;

	do {
		//if (!)break;

		if (symb == '\n' && pLastSymb != nullptr) break;

		if (symb == ',')symb = '.';
		if ((symb >= '0' && symb <= '9') || symb == '.' || symb == 'e' || symb == 'E' || symb == '-' || symb == '+') flag = true;
		else if (flag) break;

		if (flag) {
			string[i] = symb;
			i++;
			if (i >= 999) break;
		}
	} while (pFile->Read(&symb, 1)); //symb != EOF);
	string[i] = '\0';

	if (string[0] == '-' && string[1] == '-')
		ret = -1;
	else
		ret = atof(string);

	if (pLastSymb != nullptr) {
		*pLastSymb = MAKEWORD(symb, static_cast<char>(flag));
	}

	return ret;
}

//////////////////////////////////////////////////////////////////////
bool ReadValue(CStorage &file, const char *name, double *pVal) {
	if (!FindString(file, name)) {
		return false;
	}
	if (!read_dtoch(file, pVal, '<')) {
		return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////
bool ReadValue(CStorage &file, const CString &str, double *pVal) {
	if (!FindString(file, str)) {
		return false;
	}
	if (!read_dtoch(file, pVal, _T("<") )) {
		return false;
	}
	return true;
}

///////////////////////////////////////////////////////////
bool read_array(CStorage &File, char *str, std::vector<Math::C2DPoint> &Array) {
	if (!FindString(File, str)) {
		return false;
	}

	while (true) {
		Math::C2DPoint RP;
		if (!read_dtoch(File, &(RP.x), '<'))break;
		if (!read_dtoch(File, &(RP.y), '<'))break;
		Array.push_back(RP);
	}

	if (!Array.size()) {
		return false;
	}
	return true;
}
