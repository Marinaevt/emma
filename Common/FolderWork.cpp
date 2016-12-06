#include "stdafx.h"
#include "FolderWork.h"

std::vector<CString> GetFolderListing(const CString& folder, const CString& filter, bool fullpath)
{
	CFileFind finder;
	int ret = finder.FindFile(folder + CString(_T("\\")) + filter);
	std::vector<CString> files;

	while (ret)
	{
		ret = finder.FindNextFile();

		if (finder.IsDots() || finder.IsDirectory())
			continue;

		if (fullpath)
			files.push_back(finder.GetFilePath());
		else
			files.push_back(finder.GetFileName());
	}

	finder.Close();
	return files;
}

CString GenFileName(UNINT filenumber, UNINT namelength)
{
	const wchar_t str1[2] = _T("0");
	const wchar_t str2[3] = _T("00");
	const wchar_t str3[4] = _T("000");
	const wchar_t str4[5] = _T("0000");
	const wchar_t str5[6] = _T("00000");
	const wchar_t str6[7] = _T("000000");
	const wchar_t str7[8] = _T("0000000");
	const wchar_t str8[9] = _T("00000000");
	const wchar_t str9[10] = _T("000000000");
	const wchar_t str10[11] = _T("0000000000");
	const wchar_t str11[12] = _T("00000000000");
	const wchar_t str12[13] = _T("000000000000");
	const wchar_t str13[14] = _T("0000000000000");
	const wchar_t str14[15] = _T("00000000000000");
	const wchar_t str15[16] = _T("000000000000000");
	const wchar_t str16[17] = _T("0000000000000000");

	const wchar_t *_str[16] = { str1, str2, str3, str4, str5, str6, str7, str8, str9, str10, str11, str12, str13, str14, str15, str16 };

	if (namelength == 0) namelength = 1;
	if (namelength > 16) namelength = 16;

	//if ((namelength == 0) || (namelength > 16))
	//	return CString(_T("error"));

	char buf[64];
	_itoa_s(filenumber, buf, 64, 10);

	size_t l = strlen(buf);

	CString ret(_T(""));

	if (l < namelength)
		ret.SetString(_str[namelength - l], namelength - l);

	return (ret + CString(buf));
}


bool DirectoryExists(CString szPath) {

	UNLONG dwAttrib = GetFileAttributes(szPath);

	return (dwAttrib != INVALID_FILE_ATTRIBUTES &&
		(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
}
