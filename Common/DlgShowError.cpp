#include "stdafx.h"
#include "DlgShowError.h"


// диалоговое окно ошибок Dialog_show_error

CDlgShowError::CDlgShowError()
{
	m_StaticText.Empty();
}

CDlgShowError::CDlgShowError(CStringW s)
{
	m_StaticText = s;
	ShowErrorMessage();
}

CDlgShowError::CDlgShowError(bool bIf, CStringW s) {
	if (bIf) {
		m_StaticText = s;
		ShowErrorMessage();
	}
}

CDlgShowError::CDlgShowError(UNINT s) {
	int nRes = m_StaticText.LoadStringW(s);
	if (nRes != 1) {
		m_StaticText = _T("Error While Loading From Resources!");
	}
	ShowErrorMessage();
}

CDlgShowError::~CDlgShowError()
{
}


void CDlgShowError::InsertTextAndShowError(CStringW string) {
	m_StaticText = string;
	ShowErrorMessage();
}

void CDlgShowError::ShowErrorMessage() {
	
	//! VERSION FOR WINDOWS

	//Создаём модальное сообщение
	int msgboxID = MessageBoxW(
		nullptr,				//windows
		m_StaticText,			//text
		nullptr,				//caption
		MB_ICONSTOP | MB_OK | MB_DEFBUTTON1	//styles
		);

	//Создаём запись в TRACE
	TRACE(m_StaticText);

	//Создаём запись в лог


	/* !!! возможно на будущее, если нужны обработчики 
	switch (msgboxID)
	{
	case IDOK:	//OK button
		// TODO: add code
		break;
	case ...
	//*/

}