#pragma once
#include "AllTypes.h"	//для вытаскивания из ресурсов по UNINT ID 
#include "Strings.h"	//для преобразования AllToString

// окно Dialog_show_error
class CDlgShowError
{

public:
	// стандартный конструктор
	CDlgShowError();
	// сразу вызывает модальное окно с текстом ошибки 
	CDlgShowError(CStringW s);
	// вызывает модальное окно с текстом ошибки при выполнении условия bIf
	CDlgShowError(bool bIf, CStringW s);
	// вызывает модальное окно с текстом ошибки по ID ресурса (EMMAPrePost.rc, EMMASrv.rc)
	CDlgShowError(UNINT s);
	virtual ~CDlgShowError();
	
	// сохраняет текст ошибки и сразу вызывает ShowErrorMessage()
	void InsertTextAndShowError(CStringW string);
	// отображение текста ошибки в модальном окне и запись текста в вывод отладчика
	void ShowErrorMessage();

public:
	CStringW m_StaticText;

};