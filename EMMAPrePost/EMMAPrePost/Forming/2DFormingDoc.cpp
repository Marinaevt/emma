#include "stdafx.h"
#include "2DFormingDoc.h"
#include "2DFormingOperation.h"

IOIMPL (C2DFormingDoc, C2D_FORMING_DOC)

C2DFormingDoc::C2DFormingDoc(void){
	
	CString strName;
	int nNum = strName.LoadString(IDS_OPERATION_FORMING);	//Название для операции (Формование)
	InsertSubDoc(new C2DFormingOperation(strName));
	SetIconName(IDS_ICON_FORMING);
}


C2DFormingDoc::~C2DFormingDoc(void){
}

//! Добавление новой операции через "+" (и диалоговое окно)
CEMMADoc *C2DFormingDoc::Add() {
	CStringArray Names; //Types; //не используется
	size_t nSubDocNum = GetSubDocNum();
	for (size_t i = 1; i<nSubDocNum; i++) {
		Names.Add(SubDoc(i)->GetName());
	}

	CString strName;
	int nNum = strName.LoadString(IDS_OPERATION_FORMING);

	DASD_INITSTRUCT IS = { &Names, nullptr, &strName, nSubDocNum };

	CDlgAddSubDoc Dlg(IS);

	if (Dlg.DoModal() == IDOK) {
		C2DFormingOperation* pGenOperation = new C2DFormingOperation(Dlg.m_strName);
		return SubDoc(InsertSubDoc(pGenOperation));
	}
	return nullptr;
}