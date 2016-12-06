#include "StdAfx.h"
#include "2DDraggingDoc.h"
#include "2DDraggingOperation.h"

IOIMPL(C2DDraggingDoc, C2D_DRAGGING_DOC)

C2DDraggingDoc::C2DDraggingDoc(void) {
	SetIconName(IDS_ICON_DRAGGING);
	CString strName;
	int nNum = strName.LoadString(IDS_OPERATION_DRAGGING);	//Название для операции (Волочение)
	InsertSubDoc(new C2DDraggingOperation(strName));
	
}

C2DDraggingDoc::~C2DDraggingDoc(void) {
}

//! Добавление новой операции через "+" (и диалоговое окно)
CEMMADoc *C2DDraggingDoc::Add() {
	CStringArray Names; //Types; //не используется
	size_t nSubDocNum = GetSubDocNum();
	for (size_t i = 1; i<nSubDocNum; i++) {
		Names.Add(SubDoc(i)->GetName());
	}

	CString strName;
	int nNum = strName.LoadString(IDS_OPERATION_DRAGGING);

	DASD_INITSTRUCT IS = { &Names, nullptr, &strName, nSubDocNum };

	CDlgAddSubDoc Dlg(IS);

	if (Dlg.DoModal() == IDOK) {
		C2DDraggingOperation* pGenOperation = new C2DDraggingOperation(Dlg.m_strName);
		return SubDoc(InsertSubDoc(pGenOperation));
	}
	return nullptr;
}