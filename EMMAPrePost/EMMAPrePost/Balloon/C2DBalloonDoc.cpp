#include "stdafx.h"
#include "C2DBalloonDoc.h"
#include "C2DBalloonOperation.h"
IOIMPL(C2DBalloonDoc, C2D_BALLOON_DOC)
C2DBalloonDoc::C2DBalloonDoc()
{
	CString strName;
	int nNum = strName.LoadString(IDS_OPERATION_BALLOON);	//Название для операции (Формование)
	InsertSubDoc(new C2DBalloonOperation(strName));
	///////////////////////////////////////////////////////////////////////////
	// ATTENTION
	////////////////////////////////////////////////////////////////////
	SetIconName(IDS_ICON_FORMING);
}


C2DBalloonDoc::~C2DBalloonDoc()
{
}

//! Добавление новой операции через "+" (и диалоговое окно)
CEMMADoc *C2DBalloonDoc::Add() {
	CStringArray Names; //Types; //не используется
	size_t nSubDocNum = GetSubDocNum();
	for (size_t i = 1; i<nSubDocNum; i++) {
		Names.Add(SubDoc(i)->GetName());
	}

	CString strName;
	int nNum = strName.LoadString(IDS_OPERATION_BALLOON);

	DASD_INITSTRUCT IS = { &Names, nullptr, &strName, nSubDocNum };

	CDlgAddSubDoc Dlg(IS);

	if (Dlg.DoModal() == IDOK) {
		C2DBalloonOperation* pGenOperation = new C2DBalloonOperation(Dlg.m_strName);
		return SubDoc(InsertSubDoc(pGenOperation));
	}
	return nullptr;
}