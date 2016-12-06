#include "StdAfx.h"
#include "RollingDoc.h"



IOIMPL (CRollingDoc, CROLLINGDOC)

CRollingDoc::CRollingDoc(void){
	C2DSpecimen *pSpecimen = new C2DSpecimen(IDS_ICON_PODKAT);
	pSpecimen->SetName(IDS_ROLLING_SPECIMEN);
	InsertSubDoc(pSpecimen);

	SetIconName(IDS_ICON_ROLLING);

	CString strName;
	int nNumR = strName.LoadString(IDS_OPERATION_ROLLING);
	InsertRollOperation(strName);
}

CRollingDoc::~CRollingDoc(void){
}

size_t CRollingDoc::InsertRollOperation(const CString &strName, ptrdiff_t nIndex/* = -1*/){
	CRRollingOperation *Operation = new CRRollingOperation;
	Operation->SetName(strName);
	return InsertSubDoc(Operation, nIndex);
}

CEMMADoc *CRollingDoc::Add(HTREEITEM hCurItem){
	size_t n = FindBranch(hCurItem);
	if (n == -2) {
		n = GetSubDocNum();
	}

	CStringArray Names, Types;
	for(size_t i=1; i<GetSubDocNum(); i++){
		Names.Add(SubDoc(i)->GetName());
	}
	Types.SetSize(2);
	int nNumR = Types[0].LoadString(IDS_ROLLING_R_OP);
	int nNumT = Types[1].LoadString(IDS_ROLLING_T_OP);

	CString strName;
	int nNumDEF = strName.LoadString(IDS_OPERATION_ROLLING);

	DASD_INITSTRUCT IS = {&Names, &Types, &strName, n};

	CDlgAddSubDoc Dlg(IS);

	if(Dlg.DoModal()==IDOK){
		ptrdiff_t nNew=0;
		if(Dlg.m_nType == 0){
			nNew = InsertRollOperation(Dlg.m_strName, Dlg.m_nIndex+1);
		}else{
			//nNew = InsertTRollOperation(Dlg.m_strName, Dlg.m_nIndex+1);
		}
		return SubDoc(nNew);
	}
	
	return nullptr;
}

CEMMADoc *CRollingDoc::Add() {

	CStringArray Names, Types;
	size_t nSubDocNum = GetSubDocNum();
	for (size_t i = 1; i<nSubDocNum; i++) {
		Names.Add(SubDoc(i)->GetName());
	}
	Types.SetSize(2);
	int nNumR = Types[0].LoadString(IDS_ROLLING_R_OP);
	int nNumT = Types[1].LoadString(IDS_ROLLING_T_OP);

	CString strName;
	int nNumDEF = strName.LoadString(IDS_OPERATION_ROLLING);

	DASD_INITSTRUCT IS = { &Names, &Types, &strName, nSubDocNum };

	CDlgAddSubDoc Dlg(IS);

	if (Dlg.DoModal() == IDOK) {
		size_t nNew = InsertRollOperation(Dlg.m_strName, Dlg.m_nIndex + 1);
		return SubDoc(nNew);
	}

	return nullptr;
}

bool CRollingDoc::CanRemove(HTREEITEM hCurItem){
	return FindBranch(hCurItem)>0;
}