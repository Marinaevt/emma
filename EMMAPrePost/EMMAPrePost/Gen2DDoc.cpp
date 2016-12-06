#include "StdAfx.h"
#include "Gen2DDoc.h"


IOIMPL (CGen2DDoc, CGEN2DDOC)

CGen2DDoc::CGen2DDoc(void){
	m_pSpecimen = new C2DSpecimen;		//Класс для 2D задачи (Заготовка + КЭ сетка + Материал)
	m_pSpecimen->SetName(IDS_GEN2D_SPECIMEN);
	InsertSubDoc(m_pSpecimen);
	SetIconName(IDS_ICON_GEN2D);
}

CGen2DDoc::CGen2DDoc(int Icon) {
	m_pSpecimen = new C2DSpecimen;		//Класс для 2D задачи (Заготовка + КЭ сетка + Материал)
	m_pSpecimen->SetName(IDS_GEN2D_SPECIMEN);
	InsertSubDoc(m_pSpecimen);
	SetIconName(Icon);
}

CGen2DDoc::~CGen2DDoc(void){
}

//! Добавление новой операции через "+" (и диалоговое окно)
CEMMADoc *CGen2DDoc::Add(){
	CStringArray Names; //Types; //не используется
	size_t nSubDocNum = GetSubDocNum();
	for(size_t i=1; i<nSubDocNum; i++){
		Names.Add(SubDoc(i)->GetName());
	}

	CString strName;
	int nNum = strName.LoadString(IDS_OPERATION_NEW);

	DASD_INITSTRUCT IS = {&Names, nullptr, &strName, nSubDocNum };

	CDlgAddSubDoc Dlg(IS);

	if(Dlg.DoModal()==IDOK){
		CGen2DOperation* pGenOperation = new CGen2DOperation(Dlg.m_strName);
		return SubDoc(InsertSubDoc(pGenOperation));
	}
	return nullptr;
}

//! Переопределение загрузки, выделяем память.
bool CGen2DDoc::LoadBody(CStorage &file){
	
	//C2DSpecimen *pSpecimen = new C2DSpecimen;		//Класс для 2D задачи
	//pSpecimen->SetName(IDS_GEN2D_SPECIMEN);
	bool bLoad = CEMMADoc::LoadBody(file);
	m_pSpecimen = FindSpecimen();
	
	if (bLoad && m_pSpecimen) {
		return true;
	}
	
	return false;
}

C2DSpecimen *CGen2DDoc::FindSpecimen(){

	//! Ищем указатель на Заготовку (обычно индекс 0)
	for(size_t i=0; i<GetSubDocNum(); i++){
		C2DSpecimen *pSpecimen = dynamic_cast<C2DSpecimen *>(GetRoot()->SubDoc(i));
		if (pSpecimen) {
			return pSpecimen;	//! Нашли Заготовку, возвращаем.
		}
	}
	return nullptr;	//! Заготовки нет
}

C2DSketch *CGen2DDoc::GetSketch(){

	C2DSpecimen *pSpec = GetSpecimen();
	if (!pSpec) {
		CDlgShowError cError(ID_ERROR_GEN2DDOC_NO_SPECIMEN); //_T("Указатель на заготовку не найден"));
		return nullptr;
	}

	//! Проходим по Заготовке, ищем раздел КЭ сетки (обычно его индекс 1)
	for (size_t i=0; i<pSpec->GetSubDocNum(); i++){
		C2DSketch *pSketch = dynamic_cast<C2DSketch *>(pSpec->SubDoc(i));
		if (pSketch){
			return pSketch;	//! Нашли чертеж, возвращаем.
		}
	}
	return nullptr;	//! Чертежа нет
}

C2DMesher *CGen2DDoc::GetMesher(){
	
	C2DSpecimen *pSpec = GetSpecimen();
	if (!pSpec) {
		CDlgShowError cError(ID_ERROR_GEN2DDOC_NO_SPECIMEN); //_T("Указатель на заготовку не найден"));
		return nullptr;
	}

	//! Проходим по Заготовке, ищем раздел КЭ сетки (обычно его индекс 1)
	for (size_t i=0; i<pSpec->GetSubDocNum(); i++){
		C2DMesher *pMesher = dynamic_cast<C2DMesher *>(pSpec->SubDoc(i));
		if (pMesher){
			return pMesher;	//! Нашли создателя сетки, возвращаем.
		}
	}
	return nullptr;	//! Создателя сетки нет
}

const CMesh *CGen2DDoc::GetMesh(){

	if (GetMesher()){
		return GetMesher()->GetMesh();	//! Нашли сетку, возвращаем.
	}
	return nullptr;	//! КЭ сетки нет
}

CMaterialsList *CGen2DDoc::GetMaterialList() {
	C2DSpecimen *pSpec = GetSpecimen();

	size_t nDocNum = pSpec->GetSubDocNum();
	//! Проходим по Заготовке, ищем раздел Материалы (обычно индекс 2, т.е. третий)
	for (size_t i = 0; i < nDocNum; i++) {
		CMaterialsList *pMatList = dynamic_cast<CMaterialsList *>(pSpec->SubDoc(i));
		if (pMatList) {
			return pMatList;	//! Нашли список материалов
		}
	}
	return nullptr;	//! Списка материалов нет

	//return pSpec->GetMaterialList();
}

CGen2DOperation *CGen2DDoc::GetOperation(ptrdiff_t n){

	ptrdiff_t nSubDocNum = GetRoot()->GetSubDocNum();
	n = n+1;	// Первый поддокумент - Заготовка (не Операция)
	if (n > nSubDocNum) return nullptr; 

	CGen2DOperation *pOperation = dynamic_cast<CGen2DOperation*>(GetRoot()->SubDoc(n));
	if (pOperation){
		return pOperation;	//! Нашли Операцию по номеру, возвращаем.
	}
	return nullptr;	//! Операции нет
}

ptrdiff_t CGen2DDoc::AddOperationOnConstructor() {
	CString strName;
	int nNum = strName.LoadString(IDS_OPERATION_NEW);	//Сразу добавляем операцию
	CGen2DOperation* pGenOperation = new CGen2DOperation(strName);
	
	return InsertSubDoc(pGenOperation);
}