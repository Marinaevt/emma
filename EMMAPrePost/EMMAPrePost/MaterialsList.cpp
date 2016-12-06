#include "StdAfx.h"
#include "MaterialsList.h"


IOIMPL (CMaterialsList, CMATERIALSLIST)

CMaterialsList::CMaterialsList(void){
	SetIconName(IDS_ICON_MATERIALS); 
}

CMaterialsList::~CMaterialsList(void){
}

CEMMADoc* CMaterialsList::Add(){
	CMaterialDoc *pNewMaterial = new CMaterialDoc;
	if(pNewMaterial->LoadMaterial()){
		return SubDoc(InsertSubDoc(pNewMaterial));
	}
	//если не удалось загрузить
	delete pNewMaterial;
	return nullptr;
}

CMaterial CMaterialsList::GetMaterial(size_t i) {
	
	if (i < GetSubDocNum()) {
		CMaterialDoc *pMatDoc = dynamic_cast<CMaterialDoc*>(SubDoc(i));
		return pMatDoc->GetMaterial();
	}

	return CMaterial();
}

CMaterialDoc *CMaterialsList::GetMaterialDoc(size_t i) {
	
	if (i < GetSubDocNum()) {
		CMaterialDoc *pMatDoc = dynamic_cast<CMaterialDoc*>(SubDoc(i));
		return pMatDoc;
	}
	return nullptr;
}

CMaterialDoc *CMaterialsList::GetFirstMaterialDoc() {

	size_t nNum = GetSubDocNum();
	for(size_t i = 0; i<nNum; i++){
		CMaterialDoc *pMatDoc = dynamic_cast<CMaterialDoc*>(SubDoc(i));
		if (pMatDoc != nullptr) return pMatDoc;
	}
	return nullptr;
}