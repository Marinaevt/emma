#include "StdAfx.h"
#include "2DSpecimen.h"

IOIMPL (C2DSpecimen, C2DSPECIMEN)

//! Заготовка
C2DSpecimen::C2DSpecimen(void){
	m_pSketch = new C2DSketch;	//Задание чертежа
	m_pSketch->SetName(IDS_SKETCH);
	size_t resSketch = InsertSubDoc(m_pSketch);	//Пункт меню
	//RegisterMember(m_pSketch);	//Возможность сохранить\загрузить

	m_pMesher = new C2DMesher;				//Задание КЭ сетки
	m_pMesher->SetName(IDS_MESHER);
	size_t resMesher = InsertSubDoc(m_pMesher);
	//RegisterMember(m_pMesher);

	m_pMaterialsList = new CMaterialsList;	//Задание материала
	m_pMaterialsList->SetName(IDS_MATERIALS);
	size_t resMatList = InsertSubDoc(m_pMaterialsList);
	//RegisterMember(m_pMaterials);
	SetIconName(IDS_ICON_SPECIMEN);
}

C2DSpecimen::C2DSpecimen(int icon) {
	m_pSketch = new C2DSketch;	//Задание чертежа
	m_pSketch->SetName(IDS_SKETCH);
	InsertSubDoc(m_pSketch);	//Пункт меню
								//RegisterMember(m_pSketch);	//Возможность сохранить\загрузить

	m_pMesher = new C2DMesher;				//Задание КЭ сетки
	m_pMesher->SetName(IDS_MESHER);
	InsertSubDoc(m_pMesher);
	//RegisterMember(m_pMesher);

	m_pMaterialsList = new CMaterialsList;	//Задание материала
	m_pMaterialsList->SetName(IDS_MATERIALS);
	InsertSubDoc(m_pMaterialsList);
	//RegisterMember(m_pMaterials);
	SetIconName(icon);
}

C2DSpecimen::~C2DSpecimen(void)
{
}

const CMesh* C2DSpecimen::GetMesh() {
	if (!m_pMesher || m_pMesher == nullptr) {
		CDlgShowError cError(ID_ERROR_2DSPECIMEN_MESHER_NULL); //_T("m_pMesher is null"));
		return nullptr;
	}
	return m_pMesher->GetMesh();
}

C2DOutline* C2DSpecimen::GetSketchOutline() {
	if (!m_pSketch || m_pSketch == nullptr) {
		CDlgShowError cError(ID_ERROR_2DSPECIMEN_SKETCH_NULL); //_T("m_pSketch is null"));
		return nullptr;
	}
	return m_pSketch->GetOutline();
}

CMaterial C2DSpecimen::GetMaterial(int i) {
	if (!m_pMaterialsList || m_pMaterialsList == nullptr) {
		CDlgShowError cError(ID_ERROR_2DSPECIMEN_MATERLIST_NULL); //_T("m_pMaterialsList is null"));
		return CMaterial();
	}
	return m_pMaterialsList->GetMaterial(i);
}