#pragma once
#include "emmadoc.h"
#include "Resource.h"
#include "MaterialDoc.h"

class CMaterialsList : public CEMMADoc
{
	GETTYPE(CMATERIALSLIST)
public:
	CMaterialsList(void);
	virtual ~CMaterialsList(void);

public:
	//! Добавляет материал
	virtual CEMMADoc *Add();
	virtual bool CanAdd()		/*inl*/ { return true; }
	//! Можем ли удалять поддокументы (т.е. материалы) 
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc) {
		if (dynamic_cast<CMaterialDoc *>(pDoc)) return true;
		return false;
	}

	//! Получаем материал из документа с номером i
	virtual CMaterial GetMaterial(size_t i);
	// Получаем документ материала с номером i
	CMaterialDoc *GetMaterialDoc(size_t i);
	// 
	CMaterialDoc *CMaterialsList::GetFirstMaterialDoc();
	
};
IODEF (CMaterialsList)