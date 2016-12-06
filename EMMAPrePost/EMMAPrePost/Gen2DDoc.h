#pragma once
#include "emmadoc.h"
#include "Gen2DOperation.h"
#include "2DSpecimen.h"
#include "2DSketch.h"
#include "2DMesher.h"
#include "MaterialsList.h"
#include "Resource.h"
#include "DlgAddSubDoc.h"

//! 
class CGen2DDoc : public CEMMADoc
{
	GETTYPE (CGEN2DDOC)
public:
	CGen2DDoc(void);
	CGen2DDoc(int Icon);
	virtual ~CGen2DDoc(void);

public:
	//! Возвращает указатель на Заготовку, если она создана, иначе nullptr
	C2DSpecimen *GetSpecimen(){return m_pSpecimen;}
	C2DSpecimen *FindSpecimen();
	
	//! Возвращает указатель на Чертёж, если он создан, иначе nullptr
	C2DSketch *GetSketch();

	//! Возвращает указатель на создателя КЭ сетки, если он создан, иначе nullptr
	C2DMesher *GetMesher();

	//! Возвращает указатель на сетку, если она создана, иначе nullptr
	const CMesh *GetMesh();

	//! Возвращает указатель на список Материалов, если он создан, иначе nullptr
	CMaterialsList *GetMaterialList();	//TODO:

	//! Возвращает указатель на Операцию с номером n, если она создана, иначе nullptr
	CGen2DOperation *GetOperation(ptrdiff_t n);

	// Добавляет "Новую операцию" после работы конструктора из-за проблем с наследованием (формовка и волочение)
	ptrdiff_t AddOperationOnConstructor();

	//! Добавляет в задачу новую операцию с параметрами, задаваемыми пользователем
	virtual CEMMADoc *Add();

	// Возможно ли добавить подэлемент?
	virtual bool CanAdd() {return true;}
	// Возможно ли удалить подэлемент?
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc){
		if (dynamic_cast<CGen2DOperation *>(pDoc)) return true;
		return false;
	}
	virtual bool LoadBody(CStorage& file);

protected:
	C2DSpecimen *m_pSpecimen;
};

IODEF (CGen2DDoc)
