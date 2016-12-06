#pragma once
#include "EmmaDoc.h"
#include "2DMesher.h"
#include "2DSketch.h"
#include "MaterialsList.h"
#include "Resource.h"

class C2DSpecimen :	public CEMMADoc
{
	GETTYPE (C2DSPECIMEN)
public:
	C2DSpecimen(void);
	C2DSpecimen(int icon);
	virtual ~C2DSpecimen(void);

public:
	C2DOutline *GetSketchOutline();
	const CMesh *GetMesh();
	CMaterial GetMaterial(int i);	//не указатель

	//2D чертёж
	C2DSketch *GetSketch(){return m_pSketch;}
	//2D КЭ сетка
	C2DMesher *GetMesher(){return m_pMesher;}
	//Список Материалов
	CMaterialsList *GetMaterialList() { return m_pMaterialsList; }

	virtual bool CanRemove(HTREEITEM hItem){return false;}

protected:
	C2DMesher *m_pMesher;
	C2DSketch *m_pSketch;
	CMaterialsList *m_pMaterialsList;
};

IODEF (C2DSpecimen)