#pragma once
#include "emmadoc.h"
#include "../../Common/Material.h"
#include "MaterialView.h"
#include "Resource.h"

//-------------------------------------------------//

class CMaterialDoc : public CEMMADoc
{
	GETTYPE (CMATERIALDOC)
public:
	CMaterialDoc(void);
	virtual ~CMaterialDoc(void);

public:
	bool LoadMaterial(void);
	virtual void FillPropList(CMFCPropertyGridCtrl *pGrid);
	virtual void UpdateProp(double *pEval, UNLONG *pId);

	CMaterial GetMaterial() { return m_Material.Material(); };
	//возвращаем имя файла
	CString GetMaterialFileName() { return m_Material.FileName()(); };
	//возвращаем содержимое файла
	CString GetMaterialFileContent() { return m_Material.FileContent()(); }
	//устанавливаем содержимое файла
	void SetMaterialFileContent(CString s) { m_Material.SetFileContent(s); }

protected:
	virtual bool IsPaneVisible(size_t nID) {return (nID == IDR_MATERIALS_PANE);}
	virtual bool CanRemove(HTREEITEM hItem) { return true; }
	virtual void DrawGL();
	void MatToPoints();
	void PointsToMat();
	
protected:
	IO::CMaterial m_Material;
	//IO::CIOArray<CGL3DSurface, CGL3DSURFACE> m_surfaces;	// поверхности по температурам

};
IODEF (CMaterialDoc)