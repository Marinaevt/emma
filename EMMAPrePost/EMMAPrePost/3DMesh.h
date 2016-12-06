#pragma once
#include "EmmaDoc.h"
#include "../../Common/Mesher/Mesh.h"
#include "Resource.h"
#include "3DElement4.h"
#include "../../Common/3DPoint.h"
#include "../../Common/AllTypes.h"
#include "../../Common/2DOutline.h"
#include "RegisteredMessages.h"

//TODO:  V а-та-та
#include "GL3DView.h"


typedef IO::CSimpleArray<CPoint3D,CIO3DNODEPTRARRAY> C3DNodeArray;
typedef IO::CSimpleArray<C3DElement4,CIO3DELEMENTPTRARRAY> C3DElementArray;
typedef IO::CSimpleArray<double,CIODOUBLE> C3DHeightArray;

class C3DMesh : public IO::CClass
{
	GETTYPE (C3DMESHER)

public:
	C3DNodeArray m_MeshPoints;
	C3DElementArray m_MeshElements;

protected:
	C3DHeightArray m_Z;
	IO::CInt m_LayersAmount;
	CRect3D m_BoundingBox;	//TODO: Надо преобразовать в IO

public:
	C3DMesh(void);
	int GetLayersAmount();
	void SetLayersAmount(int v);
	void ClearAll();
	void PrepareElements(int nEl);
	void SetCardinalityHeightArray();
	void SetCardinalityHeightArray(int m);

	void SetHeights();
	void SetHeights(double m);
	void SetHeight(int n, double h);
	double GetHeight(int n);
	void SetElements(int i, int v1, int v2, int v3, int v4);
	C3DElement4& GetElement(int n);
	CRect3D GetDraftRect(); //получить область для отрисовки
	void SetDraftRect(DBL mx, DBL Mx, DBL my, DBL My); //задать область для отрисовки
};
IODEF (C3DMesh)