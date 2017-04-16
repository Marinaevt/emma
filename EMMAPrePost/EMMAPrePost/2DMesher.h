#pragma once
#include "EmmaDoc.h"
#include "I2DSketch.h"
#include "2DSketch.h"
#include "../../Common/2DOutline.h"
#include "../../Common/Mesher/Mesh.h"
#include "Resource.h"
#include "Child2DView.h"
#include "RegisteredMessages.h"


/*
	Класс документа для КЭ сетки в 2D
*/

class C2DMesher : public CEMMADoc, public I2DSketch
{
	GETTYPE (C2DMESHER)
public:
	C2DMesher(void);
	C2DMesher(CRuntimeClass *pViewClass, size_t nID);
	virtual ~C2DMesher(void);

	void SetOutline(C2DOutline *pOutline);
	void SetContours();
	C2DOutline *GetOutline(){return &m_Outline;}
///////////////////// Roma
	void SetActiveObject(Math::C2DPoint *ppoint){I2DSketch::SetActiveObject(ppoint);};
	void FillPropList(CMFCPropertyGridCtrl *pGrid);
	//void UpdateProp(double *pEval, UNLONG *pId);
	void UpdatePropList(CMFCPropertyGridCtrl *pGrid);	//обновляет данные документа (Doc) по новым значениям таблицы

	CMFCPropertyGridProperty * AddProp(CStringW name, double eval, size_t pos);

	void SetClearOutline();
	DBL GetInitLen();

	virtual void DrawGL(GLParam &parameter);
	virtual bool IsPaneVisible(size_t nID) {return nID == IDR_2DMESH_PANE;}
	virtual void *GetInterface(size_t nInterfaceID){
		return (nInterfaceID == IID_I2DSKETCH) ? dynamic_cast<I2DSketch *>(this) : CEMMADoc::GetInterface(nInterfaceID);
	}
	CMesh* GetMesh();

	void SetMeshType(UNINT id);
	UNINT GetMeshType();
	void AddMeshPoint(Math::C2DPoint ppoint);	//Установка точки сгущения
	void RegisterMembers();

public:
	bool GetBoundingBox(CRect2D &rect);
protected:
	C2DOutline m_Outline;
	CMesh m_Mesh;

};

IODEF (C2DMesher)