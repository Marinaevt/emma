#pragma once
#include "EmmaDoc.h"
#include "I2DSketch.h"
#include "../../Common/2DOutline.h"
#include "Resource.h"
#include "Child2DView.h"
#include "RegisteredMessages.h"


class C2DSketch : public CEMMADoc, public I2DSketch
{
public:
	GETTYPE (C2DSKETCH)
	C2DSketch(int Icon, CString strName = 0);
	C2DSketch(CRuntimeClass *pViewClass, size_t nID, int Icon, CString strName = 0);
	C2DSketch(CString strName = 0);
	C2DSketch(CRuntimeClass *pViewClass, size_t nID, CString strName = 0);
	virtual ~C2DSketch(void);

	void SetOutline(C2DOutline *pOutline);

	void SetClearOutline();
	void SetContours();
	void AddNode(Math::C2DPoint *ppoint){I2DSketch::AddNode(ppoint);};

	C2DOutline *GetOutline(){return &m_Outline;}

	virtual bool GetBoundingBox(CRect2D &rect);	

///////////////////// Roma
	void SetActiveObject(Math::C2DPoint *ppoint){I2DSketch::SetActiveObject(ppoint);};
	
	CMFCPropertyGridProperty* AddProp(CStringW name, double eval, size_t pos);
	void FillPropList(CMFCPropertyGridCtrl *pGrid);	// Заполняет таблицу свойств
	void UpdateProp(double *pEval, UNLONG *pId);		//обновляет данные документа (Doc)
	void UpdatePropList(CMFCPropertyGridCtrl *pGrid);	//обновляет все значения таблицы

	virtual void DrawGL(GLParam &parameter){I2DSketch::DrawGL(parameter);};

	//virtual void InitialClearView();

	virtual bool IsPaneVisible(size_t nID) {return nID == IDR_2DSKETCH_PANE;}
	virtual void *GetInterface(size_t nInterfaceID){
		return (nInterfaceID == IID_I2DSKETCH) ? dynamic_cast<I2DSketch *>(this) : CEMMADoc::GetInterface(nInterfaceID);
	}

protected:
	C2DOutline m_Outline;

};

IODEF (C2DSketch)
