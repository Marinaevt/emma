#pragma once
#include "Resource.h"
#include "EmmaDoc.h"
#include "I2DSketch.h"
#include "../../Common/2DMotion.h"
#include "Child2DView.h"
#include "RegisteredMessages.h"

//! Класс Траекторий
class C2DTrajectory : public CEMMADoc, public I2DSketch
{
	GETTYPE (C2DTRAJECTORY)
public:
	C2DTrajectory(void);
	virtual ~C2DTrajectory(void);
	
	CIO2DMotion* GetMotion();

	void SetOutline(C2DOutline *pOutline);
	void SetContours();
	void SetClearOutline();

	void AddNode(Math::C2DPoint *ppoint){I2DSketch::AddNode(ppoint);};

	C2DOutline *GetOutline(){return dynamic_cast<C2DOutline *>(&m_Motion);}

	bool GetBoundingBox(CRect2D &rect);	
	/////
	void SetActiveObject(Math::C2DPoint *ppoint){I2DSketch::SetActiveObject(ppoint);};
	
	CMFCPropertyGridProperty* AddProp(CStringW name, double eval, size_t pos);
	void FillPropList(CMFCPropertyGridCtrl *pGrid);	// Заполнение параметров (правая панель)
	//void UpdateProp(double *pEval, UNLONG *pId);		// Обновление параметров (правая панель)
	void UpdatePropList(CMFCPropertyGridCtrl *pGrid);	//обновляет данные документа (Doc) по новым значениям таблицы
	
	virtual void DrawGL(GLParam &parameter);

	virtual void DrawTool(GLParam &parameter);

	virtual void InitialClearView();
	
	virtual bool IsPaneVisible(size_t nID) {return nID == IDR_2DTRAJECTORY_PANE;}

	virtual void *GetInterface(size_t nInterfaceID){
		return (nInterfaceID == IID_I2DSKETCH) ? dynamic_cast<I2DSketch *>(this) : CEMMADoc::GetInterface(nInterfaceID);
	}


protected:
	CIO2DMotion m_Motion;	//Двумерная траектория движения

public:


//	int m_nState;
//	bool f_play;
//	double m_dCurTime;
	
};

IODEF (C2DTrajectory)
