#pragma once
#include "emmadoc.h"
#include "3DMesh.h"
#include "I2DSketch.h"
#include "2DMesher.h"
#include "Resource.h"
#include "Child3DView.h"
#include "../../Common/Logger/Logger.h"
#include "RegisteredMessages.h"

class C3DMesher : public CEMMADoc
{
	GETTYPE (C3DMESHER)

public:
	C3DMesher(void);
	virtual ~C3DMesher(void);

	virtual void DrawGL(GLParam &parameter);
	void Construct3DMesh();
	bool GetTetra();
	Math::C2DPoint* Get2DEPoint(int e, int n);
////////////////////////////////////////////////
//Ромино для работы с панелькой
	void SetActiveObject(Math::C2DPoint *ppoint);
	void FillPropList(CMFCPropertyGridCtrl *pGrid);
	void UpdateProp(double *pEval, UNLONG *pId);

	void GetBoundingBox();
	CMFCPropertyGridProperty * AddProp(CStringW name, double eval, size_t pos);
	virtual bool IsPaneVisible(size_t nID) {return nID == IDR_3DMESHER_PANE;}

////////////////////////////////////////////////
public:
	IO::CBool m_flag;
protected:
	C3DMesh m_3DMesh;
};

IODEF (C3DMesher)