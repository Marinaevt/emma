#pragma once
#include "calcstage.h"
#include "..\..\Common\2DTask.h"
#include "resource.h"

class C3DCalcStage : public CCalcStage{
	GETTYPE (C3DCALCSTAGE)
public:
	IO::CString m_3dstage_path;
	
public:
	C3DCalcStage(void);
	virtual ~C3DCalcStage(void);

	virtual bool LoadPathToTask(const CString &strFileName);

	virtual bool IsPaneVisible(size_t nID) /*inl*/ {return nID == IDR_3DCALCSTAGE_PANE;}
	virtual bool IsParamNeeded() { return true; }

	C2DTask *GetTask(){return dynamic_cast<C2DTask *>(m_pTask);}

	bool GetBoundingBox(CRect2D &rect);	//TODO: 3D rect?
};

IODEF (C3DCalcStage)