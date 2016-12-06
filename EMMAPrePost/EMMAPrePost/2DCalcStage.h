#pragma once
#include "CalcStage.h"
#include "../../Common/2DTask.h"
#include "Resource.h"
#include "RegisteredMessages.h"
#include <vector>
#include "CalcStage2DView.h"

class C2DCalcStage : public CCalcStage{
	GETTYPE (C2DCALCSTAGE)
public:
	IO::CString m_2dstage_path;
	
public:
	C2DCalcStage(void);
	virtual ~C2DCalcStage(void);

	virtual bool LoadPathToTask(const CString &strFileName);

	virtual bool IsPaneVisible(size_t nID) /*inl*/ {return nID == IDR_2DCALCSTAGE_PANE || nID == IDR_3DCALCSTAGE_PANE;}

	C2DTask *GetTask(){return dynamic_cast<C2DTask *>(m_pTask);}

	bool GetBoundingBox(CRect2D &rect);
};

IODEF (C2DCalcStage)