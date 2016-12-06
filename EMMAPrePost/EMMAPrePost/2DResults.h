#pragma once
#include "Results.h"
#include "2DCalcStage.h"
#include "CalcStage2DView.h"
#include "CalcStage3DView.h"
#include "resource.h"

class C2DResults :	public CResults
{
	GETTYPE (C2DRESULTS)
public:
	C2DResults(void);
	virtual ~C2DResults(void);
	
	virtual bool LoadStage(const CString &strFileName);

	//! Все стадии вычислений (Results->CalcStage) лежат в этом массиве
 // std::vector<C2DCalcStage *> m_C2DCalcStageArr;
};

IODEF (C2DResults)
