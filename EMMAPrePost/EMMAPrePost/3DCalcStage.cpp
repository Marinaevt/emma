#include "StdAfx.h"
#include "3DCalcStage.h"
#include "CalcStage3DView.h"
#include "resource.h"
#include "RegisteredMessages.h"
#include <vector>

IOIMPL (C3DCalcStage, C3DCALCSTAGE)

C3DCalcStage::C3DCalcStage(void) : CCalcStage(RUNTIME_CLASS(CCalcStage3DView), IDR_EMMAPrePostTYPE)
{
	m_3dstage_path() = _T("");
	RegisterMember(&m_3dstage_path);
}

C3DCalcStage::~C3DCalcStage(void){
}

bool C3DCalcStage::LoadPathToTask(const CString &strFileName){
	if(!CCalcStage::LoadPathToTask(strFileName)) return false;

	m_3dstage_path() = strFileName;
	return true;
}

//TODO: 3D rect?
bool C3DCalcStage::GetBoundingBox(CRect2D &rect) {
	if (!GetTask()) return false;
	rect = GetTask()->GetBoundingBox();
	return true;
}