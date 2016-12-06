#include "StdAfx.h"
#include "2DCalcStage.h"


IOIMPL (C2DCalcStage, C2DCALCSTAGE)

C2DCalcStage::C2DCalcStage(void) : CCalcStage(RUNTIME_CLASS(CCalcStage2DView), IDR_EMMAPrePostTYPE)
{
	m_2dstage_path() = _T("");	//Оставил, чтобы результаты не переделывать
	RegisterMember(&m_2dstage_path);	//TODO: закомментировать\убрать обе строки
}

C2DCalcStage::~C2DCalcStage(void){
}

// Достаточно вызвать LoadPathToTask() родителя (CCalcStage)
bool C2DCalcStage::LoadPathToTask(const CString &strFileName){
	if(!CCalcStage::LoadPathToTask(strFileName)) return false;

	m_2dstage_path() = strFileName;
	return true;
}

bool C2DCalcStage::GetBoundingBox(CRect2D &rect) {
	if (!GetTask()) return false;
	rect = GetTask()->GetBoundingBox();
	return true;
}