#include "StdAfx.h"
#include "2DResults.h"

IOIMPL (C2DResults, C2DRESULTS)

C2DResults::C2DResults(void){
	Init2DChildViewDescriptor(RUNTIME_CLASS(CCalcStage2DView), IDR_EMMAPrePostTYPE);
	Init3DChildViewDescriptor(RUNTIME_CLASS(CCalcStage3DView), IDR_EMMAPrePostTYPE);
	//Загрузка списка результирующих файлов
	//LoadResults();
	SetIconName(IDS_ICON_RESULT);
}

C2DResults::~C2DResults(void){
	//m_C2DCalcStageArr.clear();
}

bool C2DResults::LoadStage(const CString &strFileName){
	C2DCalcStage *p2DStage = new C2DCalcStage;
	if (!p2DStage->LoadPathToTask(strFileName)) {
		delete p2DStage;
		return false;
	}
	p2DStage->LoadTask();
	p2DStage->SetName(MakeName(p2DStage->m_2dstage_path()));
	InsertSubDoc(p2DStage);
	return true;
}

