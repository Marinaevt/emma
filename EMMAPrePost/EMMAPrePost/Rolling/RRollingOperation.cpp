#include "StdAfx.h"
#include "RRollingOperation.h"

IOIMPL (CRRollingOperation, CR_ROLLING_OPERATION)

CRRollingOperation::CRRollingOperation(void){
	
	//Калибры (валки)
	CGrooves *pGrooves = new CGrooves;
	pGrooves->SetName(IDS_GROOVES);
	InsertSubDoc(pGrooves);

	//Сборка решения
	C2DBuilder *pBuilder = new C2DBuilder;
	pBuilder->SetName(IDS_BUILDER);
	InsertSubDoc(pBuilder);

	//Результаты
	C2DResults *pResults = new C2DResults;
	pResults->SetName(IDS_RESULTS);
	InsertSubDoc(pResults);

	SetIconName(IDS_ICON_ROLLINGOPERATION);
}

CRRollingOperation::~CRRollingOperation(void)
{
}
