#include "StdAfx.h"
#include "Gen2DOperation.h"
#include "2DBuilder.h"	//исключение

IOIMPL (CGen2DOperation, CGEN2DOPERATION)


//! Создаёт меню операций для 2D
CGen2DOperation::CGen2DOperation(CString strName):CEMMADoc(strName){
	
	SetIconName(IDS_ICON_OPERATION_2D);
	//Первый инструмент
	AddTool();

	//Сборка решения
	C2DBuilder *pBuilder = new C2DBuilder;
	pBuilder->SetName(IDS_BUILDER);
	InsertSubDoc(pBuilder);

	//Результаты
	C2DResults *pResults = new C2DResults;
	pResults->SetName(IDS_RESULTS);
	InsertSubDoc(pResults);

	m_nInstruments = 0;
	RegisterMember(&m_nInstruments);
}


CGen2DOperation::~CGen2DOperation(void){
}

//! Создаёт новый инструмент 2D
CEMMADoc *CGen2DOperation::Add(){
	return SubDoc(AddTool());
}

ptrdiff_t CGen2DOperation::GetToolNum(){
	
	ptrdiff_t count = 0;
	for (size_t i=0; i<GetSubDocNum(); i++){
		C2DTool *pTool = dynamic_cast<C2DTool *>(SubDoc(i));
		if (pTool){
			count++;
		}
	}
	return count;
}

C2DTool* CGen2DOperation::GetTool(ptrdiff_t i){

	ptrdiff_t count = 0;
	for (size_t j=0; j<GetSubDocNum(); j++){
		C2DTool *pTool = dynamic_cast<C2DTool *>(SubDoc(j));
		if (pTool){
			if (count==i){
				return pTool;
			}
			count++;
		}
	}

	CDlgShowError cError(ID_ERROR_GEN2DOPERATION_NO_TOOL); //_T("Инструмент не найден (GetTool)"));
	return 0;
}

ptrdiff_t CGen2DOperation::AddTool(CString strName){
	C2DTool *pTool = new C2DTool(strName);

	ptrdiff_t nInstruments = GetToolNum();

	//SetIconName(IDS_ICON_TOOL);

	if(strName.IsEmpty()){
		pTool->SetName(IDS_TOOL_DEF_NAME);
		pTool->SetName(pTool->GetName()+AllToString(" ")+AllToString(nInstruments+1));
	}
	return InsertSubDoc(pTool, nInstruments);
}