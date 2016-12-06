#include "stdafx.h"
#include "2DFormingOperation.h"

IOIMPL (C2DFormingOperation, C2D_FORMING_OPERATION)

C2DFormingOperation::C2DFormingOperation(CString strName):CGen2DOperation(strName){
	SetIconName(IDS_ICON_OPERATION_FORMING);
	//первый инструмент уже добавляется в CGen2DOperation
	AddTool();	//добавляем второй (крепёж снизу)
}


C2DFormingOperation::~C2DFormingOperation(void)
{
}
