#include "stdafx.h"
#include "C2DBalloonOperation.h"


IOIMPL(C2DBalloonOperation, C2D_BALLOON_OPERATION)

C2DBalloonOperation::C2DBalloonOperation(CString strName) : CGen2DOperation(strName){
	SetIconName(IDS_ICON_OPERATION_FORMING);
	//первый инструмент уже добавляется в CGen2DOperation
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//AddTool();	//добавляем второй (крепёж снизу)
}


C2DBalloonOperation::~C2DBalloonOperation()
{
}
