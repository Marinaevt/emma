#include "stdafx.h"
#include "2DDraggingOperation.h"

IOIMPL(C2DDraggingOperation, C2D_DRAGGING_OPERATION)

C2DDraggingOperation::C2DDraggingOperation(CString strName):CGen2DOperation(strName) {
	SetIconName(IDS_ICON_OPERATION_DRAGGING);
	//первый инструмент уже добавляется в CGen2DOperation
	AddTool();	//Инструмент для протаскивания, который "тянет"
}


C2DDraggingOperation::~C2DDraggingOperation(void)
{
}
