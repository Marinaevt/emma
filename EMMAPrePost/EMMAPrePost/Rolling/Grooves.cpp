#include "StdAfx.h"
#include "Grooves.h"

IOIMPL (CGrooves, CGROOVES)

CGrooves::CGrooves(void){
	//Верхний валок
	C2DSketch *pSkretch_top = new C2DSketch(IDS_ICON_TOPROLLER);
	pSkretch_top->SetName(IDS_TOPROLLER);
	InsertSubDoc(pSkretch_top);

	//Нижний валок
	C2DSketch *pSkretch_bot = new C2DSketch(IDS_ICON_BOTTOMROLLER);
	pSkretch_bot->SetName(IDS_BOTTOMROLLER);
	InsertSubDoc(pSkretch_bot);

	SetIconName(IDS_ICON_CALIBERS);
}

CGrooves::~CGrooves(void){
}
