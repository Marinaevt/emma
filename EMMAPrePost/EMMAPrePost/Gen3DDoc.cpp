#include "StdAfx.h"
#include "Gen3DDoc.h"


IOIMPL (CGen3DDoc, CGEN3DDOC)

CGen3DDoc::CGen3DDoc(void):CEMMADoc(RUNTIME_CLASS(CChild3DView), IDR_EMMAPrePostTYPE){
	SetIconName(IDS_ICON_GEN3D);
}

CGen3DDoc::CGen3DDoc(int Icon) : CEMMADoc(RUNTIME_CLASS(CChild3DView), IDR_EMMAPrePostTYPE) {
	SetIconName(Icon);
}

CGen3DDoc::~CGen3DDoc(void)
{
}
