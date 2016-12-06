#include "StdAfx.h"
#include "ViewDescriptor.h"



int CViewDescriptor::s_nCounter = 0;
std::vector<CVDResource> s_arRes;

CViewDescriptor::CViewDescriptor(void)
{
	m_pView = nullptr;
	m_pViewClass = nullptr;
	s_nCounter++;
}

CViewDescriptor::CViewDescriptor(CRuntimeClass *pViewClass, UNINT nID)
{
	m_pView = nullptr;
	m_pViewClass = pViewClass;
	SetResource(nID);
	s_nCounter++;
}

CViewDescriptor::~CViewDescriptor(void){
	s_nCounter--;
	if(s_nCounter == 0){
		for(size_t i=0; i < s_arRes.size(); i++){
			FreeResource(m_Res.hMenu);
			FreeResource(m_Res.hAccel);
		}
		//s_arRes.SetSize(0);
		s_arRes.clear();
	}
}

void CViewDescriptor::SetResource(UNINT nID){
	m_Res.nIDR = nID;

	size_t i = s_arRes.size();
	for (size_t k = 0; k < i; k++) {	//раньше (до 07.03.2016) был while(i--){}
		if(nID == s_arRes[k].nIDR){
			m_Res = s_arRes[k];
			//break;
			return;	//нашли, записали, вышли
		}
	}

	//А если ничего не нашли
	//if(i<0){
		HINSTANCE hInst = AfxGetResourceHandle();
		m_Res.hMenu  = ::LoadMenu(hInst, MAKEINTRESOURCE(nID));
		m_Res.hAccel = ::LoadAccelerators(hInst, MAKEINTRESOURCE(nID));
		s_arRes.push_back(m_Res);
	//}
	
}