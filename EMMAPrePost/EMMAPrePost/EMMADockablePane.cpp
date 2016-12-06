#include "StdAfx.h"
#include "EMMADockablePane.h"

IMPLEMENT_DYNAMIC(CEMMADockablePane, CDockablePane)

CEMMADockablePane::CEMMADockablePane(size_t nID) : m_nResID(nID), m_pDoc(nullptr)
{
}

CEMMADockablePane::~CEMMADockablePane(void)
{
}

void CEMMADockablePane::SetDocument(CEMMADoc *pDoc){
	m_pDoc = pDoc;
}