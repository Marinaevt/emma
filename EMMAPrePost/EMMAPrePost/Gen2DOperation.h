#pragma once
#include "EmmaDoc.h"
#include "Resource.h"
#include "2DTool.h"
#include "2DResults.h"


class CGen2DOperation : public CEMMADoc
{
	GETTYPE (CGEN2DOPERATION)
public:
	CGen2DOperation(CString strName = 0);
	virtual ~CGen2DOperation(void);

public:
	virtual CEMMADoc *Add();
	virtual bool CanAdd()			/*inl*/ {return true;}
	//! Можем ли удалять поддокументы (т.е. инструменты) 
	virtual bool IsSubDocRemovable(CEMMADoc *pDoc){
		if (dynamic_cast<C2DTool *>(pDoc)) return true;
		return false;
	}
	
	ptrdiff_t GetToolNum();			//Кол-во Инструментов в Операции
	C2DTool* GetTool(ptrdiff_t i);	//Указатель на Инструмент 

protected:
	ptrdiff_t AddTool(CString strName = _T(""));
	IO::CInt m_nInstruments;			//ID созданных инструментов

};

IODEF (CGen2DOperation)
