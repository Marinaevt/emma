#pragma once
#include <afxtempl.h>
#include "Resource.h"
#include <vector>
#include "../../Common/AllTypes.h"
#include "../../Common/DlgShowError.h"

class CVDResource{
public:
	CVDResource(void):hMenu(nullptr), hAccel(nullptr), nIDR(0) {};
	~CVDResource(void){};

	HMENU	hMenu;			// Динамическое меню
	HACCEL	hAccel;			// Акселератор
	UNINT	nIDR;			// Идентификатор ресурсов связанный с дочерним фреймом документа, меню и акселератором
};

class CViewDescriptor
{
public:
	CViewDescriptor(void);
	CViewDescriptor(CRuntimeClass *pViewClass, UNINT nID);
	~CViewDescriptor(void);

public:
	CWnd *GetView(){return m_pView;}
	void SetView(CWnd *pView){m_pView = pView;}
	CRuntimeClass *GetViewClass(){return m_pViewClass;}
	HMENU	GetHMenu(){return m_Res.hMenu;}
	HACCEL	GetHAccel(){return m_Res.hAccel;}
	UNINT	GetResID(){return m_Res.nIDR;}

	void SetResource(UNINT nID);

protected:
	CWnd	*m_pView;           // Указатель на представление, ассоциированное с документом
	CRuntimeClass *m_pViewClass;// Класс представления
	CVDResource	m_Res;

	// static

	static int s_nCounter;
};
