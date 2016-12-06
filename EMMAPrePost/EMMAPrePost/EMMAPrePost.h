
// EMMAPrePost.h : main header file for the EMMAPrePost application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "Resource.h"       // main symbols
#include "afxwinappex.h"
#include "PrePostStartup.h"
#include "../../Common/Logger/Logger.h"



// CEMMAPrePostApp:
// See EMMAPrePost.cpp for the implementation of this class
//

class CEMMAPrePostApp : public CWinAppEx
{
public:
	CEMMAPrePostApp();

public:
	void UpdateAllViews();

// Overrides
public:
	virtual int InitInstance();
	virtual int ExitInstance();

// Implementation
protected:
//	HMENU  m_hMDIMenu;
//	HACCEL m_hMDIAccel;

public:
	UNINT  m_nAppLook;	//Стиль окон, кнопок и прочего
	int  m_bHiColorIcons;
	HACCEL  m_hAccelTable;	//Таблица быстрых клавиш

	virtual void PreLoadState();
	virtual void LoadCustomState();
	virtual void SaveCustomState();

	afx_msg void OnAppAbout();
//	afx_msg void OnFileNew();
	DECLARE_MESSAGE_MAP()

	int PreTranslateMessage(MSG* pMsg);
};

CEMMAPrePostApp *GetApp();

extern CEMMAPrePostApp theApp;
