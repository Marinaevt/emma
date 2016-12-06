// EMMASrv.h : main header file for the EMMASrv application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "Resource.h"       // main symbols
#include "afxwinappex.h"
#include "TaskManager.h"
#include "MainFrm.h"

// для регистрации интерактивных объектов
#include "..\..\Common\2DPlaneFem.h"
#include "..\..\Common\2DRigid.h"

// для регистрации задач
#include "..\..\Common\2DTask.h"

// проверка FolderWork
#include "..\..\Common\FolderWork.h"

// проверка CLogger
#include "..\..\Common\Logger\Logger.h"

// проверка IOInterface
#include "..\..\Common\AllTypes.h"
#include "..\..\Common\IOInterface.h"
#include "..\..\Common\Strings.h"

// проверка CFunction
#include "..\..\Common\Function.h"

// проверка CIO2DOutline
#include "..\..\Common\2DOutline.h"
#include "..\..\Common\Startup.h"

#include "..\..\Common\DlgShowError.h"


// CEMMASrvApp:
// See EMMASrv.cpp for the implementation of this class
//

class CEMMASrvApp : public CWinApp
{
public:
	CEMMASrvApp();
	CTaskManager m_taskman;	// менеджер заданий

// Overrides
public:
	virtual int InitInstance();

// Implementation

public:
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CEMMASrvApp theApp;
