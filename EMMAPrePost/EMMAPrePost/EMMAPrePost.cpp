
// EMMAPrePost.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
//#include "vld.h"	//Visual Leak Detector (https://vld.codeplex.com/)
#include "EMMAPrePost.h"

#include "MainFrm.h"
#include "ChildFrm.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CEMMAPrePostApp

BEGIN_MESSAGE_MAP(CEMMAPrePostApp, CWinAppEx)
	ON_COMMAND(ID_APP_ABOUT, &CEMMAPrePostApp::OnAppAbout)
	//ON_COMMAND(ID_FILE_NEW, &CEMMAPrePostApp::OnFileNew)
END_MESSAGE_MAP()


// CEMMAPrePostApp construction

CEMMAPrePostApp::CEMMAPrePostApp()
{
	EnableHtmlHelp();
	
	/////// ROMA   Отключает в таскбаре Win7 вывод разных окон. (на будущее)
	////// Работает, начиная с VS2010!!!
	EnableTaskbarInteraction(0);

	m_bHiColorIcons = 1;
	m_nAppLook = 0;
	m_hAccelTable = nullptr;
	// Add construction code here,
	// Place all significant initialization in InitInstance
}

// The one and only CEMMAPrePostApp object

CEMMAPrePostApp theApp;


CEMMAPrePostApp *GetApp(){return &theApp;}

void CEMMAPrePostApp::UpdateAllViews(){
	CMainFrame *pMainFrame = dynamic_cast<CMainFrame *>(m_pMainWnd);
	pMainFrame->UpdateAllViews();
}

// CEMMAPrePostApp initialization
int CEMMAPrePostApp::InitInstance(){
	
	/* for DEBUG and memory leaks
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	_CrtSetBreakAlloc(9554);
	_CrtSetBreakAlloc(9553);
	_CrtSetBreakAlloc(9552);
	//*/

	PrePostStartup();
	//LOGGER.Init(CString(_T("..\\..\\Logs\\emmaprepost_initinstance.txt")));

	// InitCommonControlsEx() is required on Windows XP if an application
	// manifest specifies use of ComCtl32.dll version 6 or later to enable
	// visual styles.  Otherwise, any window creation will fail.
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// Set this to include all the common control classes you want to use
	// in your application.
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinAppEx::InitInstance();

	// Initialize OLE libraries
	if (!AfxOleInit())
	{
		CDlgShowError cError(IDP_OLE_INIT_FAILED);
		return 0;
	}
	AfxEnableControlContainer();
	// Standard initialization
	// If you are not using these features and wish to reduce the size
	// of your final executable, you should remove from the following
	// the specific initialization routines you do not need
	// Change the registry key under which our settings are stored
	// TODO: You should modify this string to be something appropriate
	// such as the name of your company or organization
	SetRegistryKey(_T("Local AppWizard-Generated Applications"));
	// LoadStdProfileSettings(4);  // Load standard INI file options (including MRU)  
   
	//Just add the following unregister() after LoadStfProfileSetting()  
	//Unregister(); 

	InitContextMenuManager();

	InitKeyboardManager();

	InitTooltipManager();
	CMFCToolTipInfo ttParams;
	ttParams.m_bVislManagerTheme = 1;
	theApp.GetTooltipManager()->SetTooltipParams(AFX_TOOLTIP_TYPE_ALL,
		RUNTIME_CLASS(CMFCToolTipCtrl), &ttParams);

	// To create the main window, this code creates a new frame window
	// object and then sets it as the application's main window object
	CMDIFrameWnd *pFrame = new CMainFrame;

	m_pMainWnd = pFrame;
	// create main MDI frame window
	if (!pFrame->LoadFrame(IDR_MAINFRAME)) {
		return 0;
	}
	
	//! Загрузка таблицы быстрых клавиш
	m_hAccelTable = LoadAccelerators(AfxGetInstanceHandle(), MAKEINTRESOURCE(IDR_ACCELERATOR_EMMA));

	// The main window has been initialized, so show and update it
	pFrame->ShowWindow(m_nCmdShow);
	pFrame->UpdateWindow();

	return 1;
}


// CEMMAPrePostApp message handlers

int CEMMAPrePostApp::ExitInstance() 
{
	return CWinAppEx::ExitInstance();
}

// CAboutDlg dialog used for App About

class CDlgAbout : public CDialog
{
public:
	CDlgAbout();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CDlgAbout::CDlgAbout() : CDialog(CDlgAbout::IDD)
{
}

void CDlgAbout::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CDlgAbout, CDialog)
END_MESSAGE_MAP()

// App command to run the dialog
void CEMMAPrePostApp::OnAppAbout()
{
	CDlgAbout aboutDlg;
	aboutDlg.DoModal();
}

// CEMMAPrePostApp customization load/save methods

void CEMMAPrePostApp::PreLoadState()
{
	int bNameValid;
	CString strName;
	bNameValid = strName.LoadString(IDS_EDIT_MENU);

	if (!bNameValid) {
		CDlgShowError cError(ID_ERROR_EMMAPREPOST_STRING_NOT_LOADED); //_T("Имя не загружено IDS_EDIT_MENU"));
		return;
	}

	GetContextMenuManager()->AddMenu(strName, IDR_POPUP_EDIT);
}

void CEMMAPrePostApp::LoadCustomState()
{
}

void CEMMAPrePostApp::SaveCustomState()
{
}

// CEMMAPrePostApp message handlers



int CEMMAPrePostApp::PreTranslateMessage(MSG* pMsg)
{
	// TODO: добавьте специализированный код или вызов базового класса
	if (m_hAccelTable) {
		if (::TranslateAccelerator(*m_pMainWnd, m_hAccelTable, pMsg)) {
			return 1;
		}
	}
	return CWinAppEx::PreTranslateMessage(pMsg);
}
