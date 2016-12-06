// EMMASrv.cpp : Defines the class behaviors for the application.
//
#include "stdafx.h"
//#include "vld.h"	//Visual Leak Detector (https://vld.codeplex.com/)
#include "EMMASrv.h"
#include "../../Common/Strings.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CEMMASrvApp

BEGIN_MESSAGE_MAP(CEMMASrvApp, CWinApp)
	ON_COMMAND(ID_APP_ABOUT, &CEMMASrvApp::OnAppAbout)
END_MESSAGE_MAP()


// CEMMASrvApp construction

CEMMASrvApp::CEMMASrvApp()
{

	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}

// The one and only CEMMASrvApp object

CEMMASrvApp theApp;


// CEMMASrvApp initialization

BOOL CEMMASrvApp::InitInstance(){
	// InitCommonControlsEx() is required on Windows XP if an application
	// manifest specifies use of ComCtl32.dll version 6 or later to enable
	// visual styles.  Otherwise, any window creation will fail.
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// Set this to include all the common control classes you want to use
	// in your application.
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

	// Initialize OLE libraries
	if (!AfxOleInit())
	{
		CDlgShowError cError(ID_ERROR_EMMASRV_OLE_INIT_FAILED);	//load text from resource id
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
	SetRegistryKey(_T("Emma computation server"));

	// To create the main window, this code creates a new frame window
	// object and then sets it as the application's main window object
	CMainFrame* pFrame = new CMainFrame;

	m_pMainWnd = pFrame;
	// create and load the frame with its resources
	pFrame->LoadFrame(IDR_MAINFRAME, WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE, nullptr,	nullptr);
	
	// The one and only window has been initialized, so show and update it
	//pFrame->ShowWindow(SW_SHOW);
	//pFrame->UpdateWindow();
	// call DragAcceptFiles only if there's a suffix
	//  In an SDI app, this should occur after ProcessShellCommand


	Startup();

	//Aksenov считываем Задание из файла

	CString strTaskFName, strTaskFPath;// = _T("..\\..\\Dump\\test_stamping1\\file00000.emma_task");
	
	

#ifndef __LOCAL
	CFileDialog fileDlg(TRUE, _T("emma_task"),0, 4|2, _T("task file (*.emma_task)|*.emma_task|"));
	if(fileDlg.DoModal() == IDOK){
		strTaskFName = fileDlg.GetPathName();
		strTaskFPath = strTaskFName;
	}else{
		delete pFrame;
		return 0;
	}
#else
	strTaskFPath = "Calculations\\forming.emma_task";
	strTaskFName = "Calculations\\forming.emma_task";
#endif

	CStorage TaskFile;
	ITask* pTask;

	if(TaskFile.Open(strTaskFName, CStorage::modeRead | CStorage::typeBinary)){

		pTask = dynamic_cast<ITask *>(IO::LoadObject(TaskFile));
		
		if (pTask == nullptr) {
			
			CDlgShowError cError(ID_ERROR_EMMASRV_CANT_LOAD_TASKFILE); //_T("File can't be loaded (.emma_task)"));
			PostMessage(pFrame->m_nid.hWnd, WM_CLOSE, 0, 0);
			delete pFrame;
			return 0;
		}

		int nPosSlash = strTaskFName.ReverseFind( '\\' ); 
		strTaskFName = strTaskFName.Right(strTaskFName.GetLength() - nPosSlash - 1); //2 слеша
		strTaskFPath.Truncate( nPosSlash+1 ); 

		pTask->SetFileName(strTaskFName);		//устанавливаем имя файла
		pTask->SetFilePath(strTaskFPath);		//устанавливаем путь до директории
	
		TaskFile.Close();
	}else{
		//Файл не открылся. Нужен обработчик ошибок.
		CDlgShowError cError(ID_ERROR_EMMASRV_NO_EMMATASKFILE); //_T("emma_task not found (Calculations\\forming.emma_task)"));
		PostMessage(pFrame->m_nid.hWnd, WM_CLOSE, 0, 0);
		delete pFrame;
		return 0;
	}
	
	if(!pTask->Init()){
		CDlgShowError cError(ID_ERROR_EMMASRV_INIT_EMMATASKFILE);
		PostMessage(pFrame->m_nid.hWnd, WM_CLOSE, 0, 0);
		delete pFrame;
		return 0;
	};
	
	//TODO: Где-то тут должен поселиться TaskManager
	
	C2DTask* p2DTask = dynamic_cast<C2DTask* >(pTask);
	if (p2DTask){
		p2DTask->RunCalc();	// расчёт 2D задания
	}

	
	Beep( 750, 1500 );
	PostMessage(pFrame->m_nid.hWnd, WM_CLOSE, 0, 0);	//закрываем иконку
	delete pFrame;
	return 0;	//Задание расчитано, закрываем сервер
}


// CEMMASrvApp message handlers




// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()

// App command to run the dialog
void CEMMASrvApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}

// CEMMASrvApp message handlers



