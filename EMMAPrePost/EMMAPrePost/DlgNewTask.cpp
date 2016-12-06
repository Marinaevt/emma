// DlgNewTask.cpp : implementation file
//

#include "stdafx.h"
#include "DlgNewTask.h"


// CDlgNewTask dialog

IMPLEMENT_DYNAMIC(CDlgNewTask, CDialog)

CDlgNewTask::CDlgNewTask(CWnd* pParent /*=nullptr*/)
	: CDialog(CDlgNewTask::IDD, pParent),
	m_nSelect(0)
{

}

CDlgNewTask::~CDlgNewTask()
{
}

/*
void CDlgNewTask::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_LBIndex(pDX, IDC_LIST1, m_nSelect);
}
*/


BEGIN_MESSAGE_MAP(CDlgNewTask, CDialog)
	
	//ListBOX (not ListCtrl)
	ON_LBN_SELCHANGE(IDC_LIST1, &CDlgNewTask::OnItemchangedList1)
END_MESSAGE_MAP()


// CDlgNewTask message handlers

int CDlgNewTask::OnInitDialog() {
	CDialog::OnInitDialog();

	FillList();

	return 1;	// return 1 unless you set the focus to a control
				// EXCEPTION: OCX Property Pages should return 0
}

//При смене выбора задачи меняем описание к ней
void CDlgNewTask::OnItemchangedList1(){

	CListBox *pList = (CListBox *)GetDlgItem(IDC_LIST1);
	m_nSelect = pList->GetCurSel();	//получаем номер выбранной задачи
	if (m_nSelect < 0 || m_nSelect >= m_nCountTasks) {
		CDlgShowError cError(ID_ERROR_DLGNEWTASK_NSELECT_WRONG); //_T("m_nSelect is wrong"));
		return;
	}

	//Загружаем текст описаний по ID и устанавливаем в область текста диалога
	CString strTmp;
	int nRes = strTmp.LoadStringW(m_nTasksDescript[m_nSelect]);
	GetDlgItem(IDC_STATIC)->SetWindowText(strTmp);
	
}

void CDlgNewTask::FillList() {
	CListBox *pList = (CListBox *)GetDlgItem(IDC_LIST1);
	
	if (pList == nullptr) {
		CDlgShowError cError(ID_ERROR_DLGNEWTASK_PLIST_NULL); //_T("pList is null"));
		return;
	}

	CString str;

	for (int i = 0; i < m_nCountTasks; i++) {
		int nNum = str.LoadString(m_nTasksID[i]);
		pList->AddString(str);
	}
}

size_t CDlgNewTask::GetSel() {

	if (m_nSelect == -1)
		m_nSelect = 0;
	
	if (m_nSelect < 0 || m_nSelect >= m_nCountTasks) {
		CDlgShowError cError(ID_ERROR_DLGNEWTASK_NSELECT_WRONG); //_T("m_nSelect is wrong"));
		return 0;
	}

	return m_nTasksID[m_nSelect];
}
