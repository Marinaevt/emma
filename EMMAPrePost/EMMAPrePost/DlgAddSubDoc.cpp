// DlgAddSubDoc.cpp : implementation file
//

#include "stdafx.h"
#include "DlgAddSubDoc.h"


// CDlgAddSubDoc dialog

IMPLEMENT_DYNAMIC(CDlgAddSubDoc, CDialog)

CDlgAddSubDoc::CDlgAddSubDoc(const DASD_INITSTRUCT &IS, CWnd* pParent /*=nullptr*/)
	: CDialog(CDlgAddSubDoc::IDD, pParent)
	, m_nType(0)
{
	m_pList = IS.pList;
	m_pTypes = IS.pTypes;
	m_nIndex = IS.nIndex;
	m_strName = *IS.pDefaultName;
}

CDlgAddSubDoc::~CDlgAddSubDoc()
{
}

void CDlgAddSubDoc::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, m_strName);
	DDX_Control(pDX, IDC_LIST1, m_ListBox);
	DDX_Control(pDX, IDC_COMBO1, m_ComboBox);
	DDX_CBIndex(pDX, IDC_COMBO1, m_nType);
	DDX_Control(pDX, IDC_S_TYPE, m_staticType);
}


BEGIN_MESSAGE_MAP(CDlgAddSubDoc, CDialog)
	ON_LBN_SELCHANGE(IDC_LIST1, &CDlgAddSubDoc::OnLbnSelchangeList1)
	ON_EN_CHANGE(IDC_EDIT1, &CDlgAddSubDoc::OnEnChangeEdit1)
	ON_CBN_SELCHANGE(IDC_COMBO1, &CDlgAddSubDoc::OnCbnSelchangeCombo1)
END_MESSAGE_MAP()


// CDlgAddSubDoc message handlers
int CDlgAddSubDoc::OnInitDialog(){
	CDialog::OnInitDialog();

	for(ptrdiff_t i=0; i<m_pList->GetSize(); i++){
		m_ListBox.AddString(m_pList->GetAt(i));
	}
	if(m_nIndex<0 || m_nIndex > m_pList->GetSize()){
		m_nIndex = m_pList->GetSize();
	}
	//m_ListBox.InsertString(m_nIndex, m_strName);
	m_ListBox.AddString(m_strName);
	m_ListBox.SetCurSel(m_nIndex);

	if(m_pTypes && m_pTypes->GetSize()){
		for(ptrdiff_t i=0; i<m_pTypes->GetSize(); i++){
			m_ComboBox.AddString(m_pTypes->GetAt(i));
		}
		m_ComboBox.SetCurSel(m_nType);
	}else{
		m_ComboBox.ShowWindow(SW_HIDE);
		m_staticType.ShowWindow(SW_HIDE);
	}

	return 1;  
}

void CDlgAddSubDoc::OnLbnSelchangeList1(){
	int nSel = m_ListBox.GetCurSel();
	if(nSel == -1 || nSel == m_nIndex)return;
	
	m_ListBox.DeleteString(m_nIndex);
	m_nIndex = nSel;
	m_ListBox.InsertString(m_nIndex, m_strName);
	m_ListBox.SetCurSel(m_nIndex);
}

void CDlgAddSubDoc::OnEnChangeEdit1(){
	UpdateData();
	m_ListBox.DeleteString(m_nIndex);
	m_ListBox.InsertString(m_nIndex, m_strName);
	m_ListBox.SetCurSel(m_nIndex);
}


void CDlgAddSubDoc::OnCbnSelchangeCombo1()
{
	// TODO: добавьте свой код обработчика уведомлений
}
