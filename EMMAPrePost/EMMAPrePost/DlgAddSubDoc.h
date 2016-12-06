#pragma once
#include "afxwin.h"
#include "EMMAPrePost.h"

// CDlgAddSubDoc dialog


typedef struct DASD_INITSTRUCT{
	CStringArray *pList;
	CStringArray *pTypes;
	CString *pDefaultName;
	size_t nIndex;
}DASD_INITSTRUCT;

class CDlgAddSubDoc : public CDialog
{
	DECLARE_DYNAMIC(CDlgAddSubDoc)

public:
	ptrdiff_t m_nIndex;
	CString m_strName;
	CStringArray *m_pList;
	CStringArray *m_pTypes;

public:
	CDlgAddSubDoc(const DASD_INITSTRUCT &IS, CWnd* pParent = nullptr);   // standard constructor
	virtual ~CDlgAddSubDoc();

// Dialog Data
	enum { IDD = IDD_ADDSUBDOC };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual int OnInitDialog();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnLbnSelchangeList1();
	CListBox m_ListBox;
	CComboBox m_ComboBox;
	int m_nType;
	CStatic m_staticType;
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnCbnSelchangeCombo1();
};
