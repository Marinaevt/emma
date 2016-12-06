#pragma once
#include "afxwin.h"
#include "EMMAPrePost.h"
#include "../../Common/DlgShowError.h"

// CDlgNewTask dialog

class CDlgNewTask : public CDialog
{
	DECLARE_DYNAMIC(CDlgNewTask)

public:
	CDlgNewTask(CWnd* pParent = nullptr);  // standard constructor
	virtual ~CDlgNewTask();

// Dialog Data
	enum { IDD = IDD_NEW_TASK };

	size_t GetSel();

protected:
	int m_nSelect;
	//virtual void DoDataExchange(CDataExchange* pDX);  // DDX/DDV support
	void FillList();

	static const int m_nCountTasks = 6;	//Количество Задач для выбора в диалоговом окне
	UNINT m_nTasksID[m_nCountTasks] = { IDS_GEN2D, IDS_FORMING, IDS_ROLLING, IDS_DRAGGING, IDS_GEN3D, IDS_BALLOON };
	UNINT m_nTasksDescript[m_nCountTasks] = { IDS_GEN2D_DESCR, IDS_FORMING_DESCR, IDS_ROLLING_DESCR, IDS_DRAGGING_DESCR, IDS_GEN3D_DESCR, IDS_BALLOON_DESCR };

	DECLARE_MESSAGE_MAP()
	//afx_msg void CDlgNewTask::OnItemchangedList1(NMHDR* pNMHDR, LRESULT* pResult);
	afx_msg void CDlgNewTask::OnItemchangedList1();
public:
	virtual int OnInitDialog();
};

