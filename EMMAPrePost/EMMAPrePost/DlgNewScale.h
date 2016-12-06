#pragma once
#include "Resource.h"
#include "../../Common/2DScale.h"
#include "afxwin.h"

class CDlgNewScale : public CDialog
{
public:
	CDlgNewScale(CWnd* pParent = nullptr);
	CDlgNewScale(CString path, CWnd* pParent = nullptr);
	virtual ~CDlgNewScale(void);

public:
	C2DScale *m_Scale;
	CString m_dlg_path;	//путь до папки со шкалами

	enum { IDD = IDD_NEW_SCALE };

	virtual void DoDataExchange(CDataExchange* pDX); 
	virtual int OnInitDialog();

	CString GetString(CEdit &control);
	CString GetScaleName();

	DECLARE_MESSAGE_MAP()

	afx_msg void OnButtonSaveClick();

protected:

	CEdit m_ScaleName;
	CEdit m_ScaleCount;
	CButton m_IsDiscreteChk;
	CButton m_IsLogChk;

	bool SaveScale();
	void GenerateScale();
	void DrawScale();
	void DrawGradientRect(CRect &r, C2DColor &cLeft, C2DColor &cRight);
};

