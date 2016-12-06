#pragma once
#include "emmadockablepane.h"
#include "emmarightpane.h"
#include "3DCalcStage.h"
#include "CalcStage3DView.h"
#include "Resource.h"
#include "RegisteredMessages.h"
#include <list>
#include <map>

class C3DCalcStagePane : public CEMMARightPane
{
	public:
	C3DCalcStagePane (size_t nID):CEMMARightPane(nID){
//		m_cursor = 0;
		m_AnglePos = 0;
		m_StepPos = 0;
		m_angle = 0;
		m_step = 0;
		m_AngleSliderPos = 0;
		m_StepSliderPos = 0;
	}
	virtual ~C3DCalcStagePane (void);

public:
	static const int m_DefaultControlHeight = 25;

public:

	//void UpdatePane();
	virtual void DoDataExchange(CDataExchange* pDX);
	void AdjustLayout();	//! Устанавливаем размеры Тулбара и Таблицы

	CBasicView* GetView();

	//CString GetString(CEdit &control);
	//CString GetSAngleValue();
	//CString GetStepValue();

	DECLARE_MESSAGE_MAP()

	////TOOLBAR
	afx_msg void OnSliderAngle3D(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnSliderStep3D(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDrawInstrumentChanged();
	afx_msg void OnFillColorChanged();
	afx_msg void OnEditChangeStep();
	afx_msg void OnEditChangeAngle();

protected:

	void UpdateDrawInstrument();
	void UpdateFillColor();

protected:
	CStatic m_StepLabel;
	CStatic m_AngleLabel;
	
	int m_AnglePos;
	int m_StepPos;

	CSliderCtrl m_Angle;	//slider
	CSliderCtrl m_Step;

	int m_AngleSliderPos;
	int m_StepSliderPos;

	CButton m_FillColorCheckBox;
	CButton m_DrawInstrumentCheckBox;
	
	CEdit m_angle3D;
    CEdit m_step3D;
	int m_angle;
    int m_step;
	
};