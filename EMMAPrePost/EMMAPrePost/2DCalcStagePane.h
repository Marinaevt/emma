#pragma once
#include "EmmaDockablePane.h"
#include "EmmaRightPane.h"
#include "2DCalcStage.h"
#include "DlgNewScale.h"
#include "../../Common/DlgShowError.h"
#include "CalcStage2DView.h"
#include "Resource.h"
#include "RegisteredMessages.h"
#include <list>
#include <map>

class C2DCalcStagePane : public CEMMARightPane, ICallbackable
{
public:
	C2DCalcStagePane(size_t nID);
	virtual ~C2DCalcStagePane(void);

	C2DCalcStage* GetStage();

	CBasicView* GetView();

	//Получаем путь до шкал
	CString GetScalePath(){ return m_Pane_Path; }
	//Устанавливаем путь до шкал
	void SetScalePath(CString new_scale_path){ m_Pane_Path = new_scale_path; }
	
	//загружаем путь до папки со шкалами
	bool LoadPath();

	void UpdatePane();
	virtual void DoDataExchange(CDataExchange* pDX); 
	void AdjustLayout(void);

	DECLARE_MESSAGE_MAP()

	afx_msg void OnNewScale();

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

	afx_msg void OnScaleCmbChanged();
	afx_msg void OnParamCmbChanged();
	afx_msg void OnDrawLinesChanged();
	afx_msg LRESULT OnPropertyChanged( __in WPARAM wparam, __in LPARAM lparam );
	afx_msg void OnScaleOpacitySliderMove(NMHDR *pNMHDR, LRESULT *pResult);

protected:
	void Callback();
	C2DScale LoadScale(CString name);
	// Загружает шкалы в массив m_Scales
	void LoadAllScales();
	void InitFields();
	void UpdateCalcStage();
	void UpdateParamCmb();
	void UpdateScaleCmb();
	void UpdateScaleGrid();
	void UpdateScaleDrawLines();

public:
	
	//высота контролов
	static const int m_DefaultControlHeight = 25;

protected:
	//Путь до сохранённых шкал
	CString m_Pane_Path;
	C2DCalcStage *m_OldStage;

	CStatic m_ParamLabel;
	CComboBox m_ParamCmb;
	CStatic m_ScaleLabel;
	CComboBox m_ScaleCmb;
	CStatic m_ScaleOpacityLabel;
	CButton m_DrawLinesCheckBox;
	CSliderCtrl m_ScaleOpacitySlider;
	CMFCPropertyGridCtrl m_ScaleGrid;

	std::vector<C2DScale> m_Scales;
	std::map<int, CString> m_FieldMap;
	
	size_t m_SelectedScaleIndex;
	int m_SliderPos;
};
