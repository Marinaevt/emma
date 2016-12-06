#include "StdAfx.h"
#include "3DCalcStagePane.h"


C3DCalcStagePane::~C3DCalcStagePane(void)
{
}

BEGIN_MESSAGE_MAP(C3DCalcStagePane, CEMMARightPane)
	ON_WM_CREATE()
	ON_WM_SIZE()

	ON_NOTIFY(NM_CUSTOMDRAW, ID_ANGLE_SLIDER_3D, &C3DCalcStagePane::OnSliderAngle3D)
	ON_NOTIFY(NM_CUSTOMDRAW, ID_STEP_SLIDER_3D, &C3DCalcStagePane::OnSliderStep3D)

	ON_EN_CHANGE(IDC_EDIT_STEP, &C3DCalcStagePane::OnEditChangeStep)
	ON_EN_CHANGE(IDC_EDIT_ANGLE, &C3DCalcStagePane::OnEditChangeAngle)


	ON_COMMAND(ID_INSTRUMENT_CHKBX, &C3DCalcStagePane::OnDrawInstrumentChanged)
	ON_COMMAND(ID_FILL_COLOR_CHKBX, &C3DCalcStagePane::OnFillColorChanged)

END_MESSAGE_MAP()


int C3DCalcStagePane::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CEMMADockablePane::OnCreate(lpCreateStruct) == -1) {
		return -1;
	}

	CRect rectDummy;
	rectDummy.SetRectEmpty();

	//---------------------------------------------------------------------------------	
	CString str_angleLabel;
	str_angleLabel.LoadStringW(ID_3DCALCSTAGEPANE_ANGLELABEL);
	if (!m_AngleLabel.Create(str_angleLabel, WS_VISIBLE | WS_CHILD, rectDummy, this))
	{
		CDlgShowError cError(ID_ERROR_3DCALCSTAGEPANE_FAIL_ANGLELABEL); //_T("Failed to create Angle Label \n"));
		return -1;      // fail to create
	}

	//EditText
	if (!m_angle3D.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, IDC_EDIT_ANGLE))
	{
		CDlgShowError cError(ID_ERROR_3DCALCSTAGEPANE_FAIL_ANGLETEXT); //_T("Failed to create Angle EditText \n"));
		return -1;      // fail to create
	}
	m_angle3D.SetWindowTextW(_T("5"));

	//Slider
	if (!m_Angle.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, ID_ANGLE_SLIDER_3D))
	{
		CDlgShowError cError(ID_ERROR_3DCALCSTAGEPANE_FAIL_ANGLESLIDER); //_T("Failed to create Angle Slider \n"));
		return -1;      // fail to create
	}

	m_Angle.SetPos(0);
	for (int i = 0; i < 361; i += 10) {
		m_Angle.SetTic(i);
	}

	//---------------------------------------------------------------------------------

	CString str_stepLabel;
	str_stepLabel.LoadStringW(ID_3DCALCSTAGEPANE_STEPLABEL);
	if (!m_StepLabel.Create(str_stepLabel, WS_VISIBLE | WS_CHILD, rectDummy, this))
	{
		CDlgShowError cError(ID_ERROR_3DCALCSTAGEPANE_FAIL_STEPLABEL); //_T("Failed to create Step Label \n"));
		return -1;      // fail to create
	}

	//EditText
	if (!m_step3D.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, IDC_EDIT_STEP))
	{
		CDlgShowError cError(ID_ERROR_3DCALCSTAGEPANE_FAIL_STEPTEXT); //_T("Failed to create Step EditText \n"));
		return -1;      // fail to create
	}
	m_step3D.SetWindowTextW(_T("10"));

	//Slider
	if (!m_Step.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, ID_STEP_SLIDER_3D))
	{
		CDlgShowError cError(ID_ERROR_3DCALCSTAGEPANE_FAIL_STEPSLIDER); //_T("Failed to create Step Slider \n"));
		return -1;      // fail to create
	}

	m_Step.SetPos(0);
	for (int i = 0; i < 361; i += 10) {
		m_Step.SetTic(i);
	}

	//--------------------------------------------------------------------------------------	

	CString str_drawCheckbox;
	str_drawCheckbox.LoadStringW(ID_3DCALCSTAGEPANE_DRAWINSTRUM);
	if (!m_DrawInstrumentCheckBox.Create(str_drawCheckbox, WS_VISIBLE | BS_AUTOCHECKBOX, rectDummy, this, ID_INSTRUMENT_CHKBX))
	{
		CDlgShowError cError(ID_ERROR_3DCALCSTAGEPANE_FAIL_DRAWINSTRUM); //_T("Failed to create Draw Instrument \n"));
		return -1;      // fail to create
	}

	m_DrawInstrumentCheckBox.SetCheck(BST_CHECKED);

	//-----------------------------------------------------------------------------------------   

	CString str_fillColor;
	str_fillColor.LoadStringW(ID_3DCALCSTAGEPANE_FILLCOLOR);
	if (!m_FillColorCheckBox.Create(str_fillColor, WS_VISIBLE | BS_AUTOCHECKBOX, rectDummy, this, ID_FILL_COLOR_CHKBX))
	{
		CDlgShowError cError(ID_ERROR_3DCALCSTAGEPANE_FAIL_FILLCOLOR); //_T("Failed to create Fill Color\n"));
		return -1;      // fail to create
	}

	m_FillColorCheckBox.SetCheck(BST_CHECKED);

	UpdateData();
	AdjustLayout();

	return 0;

}

void C3DCalcStagePane::AdjustLayout(void) {
	CEMMADockablePane::AdjustLayout();

	if (GetSafeHwnd() == nullptr)		return;

	CRect rectClient;
	GetClientRect(rectClient);


	long SliderAngle3DLabelTop = rectClient.top;
	long SliderAngle3DEditTop = SliderAngle3DLabelTop + m_DefaultControlHeight;
	long SliderAngle3DTop = SliderAngle3DEditTop + m_DefaultControlHeight;

	long SliderStep3DLabelTop = SliderAngle3DTop + m_DefaultControlHeight;
	long SliderStep3DEditTop = SliderStep3DLabelTop + m_DefaultControlHeight;
	long SliderStep3DTop = SliderStep3DEditTop + m_DefaultControlHeight;

	long DrawInstrumentCheckBox3DLabelTop = SliderStep3DTop + m_DefaultControlHeight;
	long DrawInstrumentCheckBox3DTop = DrawInstrumentCheckBox3DLabelTop + m_DefaultControlHeight;

	long FillColorCheckBox3DTop = DrawInstrumentCheckBox3DTop;

	m_AngleLabel.SetWindowPos(nullptr, rectClient.left, SliderAngle3DLabelTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	m_angle3D.SetWindowPos(nullptr, rectClient.left, SliderAngle3DEditTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	m_Angle.SetWindowPos(nullptr, rectClient.left, SliderAngle3DTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);

	m_StepLabel.SetWindowPos(nullptr, rectClient.left, SliderStep3DLabelTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	m_step3D.SetWindowPos(nullptr, rectClient.left, SliderStep3DEditTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	m_Step.SetWindowPos(nullptr, rectClient.left, SliderStep3DTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);

	m_FillColorCheckBox.SetWindowPos(nullptr, rectClient.left, FillColorCheckBox3DTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	m_DrawInstrumentCheckBox.SetWindowPos(nullptr, rectClient.left, DrawInstrumentCheckBox3DLabelTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);

	//m_ScaleOpacityLabel.SetWindowPos(nullptr, rectClient.left, scaleOpacityLabelTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	//m_ScaleOpacitySlider.SetWindowPos(nullptr, rectClient.left, scaleOpacitySliderTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);	


}

void C3DCalcStagePane::DoDataExchange(CDataExchange* pDX)
{
	//m_angle = m_AngleSliderPos/ 3.6;
	//  m_step = m_StepSliderPos;

	CEMMADockablePane::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_ANGLE, m_AngleSliderPos);
	DDX_Slider(pDX, ID_ANGLE_SLIDER_3D, m_AngleSliderPos);
	DDV_MinMaxInt(pDX,  /*m_angle*/ m_AngleSliderPos, 0, 360);

	//DDX_Control(pDX, IDC_EDIT_ANGLE, m_angle);
	DDX_Text(pDX, IDC_EDIT_STEP, m_StepSliderPos);
	DDX_Slider(pDX, ID_STEP_SLIDER_3D, m_StepSliderPos);
	DDV_MinMaxInt(pDX,  /* m_step*/  m_StepSliderPos, 0, 360);
	//DDX_Control(pDX, IDC_EDIT_STEP, m_step);
}


void C3DCalcStagePane::OnSliderAngle3D(NMHDR *pNMHDR, LRESULT *pResult)
{
	UpdateData(true);


	CBasicView *view = GetView();
	if (view && view->GetParameter()) {

		view->GetParameter()->SetAngle(m_AngleSliderPos);

		//view->Invalidate();
	}

	*pResult = 0;
}

void C3DCalcStagePane::OnSliderStep3D(NMHDR *pNMHDR, LRESULT *pResult)
{
	UpdateData(true);

	CBasicView *view = GetView();
	if (view && view->GetParameter()) {

		view->GetParameter()->SetStep3D(m_StepSliderPos);

		//view->Invalidate();
	}

	*pResult = 0;
}



void C3DCalcStagePane::OnFillColorChanged()
{
	UpdateFillColor();
}

void C3DCalcStagePane::OnDrawInstrumentChanged()
{
	UpdateDrawInstrument();
}


void C3DCalcStagePane::UpdateDrawInstrument()
{
	//SetDrawInstrument();
	UpdateData(true);

	CBasicView *view = GetView();
	if (view && view->GetParameter()) {

		view->GetParameter()->SetDrawInstrument(m_DrawInstrumentCheckBox.GetCheck() == BST_CHECKED);
	}
}
void C3DCalcStagePane::UpdateFillColor()
{
	UpdateData(true);

	CBasicView *view = GetView();
	if (view && view->GetParameter()) {

		view->GetParameter()->SetFillColor(m_FillColorCheckBox.GetCheck() == BST_CHECKED);
	}
}

void C3DCalcStagePane::OnEditChangeStep()
{
	//UpdateData(true);

	//CBasicView *view = GetView();
	//if (view && view->GetParameter()){
	////m_step=m_StepSliderPos;
	//	view->GetParameter()->SetStep3D(m_StepSliderPos);
	//	//view->Invalidate();
	//
	//}

}

void C3DCalcStagePane::OnEditChangeAngle()
{
	//UpdateData(true);

	//CBasicView *view = GetView();
	//if (view && view->GetParameter()){
	////m_angle=m_AngleSliderPos;
	//	view->GetParameter()->SetEditAngle(m_AngleSliderPos);
	//	//view->Invalidate();
	//}

}

CBasicView* C3DCalcStagePane::GetView()
{
	if (!m_pDoc) return nullptr;

	return dynamic_cast<CBasicView*>(m_pDoc->GetView());
}