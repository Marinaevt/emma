#include "StdAfx.h"
#include "2DCalcStagePane.h"

C2DCalcStagePane::C2DCalcStagePane(size_t nID) :CEMMARightPane(nID) {
	m_Pane_Path = CString(_T(""));
	m_SelectedScaleIndex = 0;
	m_SliderPos = 0;
	m_OldStage = nullptr;
	InitFields();
}

C2DCalcStagePane::~C2DCalcStagePane(void)
{
}

BEGIN_MESSAGE_MAP(C2DCalcStagePane, CEMMARightPane)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_REGISTERED_MESSAGE(AFX_WM_PROPERTY_CHANGED,C2DCalcStagePane::OnPropertyChanged)
	ON_NOTIFY(NM_CUSTOMDRAW, ID_SCALE_OP_SLIDER, &C2DCalcStagePane::OnScaleOpacitySliderMove)

    ON_CBN_SELENDOK(ID_SCALE_CMB, &C2DCalcStagePane::OnScaleCmbChanged)
	ON_CBN_SELENDOK(ID_PARAM_CMB, &C2DCalcStagePane::OnParamCmbChanged)

	ON_COMMAND(ID_NEW_SCALE, &C2DCalcStagePane::OnNewScale)
	ON_COMMAND(ID_LINES_CHKBX, &C2DCalcStagePane::OnDrawLinesChanged)
END_MESSAGE_MAP()

C2DCalcStage* C2DCalcStagePane::GetStage() {
	if (!m_pDoc) return nullptr;

	// Ничего лишнего
	// Из документа Панели получаем Представление, а из него - текущий документ
	CBasicView *pView = dynamic_cast<CBasicView*>(m_pDoc->GetView());

	if (!pView) return nullptr;

	return dynamic_cast<C2DCalcStage*>(pView->GetDocument());
}

CBasicView* C2DCalcStagePane::GetView()
{
	if (!m_pDoc) return nullptr;

	return dynamic_cast<CBasicView*>(m_pDoc->GetView());
}

#pragma region

int C2DCalcStagePane::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CEMMADockablePane::OnCreate(lpCreateStruct) == -1)
		return -1;

	CRect rectDummy;
	rectDummy.SetRectEmpty();

	CString str4lables;	//для надписей (тех, что ниже)
	int nLoadStr = -1;	//для возврата результата загрузки
	nLoadStr = str4lables.LoadStringW(IDS_SCALE_LABEL_SCALE);	//шкалы

	if (!m_ScaleLabel.Create(str4lables, WS_VISIBLE | WS_CHILD | SS_CENTERIMAGE, rectDummy, this))
	{
		CDlgShowError cError(ID_ERROR_2DCALCSTAGEPANE_FAIL_SCALELABEL); //_T("Failed to create ScaleLabel \n"));
		return -1;      // fail to create
	}

	if (!m_ScaleCmb.Create(WS_VISIBLE | WS_CHILD | CBS_DROPDOWNLIST, rectDummy, this, ID_SCALE_CMB))
	{
		CDlgShowError cError(ID_ERROR_2DCALCSTAGEPANE_FAIL_SCALECMB); //_T("Failed to create ScaleCmb \n"));
		return -1;      // fail to create
	}

	nLoadStr = str4lables.LoadStringW(IDS_SCALE_LABEL_FIELDS);	//какие значения выводим (поля)
	if (!m_ParamLabel.Create(str4lables, WS_VISIBLE | WS_CHILD | SS_CENTERIMAGE , rectDummy, this))
	{
		CDlgShowError cError(ID_ERROR_2DCALCSTAGEPANE_FAIL_PARAMLABEL); //_T("Failed to create ParameterLabel \n"));
		return -1;      // fail to create
	}

	if (!m_ParamCmb.Create(WS_VISIBLE | WS_CHILD | CBS_DROPDOWNLIST | CBS_SORT, rectDummy, this, ID_PARAM_CMB))
	{
		CDlgShowError cError(ID_ERROR_2DCALCSTAGEPANE_FAIL_PARAMCMB); //_T("Failed to create ParameterCmb \n"));
		return -1;      // fail to create
	}

	nLoadStr = str4lables.LoadStringW(IDS_SCALE_LABEL_OPACITY);	//прозрачность
	if (!m_ScaleOpacityLabel.Create(str4lables, WS_VISIBLE | WS_CHILD | SS_CENTERIMAGE , rectDummy, this))
	{
		CDlgShowError cError(ID_ERROR_2DCALCSTAGEPANE_FAIL_OPACITYLABEL); //_T("Failed to create ScaleOpacityLabel \n"));
		return -1;      // fail to create
	}

	if (!m_ScaleOpacitySlider.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, ID_SCALE_OP_SLIDER))
	{
		CDlgShowError cError(ID_ERROR_2DCALCSTAGEPANE_FAIL_OPACITYSLIDER); //_T("Failed to create ScaleOpacitySlider \n"));
		return -1;      // fail to create
	}

	m_ScaleOpacitySlider.SetPos(60);	

	for (int i = 0; i < 101; i += 20){
		m_ScaleOpacitySlider.SetTic(i);
	}

	nLoadStr = str4lables.LoadStringW(IDS_SCALE_LABEL_DRAWLINES);	//отрисовка границ
	if (!m_DrawLinesCheckBox.Create(str4lables, WS_VISIBLE | BS_AUTOCHECKBOX, rectDummy, this, ID_LINES_CHKBX))
	{
		CDlgShowError cError(ID_ERROR_2DCALCSTAGEPANE_FAIL_DRAWCHECKBOX); //_T("Failed to create DrawLinesCheckBox \n"));
		return -1;      // fail to create
	}

	m_DrawLinesCheckBox.SetCheck(BST_CHECKED);
	
	if (!m_ScaleGrid.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, ID_SCALE_GRID))
	{
		CDlgShowError cError(ID_ERROR_2DCALCSTAGEPANE_FAIL_SCALEGRID); //_T("Failed to create ScaleGrid \n"));
		return -1;      // fail to create
	}

	//Сетка значения-цвета
	CString str_value, str_color;
	int nLoadStr_val = -1, nLoadStr_color = -1;
	nLoadStr_val = str_value.LoadStringW(IDS_SCALE_GRID_VALUE);
	nLoadStr_color = str_color.LoadStringW(IDS_SCALE_GRID_COLOR);

	m_ScaleGrid.EnableHeaderCtrl(true, str_value, str_color);
	UpdateData();
	AdjustLayout();

	return 0;
}

void C2DCalcStagePane::AdjustLayout(void){
	CEMMADockablePane::AdjustLayout();
	if (GetSafeHwnd() == nullptr)		return;

	CRect rectClient;
	GetClientRect(rectClient);	//получаем габариты клиентского окна

	//TODO: <long> is different between LLP64/LP64 compilers
	long scaleLabelTop = rectClient.top;
	long scaleCmbTop = scaleLabelTop + m_DefaultControlHeight;

	long paramLabelTop = scaleCmbTop + m_DefaultControlHeight;
	long paramCmbTop = paramLabelTop + m_DefaultControlHeight;
	long scaleOpacityLabelTop = paramCmbTop + m_DefaultControlHeight;
	long scaleOpacitySliderTop = scaleOpacityLabelTop + m_DefaultControlHeight;
	long drawLinesCheckBoxTop = scaleOpacitySliderTop + m_DefaultControlHeight;
	long scaleGridTop = drawLinesCheckBoxTop + m_DefaultControlHeight;
	
	m_ScaleLabel.SetWindowPos(nullptr, rectClient.left, scaleLabelTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	m_ScaleCmb.SetWindowPos(nullptr, rectClient.left, scaleCmbTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);

	m_ParamLabel.SetWindowPos(nullptr, rectClient.left, paramLabelTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	m_ParamCmb.SetWindowPos(nullptr, rectClient.left, paramCmbTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	
	m_ScaleOpacityLabel.SetWindowPos(nullptr, rectClient.left, scaleOpacityLabelTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);
	m_ScaleOpacitySlider.SetWindowPos(nullptr, rectClient.left, scaleOpacitySliderTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);	
	m_DrawLinesCheckBox.SetWindowPos(nullptr, rectClient.left, drawLinesCheckBoxTop, rectClient.Width(), m_DefaultControlHeight, SWP_NOACTIVATE | SWP_NOZORDER);

	CRect scaleGridRect = CRect(rectClient.left, scaleGridTop, rectClient.right, rectClient.bottom);
	m_ScaleGrid.SetWindowPos(nullptr, scaleGridRect.left, scaleGridRect.top, scaleGridRect.Width(), scaleGridRect.Height(), SWP_NOACTIVATE | SWP_NOZORDER);
}

void C2DCalcStagePane::InitFields()
{
	m_FieldMap.clear();

	m_FieldMap[int_d] = _T("интенсивность деформации");
	m_FieldMap[int_ds] = _T("итенсивность скоростей деформации");
	m_FieldMap[avg_d] = _T("средние деформации");
	m_FieldMap[avg_ds] = _T("средние скорости деформации");
	m_FieldMap[smo_d] = _T("сглаженные деформации");
	m_FieldMap[smo_ds] = _T("сглаженные скорости деформации");
}

#pragma endregion Initialization

void C2DCalcStagePane::UpdatePane()
{
	UpdateParamCmb();
	UpdateScaleCmb();
	UpdateCalcStage();
	UpdateScaleGrid();
}

void C2DCalcStagePane::LoadAllScales()
{
	m_Scales.clear();
	CFileFind fileFound;

	//если только открыли результат
	if (m_Pane_Path.IsEmpty()) { 
		LoadPath(); 
	}

	// Выделяем все файлы с расширением .emma_scale
	int isFound = fileFound.FindFile(m_Pane_Path + _T("*.emma_scale"));

	while (isFound != 0)
	{
		isFound = fileFound.FindNextFile();									// Выделяем конкретный файл
		CString name = fileFound.GetFileName();								// Получаем его имя
		m_Scales.push_back(LoadScale(name));								// Грузим шкалу по имени
	}
}

C2DScale C2DCalcStagePane::LoadScale(CString name)
{
	if (name.Find(_T(".emma_scale")) == -1){								// Добавим расширение в имя файла, если
		name += _T(".emma_scale");											// его нет
	}

	CFileFind fileFound;
	CString fullName = m_Pane_Path + name;									// Получим полное имя файла (полный путь
	C2DScale scale;															// + имя файла)

	if (fileFound.FindFile(fullName))										// Проверяем есть ли такой файл
	{
		CStorage scaleFile;													// Если есть, то грузим
		scaleFile.Open(fullName, CStorage::modeRead | CStorage::typeBinary);// Если нет, то вернем дефолтную шкалу
		scale.Load(scaleFile);
	}
			
	return scale;
}

bool C2DCalcStagePane::LoadPath(){
	
	//Находим путь до папки с расчётами
	m_Pane_Path = GetStage()->m_2dstage_path();	
	int nPosSlash = m_Pane_Path.ReverseFind( '\\' ); 
	m_Pane_Path.Truncate( nPosSlash+1 );
	m_Pane_Path += "Scales\\";	//Поддиректория для шкал

	//Пытаемся создать директорию для шкал и проверяем существование
	if (!CreateDirectory(m_Pane_Path,nullptr) && !PathIsDirectory(m_Pane_Path)){
		//Вывод сообщения об ошибке создания директории
		CString strName;
		int nNum = strName.LoadString(ID_ERROR_2DCALCSTAGEPANE_NO_NEW_DIR);
		CDlgShowError Diag_err(strName+m_Pane_Path);	//Показываем окно ошибки
		return false;
	}
	return true;
}

void C2DCalcStagePane::UpdateParamCmb()
{
	m_ParamCmb.ResetContent();

	for (auto it = m_FieldMap.begin(), end = m_FieldMap.end(); it != end; ++it) {
		m_ParamCmb.AddString(it->second);
	}

	CBasicView *view = GetView();
	GLParam *parameter;

	if (!view) {
		return;
	}	
	parameter = view->GetParameter();

	if (parameter->GetField() == 0) {
		parameter->SetField(m_FieldMap.begin()->first);
	}
	
	m_ParamCmb.SelectString(-1, m_FieldMap[parameter->GetField()]);
	m_ParamCmb.Invalidate();
}

void C2DCalcStagePane::UpdateScaleCmb()
{
	m_ScaleCmb.ResetContent();
	LoadAllScales();

	//Create new Scale
	CString str_new_scale;
	int nLoadStr = str_new_scale.LoadString(IDS_NEW_SCALE);
	m_ScaleCmb.AddString(str_new_scale);

	for (auto it = m_Scales.begin(), end = m_Scales.end(); it != end; ++it){
		int nResultAdd = m_ScaleCmb.AddString(it->m_Name());
	}

	//if (m_Scales.size() > 0 && m_SelectedScaleIndex >= 0 && m_SelectedScaleIndex < m_Scales.size()) {
	if (m_Scales.size() > 0 && m_SelectedScaleIndex < m_Scales.size()) {
		int nResultSet1 = m_ScaleCmb.SetCurSel(m_SelectedScaleIndex + 1);
	}
	else {
		//m_SelectedScaleIndex = -1;
		m_SelectedScaleIndex = 0;
		int nResultSet2 = m_ScaleCmb.SetCurSel(-1);
	}
	
	m_ScaleCmb.Invalidate();
}

void C2DCalcStagePane::UpdateScaleGrid()
{	
	m_ScaleGrid.RemoveAll();

	//if (m_SelectedScaleIndex == -1)
	if (m_ScaleCmb.GetCurSel() == -1) {
		return;
	}

	C2DScale selectedScale = m_Scales[m_SelectedScaleIndex];

	size_t count = selectedScale.IsDescrete() ? (selectedScale.GetCount() - 1) : selectedScale.GetCount();

	for (size_t i = 1; i < count; i++)
	{
		C2DColor color = selectedScale.GetColor(selectedScale.GetValue(count - i));
		CMFCPropertyGridColorProperty *colorProperty = 
			new CMFCPropertyGridColorProperty(	selectedScale.GetStringValue(count - i),
												RGB(color.GetR() * 255, color.GetG() * 255, color.GetB() * 255),
												nullptr);
		//"Other" button
		CString str_other;
		int nLoadStr_color = str_other.LoadString(ID_BUTTON_OTHER);

		colorProperty->EnableOtherButton(str_other);
		colorProperty->SetData(count - i);
		m_ScaleGrid.AddProperty(colorProperty);
		//if (i == 0) break;
	}

	m_ScaleGrid.Invalidate();
}

void C2DCalcStagePane::UpdateCalcStage()
{
	CBasicView *view = GetView();

	if (!view || !view->GetParameter()) return;
	
	auto parameter = view->GetParameter();

	if (parameter->GetScale() != nullptr) {
		parameter->GetScale()->SetCallback(nullptr);
	}

	//if (m_SelectedScaleIndex != -1)
	if (m_ScaleCmb.GetCurSel() != -1)
	{
		m_Scales[m_SelectedScaleIndex].SetCallback(this);
		parameter->SetScale(m_Scales[m_SelectedScaleIndex]);				// Устанавливаем новую шкалу
		view->Invalidate();													// Перерисовка	
	}
}

void C2DCalcStagePane::UpdateScaleDrawLines()
{
	m_Scales[m_SelectedScaleIndex].SetDrawLines(m_DrawLinesCheckBox.GetCheck() == BST_CHECKED);
}

#pragma region

void C2DCalcStagePane::Callback()
{
	UpdateScaleGrid();
}

void C2DCalcStagePane::DoDataExchange(CDataExchange* pDX)
{
	CEMMADockablePane::DoDataExchange(pDX);
	DDX_Slider(pDX, ID_SCALE_OP_SLIDER, m_SliderPos);
	DDV_MinMaxInt(pDX, m_SliderPos, 0, 100);
}

void C2DCalcStagePane::OnDrawLinesChanged()
{
	UpdateScaleDrawLines();
	UpdateCalcStage();
}

void C2DCalcStagePane::OnNewScale()
{
	bool bLoad = LoadPath();

	CDlgNewScale dlgScale(m_Pane_Path);
	ptrdiff_t result = dlgScale.DoModal();

	UpdateScaleCmb();
	UpdateScaleGrid();

	// Если создали новую шкалу, то выбираем её в m_ScaleCmb
	if (result == IDOK) {
		CString scaleName = dlgScale.GetScaleName();

		for (size_t i = 0; i < m_Scales.size(); i++){
			if (m_Scales[i].m_Name() == scaleName) {
				m_SelectedScaleIndex = i;	//TODO: приведение типов
				int nSet = m_ScaleCmb.SetCurSel(i + 1);
				break;
			}
		}

		UpdateCalcStage();
		UpdateScaleGrid();
	}
}

void C2DCalcStagePane::OnParamCmbChanged()
{
	GLParam *parameter;
	CBasicView *view = GetView();

	if (!view) {
		return;
	}	
	parameter = view->GetParameter();

	CString selectedField;
	m_ParamCmb.GetLBText(m_ParamCmb.GetCurSel(), selectedField);

	for (auto it = m_FieldMap.begin(), end = m_FieldMap.end(); it != end; ++it) {
		if (selectedField == it->second) {
			parameter->SetField(it->first);
			break;
		}
	}
	
	UpdateCalcStage();
}


LRESULT C2DCalcStagePane::OnPropertyChanged(__in WPARAM wparam, __in LPARAM lparam)
{
	CMFCPropertyGridColorProperty * colorProperty = ( CMFCPropertyGridColorProperty * ) lparam;
	UNLONG colorIndex = colorProperty->GetData();
	C2DColor color = C2DColor(colorProperty->GetColor());
	m_Scales[m_SelectedScaleIndex].SetColor(colorIndex, color);
	bool bSave = m_Scales[m_SelectedScaleIndex].SaveToPath(m_Pane_Path);

	UpdateCalcStage();

	return 0;
}

void C2DCalcStagePane::OnScaleCmbChanged()
{
	int curSel = m_ScaleCmb.GetCurSel();

	if (curSel == 0) {		
		OnNewScale();
		return;
	}

	m_SelectedScaleIndex = curSel - 1;	

	UpdateScaleDrawLines();
	UpdateCalcStage();
	UpdateScaleGrid();
}

void C2DCalcStagePane::OnScaleOpacitySliderMove(NMHDR *pNMHDR, LRESULT *pResult)
{
	UpdateData(true);

	CBasicView *view = GetView();
	if (view && view->GetParameter()){
		view->GetParameter()->m_ScaleVisualizer.SetTrasparency(m_SliderPos);
		view->Invalidate();
	}

	*pResult = 0;
}

#pragma endregion Events
