#pragma once
#include "stdafx.h"
#include "DlgParametricPreview.h"
#include "afxdialogex.h"

IMPLEMENT_DYNAMIC(CDlgParametricPreview, CDialog)
/////////////////////////////////////////////////////////
CDlgParametricPreview::CDlgParametricPreview(CWnd* pParent)
	: CDialog(CDlgParametricPreview::IDD, pParent)
{
	m_pModelOutline = nullptr;
	m_NCurrentFile = 0;
	m_path2folder.Empty();
}

CDlgParametricPreview::~CDlgParametricPreview()
{
	//delete m_pModelOutline;
}

////////////////////////////////////////////////////////
int CDlgParametricPreview::OnInitDialog()
{

	CDialog::OnInitDialog(); //вызываем функцию инициализации от родителя
	CRect rectWindow; //чтобы можно было обращаться к размерам родительского окна
	GetClientRect(&rectWindow);

	//ИНИЦИАЛИЗАЦИЯ ОКНА ПРЕВЬЮ МОДЕЛИ

	if (m_ModelPreview.Create(nullptr,
		_T("Предпросмотр"), WS_BORDER | WS_CHILD | WS_VISIBLE,
		CRect(10, 10, rectWindow.Height() - 25, rectWindow.Height() - 35),
		this, 1,
		nullptr) == -1) {
		return 0;
	}

	//ИНИЦИАЛИЗАЦИЯ ТАБЛИЦЫ СВОЙСТВ
	m_PropertyGrid.Create(WS_BORDER | WS_CHILD | WS_VISIBLE,
		CRect(rectWindow.Width() - 300, 10, rectWindow.Width() - 200, rectWindow.Height() - 35), 
		this, 200234);
	
	//ИНИЦИАЛИЗАЦИЯ ТАБЛИЦЫ ФАЙЛОВ
	m_ListBox.Create(WS_BORDER | WS_CHILD | WS_VISIBLE | WS_VSCROLL | LBS_NOTIFY,
		//позиция X1, позиция Y1 (левая верхняя точка), позиция X2, Y2
		CRect(rectWindow.Width() - 180, 10, rectWindow.Width(), rectWindow.Height() - 60),
		this, 3);

	//устанавливаем такой же шрифт, как в таблице свойств
	m_ListBox.SetFont(m_PropertyGrid.GetFont());

	//ИНИЦИАЛИЗАЦИЯ КНОПКИ СМЕНЫ ДИРЕКТОРИИ
	m_ButtonChangeFolder.Create(_T("Изменить папку"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON,
		CRect(rectWindow.Width() - 180, rectWindow.Height() - 65, rectWindow.Width(), rectWindow.Height() - 40),
		this, 123);
	
	//устанавливаем такой же шрифт, как в таблице свойств
	m_ButtonChangeFolder.SetFont(m_PropertyGrid.GetFont());

	//устанавливаем директорию с .meta файлами
	bool bIsChoosedFolder = ChooseFolderInPreview();
	if (!bIsChoosedFolder) {
		return 0;
	}
	
	//Получаем список файлов из папки и устанавливаем их в окно
	bool bGet = GetListOfMetaAndSetIn();
	if (bGet) {
		return 1;
	}
	
	return 0;
}

bool CDlgParametricPreview::GetListOfMetaAndSetIn() {
	//получаем список файлов из папки
	m_FileList = GetFolderListing(m_path2folder, _T("*.meta"), false);

	int iAddString = 0;
	m_ListBox.ResetContent();	//очищаем всё, что было в списке
	for (size_t i = 0; i < m_FileList.size(); i++)
	{
		// выводим список файлов в таблицу справа
		iAddString = m_ListBox.AddString(m_FileList[i]);
	}

	//обнуляем номер файла, т.к. содержимое может измениться
	m_NCurrentFile = 0;
	//по умолчанию выбираем первый в списке файл
	int iSetCurrent = m_ListBox.SetCurSel(m_NCurrentFile);

	//и сразу устанавливаем модель из него
	bool bIsSetCurrentModel = SetCurrentModel();
	if (!bIsSetCurrentModel) {
		return false;
	}

	//отображаем её в превью окне слева
	bool bIsSetModelPreview = m_ModelPreview.SetModel(m_pModelOutline);
	if (!bIsSetModelPreview) {
		return false;
	}

	return true;
}

bool CDlgParametricPreview::ChooseFolderInPreview(){
	
	//Если папка по-умолчанию есть, то смотрим в неё
	CString path = _T("..\\Data\\Test\\");
	if (DirectoryExists(path)) {
		m_path2folder = path;
		return true;
	}

	//Диалоговое окно для выбора папки с .meta файлами
	ChooseAnotherFolderInPreview();

	return true;
}

afx_msg void CDlgParametricPreview::ChooseAnotherFolderInPreview() {
	CFileDialog fileDlg(1, _T("meta"), 0, 4 | 2, _T("sketch file (*.meta)|*.meta|"));
	if (fileDlg.DoModal() == IDOK) {
		m_path2folder = fileDlg.GetFolderPath();	//Получаем папку
	}
	//обновляем список файлов
	GetListOfMetaAndSetIn();
}

void CDlgParametricPreview::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


//////////////////////////////////////////////////////
//Установка модели из файла
bool CDlgParametricPreview::SetCurrentModel()
{
	CString strFN(m_path2folder);
	strFN += _T("\\") + m_FileList[m_NCurrentFile];

	//m_MMOutline.ClearAll();	//очищаем всё, что было
	if (m_MMOutline.LoadFile(strFN) == xtd::no_good)
	{
		//delete m_pModelOutline;
		//m_MMOutline.ClearAll();
		//m_pModelOutline = nullptr;
		return false;
	} else {
		//delete m_pModelOutline;
		//получаем параметры из уже сохранённых 
		m_MMOutline.GetParameters(m_ParametricParameters);
		m_pModelOutline = m_MMOutline.Instantiate();
		
		bool bDeleteProp = DeleteProperties();
		bool bShowModelProp = ShowModelProperties();
	}

	return true;
}

//////////////////////////////////////////////////////
void CDlgParametricPreview::ModelChange()
{
	m_MMOutline.SetParametersValues(m_ParametricParameters); //передаем обновленные параметры в параметрическую модель

	delete m_pModelOutline;
	m_pModelOutline = m_MMOutline.Instantiate(); //преобразовываем параметрическую модель в обычную
	bool bIsSetModelPreview = m_ModelPreview.SetModel(m_pModelOutline); //устанавливаем модель в GL окно
}


//////////////////////////////////////////////////////
bool CDlgParametricPreview::ShowModelProperties()
{
	m_PropertyGrid.EnableHeaderCtrl(0);
	CMFCPropertyGridProperty* pProp;

	for (size_t i = 0; i < m_ParametricParameters.size(); i++) //запускаем цикл вывода параметров модели в таблицу
	{
		CStringW t_string;
		switch (m_ParametricParameters[i].first[0])
		{
		case 'h':
			if (m_ParametricParameters[i].first.size() == 2)
				t_string = L"Высота выступа";
			else
				t_string = L"Высота";
			break;
		case 'w':
			if (m_ParametricParameters[i].first.size() == 2)
				t_string = L"Длина выступа";
			else
				t_string = L"Длина";
			break;
		case 'r':
			t_string = L"Радиус";
			break;
		default:
			t_string = L"Неизвестный параметр";
			break;
		}
		pProp = new CMFCPropertyGridProperty(t_string, (m_ParametricParameters[i].second), _T("Определяет параметры модели"));
		pProp->SetData(i);

		m_PropertyGrid.AddProperty(pProp); //добавляем свойство в таблицу
	}

	return true;
}



//////////////////////////////////////////////////////
bool CDlgParametricPreview::DeleteProperties()
{
	m_PropertyGrid.RemoveAll(); //удаление старых свойств из таблицы

	return true;
}
//////////////////////////////////////////////////////
afx_msg void CDlgParametricPreview::OnSelChange() //!Смена модели
{
	m_NCurrentFile = m_ListBox.GetCurSel(); //изменение номера выбранной модели
											//переменные для отладки

	bool bIsSetCurrentModel = SetCurrentModel(); //установка мета-модели
	bool bIsSetModelPreview = false;
	if (bIsSetCurrentModel) {
		bIsSetModelPreview = m_ModelPreview.SetModel(m_pModelOutline); //устанавливаем модель в GL окно
	}
	else {
		m_ModelPreview.ClearModel();
	}

	//ModelChange();
}

//////////////////////////////////////////////////////
LRESULT CDlgParametricPreview::OnPropertyChanged(__in WPARAM wparam, __in LPARAM lparam)
{
	CMFCPropertyGridProperty * pProperty = (CMFCPropertyGridProperty *)lparam;
	COleVariant v = pProperty->GetValue();
	UNLONG id = pProperty->GetData();

	v.ChangeType(VT_R8); //меняем насильно тип на double
	m_ParametricParameters[id].second = v.dblVal; //меняем параметр в массиве параметров
	ModelChange(); //исправляем модель при смене параметра
	return (0);
}


BEGIN_MESSAGE_MAP(CDlgParametricPreview, CDialog)
	ON_LBN_SELCHANGE(3, CDlgParametricPreview::OnSelChange)
	ON_REGISTERED_MESSAGE(AFX_WM_PROPERTY_CHANGED, CDlgParametricPreview::OnPropertyChanged)
	ON_COMMAND(123, CDlgParametricPreview::ChooseAnotherFolderInPreview)
END_MESSAGE_MAP()