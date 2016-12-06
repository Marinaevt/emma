#pragma once
#include "GLModelPreview.h"
#include "afxpropertygridctrl.h"
#include "../../../Common/FolderWork.h"
#include "../../../Common/2DOutline.h"
#include "../ImportExport/Parametric/Outline.h"
#include <vector>

using namespace std;
typedef std::vector<std::pair<std::string, double> > CParameters;

//////////////////////////////////////////////////////////////////////////////////////////
//ОПИСАНИЕ КЛАССА
/*! \class CDlgParametricPreview
    \brief Отображает диалог, предназначеный для предпросмотра модели и инициализирует все его элементы 
!*/
class CDlgParametricPreview : public CDialog
{
	DECLARE_DYNAMIC(CDlgParametricPreview)

	/*
	//МЕТОДЫ КЛАССА
	стандартный конструктор CDlgParametricPreview(CWnd* pParent = nullptr);
	виртуальный деструктор virtual ~CDlgParametricPreview();
	поддержка DDX и DDV virtual void DoDataExchange(CDataExchange* pDX);
	инициализация диалога  OnInitDialog();
	метод выбора файла SetCurrentModel();
	при смене выбора файла afx_msg void OnSelChange()
	при смене свойств модели void ModelChange();
	при изменении параметров в таблице вызывается сначала это метод 	afx_msg LRESULT OnPropertyChanged(__in WPARAM wparam, __in LPARAM lparam);
	загрузка списка свойств bool ShowModelProperties();
	удаление старых свойств bool DeleteProperties();
	*/


public:
	CDlgParametricPreview(CWnd* pParent = nullptr);
	virtual ~CDlgParametricPreview();

	enum { IDD = IDD_DlgParametricPreview }; 
	C2DOutline *GetOutline(){return m_pModelOutline;}

protected:
	virtual void DoDataExchange(CDataExchange* pDX); 
	virtual int OnInitDialog();
	bool SetCurrentModel();
	afx_msg void OnSelChange();
	void ModelChange();
	afx_msg LRESULT OnPropertyChanged( __in WPARAM wparam, __in LPARAM lparam );
	bool ShowModelProperties();
	bool DeleteProperties();

	/* 
		Если есть папка по-умолчанию, то смотрим .meta файлы в ней
		Или выбираем любую другую папку
	*/ 
	bool ChooseFolderInPreview();
	
	bool GetListOfMetaAndSetIn();
	afx_msg void ChooseAnotherFolderInPreview();

	DECLARE_MESSAGE_MAP()

protected:

	CGLModelPreview m_ModelPreview;			//экземпляр GL - окна 
	CMFCPropertyGridCtrl m_PropertyGrid;	//экземпляр таблицы свойств
	CParameters m_ParametricParameters;		//контейнер для свойств модели 
	CListBox m_ListBox;				//экземпляр таблицы файлов 
	std::vector<CString> m_FileList;		//массив имен файлов
	CButton m_ButtonChangeFolder;	//кнопка для выбора новой папки и загрузки списка .meta файлов

	int m_NCurrentFile;				//переменная для хранения текущего(выбранного) номера строки(имени файла) 
	C2DOutline *m_pModelOutline;	//переменная для хранения мета - модели
	MM::COutline m_MMOutline;		//переменная для хранения мета-модели в параметрическом виде
	
	CString m_path2folder;			//путь до папки с файлами мета-моделей
};
