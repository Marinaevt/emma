#include "StdAfx.h"
#include "MaterialPane.h"


CMaterialPane::~CMaterialPane(void)
{
}

BEGIN_MESSAGE_MAP(CMaterialPane, CEMMARightPane)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_WM_PAINT()

	//Buttons
	ON_BN_CLICKED(ID_MATERIALPANE_BUTTON_REFRESH, CMaterialPane::RefreshText)
	ON_UPDATE_COMMAND_UI(ID_MATERIALPANE_BUTTON_REFRESH, CMaterialPane::OnEnableButton)

	ON_BN_CLICKED(ID_MATERIALPANE_BUTTON_FIX_CHANGES, CMaterialPane::SaveTextInside)
	ON_UPDATE_COMMAND_UI(ID_MATERIALPANE_BUTTON_FIX_CHANGES, CMaterialPane::OnEnableButton)

	ON_BN_CLICKED(ID_MATERIALPANE_BUTTON_SAVE2OTHER, CMaterialPane::SaveTextToFile)
	ON_UPDATE_COMMAND_UI(ID_MATERIALPANE_BUTTON_SAVE2OTHER, CMaterialPane::OnEnableButton)

END_MESSAGE_MAP()

int CMaterialPane::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CEMMADockablePane::OnCreate(lpCreateStruct) == -1)
		return -1;

	CRect rectDummy;
	rectDummy.SetRectEmpty();

	//сначала создаём текст
	CString str4Label;
	str4Label.LoadStringW(ID_MATERIALPANE_LABELTEXT);

	if (!m_Label.Create(str4Label, WS_VISIBLE | WS_CHILD, rectDummy, this)) {
		CDlgShowError cError(ID_ERROR_MATERIALPANE_LABEL_FAIL); //_T("Failed to create Label \n"));
		return -1;      // fail to create
	};

	//сначала создаём текст
	CString str4LabelRuCodec;
	str4LabelRuCodec.LoadStringW(ID_MATERIALPANE_LABELRUTEXT);

	if (!m_Label_ru.Create(str4LabelRuCodec, WS_VISIBLE | WS_CHILD, rectDummy, this)) {
		CDlgShowError cError(ID_ERROR_MATERIALPANE_LABEL_RU_FAIL); //_T("Failed to create Label \n"));
		return -1;      // fail to create
	};

	//потом создаём кнопки
	CString str4ButtonRefreshText;
	str4ButtonRefreshText.LoadStringW(ID_MATERIALPANE_BUTTON_REFRESH);
	if (!m_button_contour.Create(str4ButtonRefreshText, WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON, rectDummy, this, ID_MATERIALPANE_BUTTON_REFRESH)) {
		CDlgShowError cError(ID_ERROR_MATERIALPANE_BUTTON_REFRESH_FAIL); //_T("Failed to create Button \n"));
		return -1;      // fail to create
	}
	
	CString str4ButtonFixChanges;
	str4ButtonFixChanges.LoadStringW(ID_MATERIALPANE_BUTTON_FIX_CHANGES);
	if (!m_button_save_inside.Create(str4ButtonFixChanges, WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON, rectDummy, this, ID_MATERIALPANE_BUTTON_FIX_CHANGES)) {
		CDlgShowError cError(ID_ERROR_MATERIALPANE_BUTTON_FIX_FAIL); //_T("Failed to create Button \n"));
		return -1;      // fail to create
	}

	CString str4ButtonSaveToOther;
	str4ButtonSaveToOther.LoadStringW(ID_MATERIALPANE_BUTTON_SAVE2OTHER);
	if (!m_button_save_to_meta.Create(str4ButtonSaveToOther, WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON, rectDummy, this, ID_MATERIALPANE_BUTTON_SAVE2OTHER)) {
		CDlgShowError cError(ID_ERROR_MATERIALPANE_BUTTON_SAVE2OTHER); //_T("Failed to create Button \n"));
		return -1;      // fail to create
	}
	
	//! Устанавливаем размеры Тулбара и Таблицы и Кнопок
	AdjustLayout();

	m_cursor = 0;

	return 0;
}

//! Устанавливаем размеры Тулбара и Таблицы и Кнопок
void CMaterialPane::AdjustLayout(void){
	if (GetSafeHwnd() == nullptr) {
		return;
	}

	CRect rectClient;
	GetClientRect(rectClient);

	int cyTlb = 25;

	//устанавливаем размеры Кнопки
	int cyPropList = rectClient.top + cyTlb;
	int nFromLeft = rectClient.left + int(rectClient.Width() / 10.0);
	int nFromRight = rectClient.Width() - int(rectClient.Width() / 10.0)*2;	//обязательное умножение, чтобы был отступ справа
	m_button_contour.SetWindowPos(nullptr, nFromLeft, cyPropList, nFromRight, 20, SWP_NOOWNERZORDER | SWP_SHOWWINDOW);
	cyPropList += cyTlb;	//увеличиваем отступ
	m_button_save_inside.SetWindowPos(nullptr, nFromLeft, cyPropList, nFromRight, 20, SWP_NOOWNERZORDER | SWP_SHOWWINDOW);
	cyPropList += cyTlb;	//увеличиваем отступ
	m_button_save_to_meta.SetWindowPos(nullptr, nFromLeft, cyPropList, nFromRight, 20, SWP_NOOWNERZORDER | SWP_SHOWWINDOW);
	cyPropList += cyTlb;	//увеличиваем отступ

	//устанавливаем размер Текста
	m_Label.SetWindowPos(nullptr, rectClient.left, cyPropList, rectClient.Width(), cyTlb*4, SWP_NOACTIVATE | SWP_NOZORDER);
	cyPropList += cyTlb*5;	//увеличиваем отступ
	m_Label_ru.SetWindowPos(nullptr, rectClient.left, cyPropList, rectClient.Width(), cyTlb * 4, SWP_NOACTIVATE | SWP_NOZORDER);
}


/* BUTTONS */

void CMaterialPane::OnEnableButton(CCmdUI *pCmdUI) {
	//Приходится активировать кнопку для нажатий (MFC style)
	pCmdUI->Enable();
}

void CMaterialPane::RefreshText() {
	
	CString sText;
	CMaterialDoc *pMat = GetMaterialDoc();	//получаем ссылку на документ (Материал)
	if (pMat) {
		CMaterialView* pView = dynamic_cast<CMaterialView*>(pMat->GetView());
		if (pView) {
			//Обновляем текстовое поле содержимым сохранённого файла материала
			sText = pMat->GetMaterialFileContent();
			sText.Replace(_T("\r"), _T(""));
			sText.Replace(_T("\n"), _T("\r\n"));
			pView->m_MatEditArea.SetWindowText(sText);
		}
	}
	//*/
}

//Сохранение всего материала в .spf файл
void CMaterialPane::SaveTextToFile() {
	CFileDialog fileDlg(0, _T("spf"), 0, 4 | 2, _T("material file (*.spf)|*.spf|"));
	if (fileDlg.DoModal() == IDOK) {
		CString fpath = fileDlg.GetPathName();

		CStorage file;
		if (file.Open(fpath, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeText | CStorage::shareDenyNone)) {
			
			CString sTextToFile;
			//получаем ссылку на документ (Материал)
			CMaterialDoc *pMat = GetMaterialDoc();
			if (pMat) {
				CMaterialView* pView = dynamic_cast<CMaterialView*>(pMat->GetView());
				if (pView) {
					//Получаем содержимое текстового поля в строку
					pView->m_MatEditArea.GetWindowText(sTextToFile);
				}
			}
			sTextToFile.Replace(_T("\r"), _T(""));
			file.WriteString(sTextToFile);
			file.Close();
		}
		else {
			//не смогли открыть файл
			CDlgShowError cError(ID_ERROR_MATERIALPANE_OPEN_FILE_FAIL); //_T("Failed to open Material File to save \n"));
		}
	}
}

//Запоминание изменений в IO:CMaterial
void CMaterialPane::SaveTextInside() {
	
	CString sTextToString;
	CMaterialDoc *pMat = GetMaterialDoc();	//получаем ссылку на документ (Материал)
	if (pMat) {
		CMaterialView* pView = dynamic_cast<CMaterialView*>(pMat->GetView());
		if (pView) {
			//Обновляем материал изменённым содержимым текстовой области
			pView->m_MatEditArea.GetWindowText(sTextToString);
			pMat->SetMaterialFileContent(sTextToString);
			//pView->m_MatEditArea.SetWindowText(sTextToString);

			CString fpath = pMat->GetMaterialFileName();
			CStorage file;
			CFileException e;

			/* позволяет сохранять текст с русской кодировкой, использует поток
			sTextToString.Replace(_T("\r"), _T(""));
			bool bSave = SaveWideChar(fpath, sTextToString);

			if(!bSave) {
				//не смогли открыть файл / открыть поток
				//CDlgShowError cError(_T("Failed to open Material File to save \n"));
			}
			//*/

			
			if (file.Open(fpath, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeText | CStorage::shareDenyNone, &e)) {

				sTextToString.Replace(_T("\r"), _T(""));
				//sTextToString += _T("\0");
				/* не работает, пишет через пробел и в 2 раза меньше
				int nTextLen = sTextToString.GetLength();
				file.Write(sTextToString.GetBuffer(nTextLen), nTextLen );
				sTextToString.ReleaseBuffer();
				//*/

				file.WriteString(sTextToString);
				file.Close();
			}
			else {
				//не смогли открыть файл
				CString strError;
				strError.LoadStringW(ID_ERROR_MATERIALPANE_OPEN_FILE_FAIL);
				CDlgShowError cError(strError + AllToString(e.m_cause));
			} 
			//*/
		}
	}
}