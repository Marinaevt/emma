#pragma once
#include "emmarightpane.h"
#include "MaterialDoc.h"
#include "Resource.h"
#include "RegisteredMessages.h"
#include "../../Common/DlgShowError.h"


class CMaterialPane : public CEMMARightPane
{
public:
	CMaterialPane(size_t nID) : CEMMARightPane(nID){
		m_cursor = 0;
	}
	virtual ~CMaterialPane(void);

protected:
	CButton m_button_contour;
	CButton m_button_save_inside;
	CButton m_button_save_to_meta;
	CStatic m_Label;
	CStatic m_Label_ru;

public:

	CMaterialPane* GetMaterialPane(){
		return m_pDoc ? dynamic_cast<CMaterialPane*> (m_pDoc) : nullptr;
	}
	CMaterialDoc* GetMaterialDoc() {
		return m_pDoc ? dynamic_cast<CMaterialDoc *>(m_pDoc) : nullptr;
	}
	//! Устанавливаем размеры Тулбара и Таблицы
	void AdjustLayout();

	DECLARE_MESSAGE_MAP()

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

	//Button
	afx_msg void OnEnableButton(CCmdUI *pCmdUI);
	//Обновляем текст (если случайно отредактировали)
	afx_msg void RefreshText();
	//Сохранение нового материала в строку IO:CMaterial
	afx_msg void SaveTextInside();
	//Сохранение нового материала в файл
	afx_msg void SaveTextToFile();
	//afx_msg void OnEnableButtonSave(CCmdUI *pCmdUI);
};
