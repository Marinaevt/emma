#pragma once
#include "emmarightpane.h"
#include "2DBuilder.h"
#include "RegisteredMessages.h"

class C2DBuilderPane :	public CEMMARightPane
{
	public:
	C2DBuilderPane(size_t nID):CEMMARightPane(nID){
//		m_cursor = 0;
	}
	virtual ~C2DBuilderPane(void);

public:

	C2DTool* GetTool(){return m_pDoc ? dynamic_cast<C2DTool*>(m_pDoc) : nullptr;}

	void AdjustLayout();	//! Устанавливаем размеры Тулбара и Таблицы

	DECLARE_MESSAGE_MAP()

	////TOOLBAR
	afx_msg void OnSave();
	afx_msg void OnUpdateSave(CCmdUI *pCmdUI);

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

};
