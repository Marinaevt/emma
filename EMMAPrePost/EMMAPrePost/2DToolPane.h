#pragma once
#include "emmarightpane.h"
#include "2DTool.h"
#include "Resource.h"
#include "RegisteredMessages.h"


class C2DToolPane :	public CEMMARightPane
{
public:
	C2DToolPane(size_t nID):CEMMARightPane(nID){
		m_cursor = 0;
	}
	virtual ~C2DToolPane(void);

public:
	
	CRightToolBar m_wndToolBar_sketch;	//дополнительный тулбар (2я строка)

	void AdjustLayout();	//! Устанавливаем размеры Тулбара и Таблицы
	
	C2DTool* GetTool(){return m_pDoc ? dynamic_cast<C2DTool*>(m_pDoc) : nullptr;}

	DECLARE_MESSAGE_MAP()

	////TOOLBAR
	afx_msg void OnClearSketch();
	afx_msg void OnUpdateClearSketch(CCmdUI *pCmdUI);

	afx_msg void OnNewNode();
	afx_msg void OnUpdateNewNode(CCmdUI *pCmdUI);

	afx_msg void OnNewCurve();
	afx_msg void OnUpdateNewCurve(CCmdUI *pCmdUI);

	afx_msg void OnNewCircleArc();
	afx_msg void OnUpdateNewCircleArc(CCmdUI *pCmdUI);

	afx_msg void OnNewFacet();
	afx_msg void OnUpdateNewFacet(CCmdUI *pCmdUI);

	afx_msg void OnNewToolNode();
	afx_msg void OnUpdateNewToolNode(CCmdUI *pCmdUI);

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

};
