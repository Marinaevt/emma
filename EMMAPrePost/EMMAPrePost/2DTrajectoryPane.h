#pragma once
#include "emmarightpane.h"
#include "2DTrajectory.h"
#include "Resource.h"
#include "RegisteredMessages.h"


class C2DTrajectoryPane :	public CEMMARightPane
{
	public:
	C2DTrajectoryPane(size_t nID):CEMMARightPane(nID){
		m_cursor = 0;
	}
	virtual ~C2DTrajectoryPane(void);

public:

	CRightToolBar m_wndToolBar_sketch;	//дополнительный тулбар (2я строка)

	C2DTrajectory* GetTrajectory(){
		return m_pDoc ? dynamic_cast<C2DTrajectory*>(m_pDoc) : nullptr;
	}

	void AdjustLayout();	//! Устанавливаем размеры Тулбара и Таблицы

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

	afx_msg void OnPlay();
	afx_msg void OnUpdatePlay(CCmdUI *pCmdUI);

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

};
