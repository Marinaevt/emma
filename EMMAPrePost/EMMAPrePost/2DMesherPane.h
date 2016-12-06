#pragma once
#include "emmarightpane.h"
#include "2DMesher.h"
#include "2DSpecimen.h"
#include "Resource.h"
#include "RegisteredMessages.h"
#include "../../Common/DlgShowError.h"


class C2DMeshPane :	public CEMMARightPane
{
public:
	C2DMeshPane(size_t nID): CEMMARightPane(nID){
		m_cursor = 0;
	}
	virtual ~C2DMeshPane(void);

	CButton m_button_contour;

public:

	C2DMesher* GetMesher(){return m_pDoc ? dynamic_cast<C2DMesher*>(m_pDoc) : nullptr; }
	C2DOutline* GetOutline();

	void AdjustLayout();	//! Устанавливаем размеры Тулбара и Таблицы

	DECLARE_MESSAGE_MAP()

	////TOOLBAR
	afx_msg void OnClearMesh();
	afx_msg void OnUpdateClearMesh(CCmdUI *pCmdUI);

	afx_msg void OnNewMesh();
	afx_msg void OnUpdateNewMesh(CCmdUI *pCmdUI);

	afx_msg void OnNewMeshPoint();
	afx_msg void OnUpdateNewMeshPoint(CCmdUI *pCmdUI);

	//afx_msg void OnImportSketch();
	//afx_msg void OnUpdateImportSketch(CCmdUI *pCmdUI);

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

	afx_msg void ShowContours(); 
	afx_msg void OnEnableButton(CCmdUI *pCmdUI);
};
