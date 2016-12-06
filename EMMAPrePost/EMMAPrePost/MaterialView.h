#pragma once
#include "BasicView.h"
#include "MaterialDoc.h"

// CMaterialView
#define	CMaterialView_Parent CBasicView	//CBasicView CEditView

class CMaterialView : public CMaterialView_Parent
{
DECLARE_DYNCREATE(CMaterialView);
public:
	CMaterialView(void);
	virtual ~CMaterialView(void);
	virtual void SetDocument(CEMMADoc *pDoc);
	//получаем из документа содержимое файла материала
	void FillEditArea(CEMMADoc *pDoc);
	
	DECLARE_MESSAGE_MAP()
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg LRESULT OnInitialUpdate(WPARAM wParam, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg void OnSize(UNINT nType, int cx, int cy);

public:
	CEdit m_MatEditArea;
};
