// Child2DView.h : interface of the CChild2DView class
//

#pragma once

#include "GL2DView.h"
//#include "GL3DView.h"
#include "RegisteredMessages.h"
#include "ParametricPreview/DlgParametricPreview.h"
#include "ImportExport/DXF/CGraph.h"
#include "ImportExport/Parametric/Outline.h"
#include "resource.h"
#include "EMMAPrePost.h"
#include "../../Common/AllTypes.h"
//#include "ChildFrm.h"

//#define	CChild2DView_Parent CGL2DView		//ненужный define, нигде кроме этого файла

class CChild2DView : public CGL2DView
{
	DECLARE_DYNCREATE(CChild2DView);
// Construction
public:
	CChild2DView();

public:
	virtual void SetDocument(CEMMADoc *pDoc){
		CGL2DView::SetDocument(pDoc);
		if(m_pDoc)m_pI2DSketch = (I2DSketch *)m_pDoc->GetInterface(IID_I2DSKETCH);
	}

// Overrides
protected:
	virtual int PreCreateWindow(CREATESTRUCT& cs);

// Implementation
public:
	virtual ~CChild2DView();

	// Generated message map functions
protected:
	virtual void DrawScene(void);
	C2DOutline *Get2DOutline(){return m_pI2DSketch ? m_pI2DSketch->GetOutline() : nullptr;}

	void Set2DOutline(C2DOutline *pOutline){ if(m_pI2DSketch) m_pI2DSketch->SetOutline(pOutline);}

	void SetClearOutline(){ if(m_pI2DSketch) m_pI2DSketch->SetClearOutline();}

protected:
	I2DSketch *m_pI2DSketch;

DECLARE_MESSAGE_MAP()
	afx_msg LRESULT OnDraftRect(WPARAM wParam, LPARAM lParam);
	afx_msg void OnImportPrm();
	afx_msg void OnUpdateImportPrm(CCmdUI *pCmdUI);
	afx_msg void OnImportDxf();
	afx_msg void OnUpdateImportDxf(CCmdUI *pCmdUI);

	afx_msg void OnTimer(size_t nIDEvent);		// функция реакции на таймер
	void StopTimerPlay(size_t nIDEvent); //функция остановки таймера и анимации

////////////////////////////////// Roma
public:
	//При движении мыши различные действия (тип курсора)
	afx_msg void OnMouseMove(UNINT nFlags, CPoint point);
	//После нажатия мыши различные действия (тип курсора)
	afx_msg void OnLButtonUp(UNINT nFlags, CPoint point);

	afx_msg LRESULT OnSetCursor(WPARAM wParam, LPARAM lParam);
	
};

