#pragma once
#include "GL3DView.h"

//Представление для CalcStage (вывод расчётов)
class CCalcStage3DView : public CGL3DView
{

	DECLARE_DYNCREATE(CCalcStage3DView);
public:
	CCalcStage3DView(void);
	~CCalcStage3DView(void);
	virtual void SetDocument(CEMMADoc *pDoc);
	//virtual void SetCentr(double x, double y,double z);	
	GLint* RectToViewport(const CRect &rect);

	DECLARE_MESSAGE_MAP()
	afx_msg LRESULT OnDraftRect(WPARAM wParam, LPARAM lParam);
	afx_msg void On2D();
	afx_msg LRESULT OnActivate(WPARAM wParam, LPARAM lParam);

};
