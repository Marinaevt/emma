#pragma once
#include "EMMADockablePane.h"
#include "Resource.h"

/*
	CRightToolBar
	Класс для Тулбаров у Панелей, находящихся справа.
*/

class CRightToolBar : public CMFCToolBar
{
	
	//Слушаем только команды Панели CEMMARightPane
	virtual void OnUpdateCmdUI(CFrameWnd* pF, int bDisableIfNoHndler) {
		CFrameWnd* pOwner = (CFrameWnd*)GetOwner();
		CMFCToolBar::OnUpdateCmdUI(pOwner, bDisableIfNoHndler);
	}
	//*/

	//Determines whether the toolbar is displayed in the list of toolbars on the Toolbars pane of the Customize dialog box.
	virtual int AllowShowOnList() const { return 0; }	//No
};


/*
	Класс предназначен для Панелей, находящихся справа.
	У них обычно присутствуют Тулбар(ы) и Таблица параметров.
	Список Панелей: 2DSketchPane, 2DMesherPane, 2DToolPane, C2DTrajectoryPane, C2DBuilderPane, C2DCalcStagePane.
*/

class CEMMARightPane : public CEMMADockablePane {
public:
	CEMMARightPane(size_t nID) :CEMMADockablePane(nID) { m_cursor = 0; };
	~CEMMARightPane(void);

	CMFCPropertyGridCtrl m_wndPropList;
	CRightToolBar m_wndToolBar;
	int m_cursor; //тип курсора

	virtual void UpdatePane();		//! Обновляем правую панель с данными
	int OnEraseBkgnd(CDC* pDC);	//! На случай, если OnPaint() не сработает

	void SetColumnNames(void);	//! Устанавливаем заголовки для колонок таблицы

	//
	DECLARE_MESSAGE_MAP()
	afx_msg void OnPaint();
	afx_msg void OnSize(UNINT nType, int cx, int cy);
	afx_msg LRESULT OnPropertyChanged(__in WPARAM wparam, __in LPARAM lparam);
	afx_msg LRESULT OnUpdateRightPanes(WPARAM wParam, LPARAM lParam);
};

