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

	CMFCPropertyGridCtrl m_wndPropList;	//таблица свойств
	CRightToolBar m_wndToolBar;
	int m_cursor; //тип курсора
	CButton m_buttonApply;	//применение изменений после обновление значений по кнопке

	virtual void UpdatePane();		//! Обновляем правую панель с данными
	

	void SetColumnNames(void);	//! Устанавливаем заголовки для колонок таблицы
	int CreateApplyButton(CRect rectDummy);	//! Создаём кнопку для таблицы свойств

	//
	DECLARE_MESSAGE_MAP()
	afx_msg void OnPaint();
	afx_msg void OnSize(UNINT nType, int cx, int cy);
	afx_msg int OnEraseBkgnd(CDC* pDC);	

	afx_msg void OnEnableButtonApply(CCmdUI *pCmdUI);	//для всех панелей
	

};

