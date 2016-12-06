// MainFrm.h : interface of the CMainFrame class
//
#pragma once
#include "EMMADoc.h"
#include "EMMAPrePost.h"
#include "ExplorerPane.h"
#include "ChildFrm.h"
#include "DlgNewTask.h"

#include "Gen3DDoc.h"
#include "Gen2DDoc.h"

#include "./Rolling/RollingDoc.h"
#include "./Forming/2DFormingDoc.h"
#include "./Dragging/2DDraggingDoc.h"
#include "./Balloon/C2DBalloonDoc.h"
#include "RegisteredMessages.h"

#include "2DSketchPane.h"
#include "2DMesherPane.h"
#include "MaterialPane.h"
#include "2DTrajectoryPane.h"
#include "2DToolPane.h"
#include "2DBuilderPane.h"
#include "2DCalcStagePane.h"
#include "3DCalcStagePane.h"


//!  Главное окно приложения 
/*!
  хранит данные приложения в переменной m_pDoc;
  управляет окнами, содержащими представления документа (CChildFrame);
  управляет панелями инструментов, меню, статус баром и другими элементами интерфейса.
*/
class CMainFrame : public CMDIFrameWndEx
{
	DECLARE_DYNAMIC(CMainFrame)
public:
	CMainFrame();
	virtual ~CMainFrame();

protected:

	//! документ - этот объект хранит все данные приложения
	CEMMADoc *m_pDoc;

	//! Список дочерних окон CChildFrame, содержащих различные представления документа
	CList<CChildFrame *> m_CFList;

	//! Все тулбары приложения лежат в этом массиве
	std::vector<CMFCToolBar *> m_ToolBars;

	//! Все плавающие управляющие панели лежат в этом массиве
	std::vector<CEMMADockablePane *> m_DockBars;

	//! Меню
	//CMFCMenuBar       m_wndMenuBar;		//Заменили на Риббон-ленту
	
	//! Строка состояния
	CMFCStatusBar     m_wndStatusBar;
	
	//! ExplorerPane - панель, отображающая структуру документа. Расположена в левой части главного окна.
	CExplorerPane	  m_wndExplorerPane;

	//! Ribbon - лента кнопок
	CMFCRibbonBar     m_wndRibbonBar;
	CMFCRibbonApplicationButton m_MainButton;

public:
	virtual int LoadFrame(UNINT nIDResource, UNLONG dwDefaultStyle = WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE, CWnd* pParentWnd = nullptr, CCreateContext* pContext = nullptr);
	void CloseDocument();
	void UpdateAllViews();

protected:
	//! Уничтожение (освобождение) документа.
	void DestroyDocument();
	//! Закрываем все дочерние фреймы. Вызывается при уничтожении документа.
	void CloseAllChild();
	//! Отображение тулбаров и всех интерфейсных панелей кроме ExplorerPane.
	void ShowControlBars();
	//! Создаем ExplorerPane.
	void CreateExplorerPane();
	//! Создаем новый дочерний фрейм.
	CMDIChildWnd* CreateNewChild(CEMMADoc *pDoc);
	//! Создаем плавающие панели.
	void CreateDockPanes();
	//! Меняем активные категории\вкладки у Ribbon панели
	void SetActiveRibbonCategory();
	//! Получаем текущий активный документ из дерева
	CEMMADoc* GetCurDoc() {
		return m_wndExplorerPane.GetCurDoc();
	}

	DECLARE_MESSAGE_MAP()
	// Generated message map functions
protected:
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnApplicationLook(UNINT id);
	afx_msg void OnUpdateApplicationLook(CCmdUI* pCmdUI);

public:
	afx_msg void OnFileNew();
	afx_msg void OnFileOpen();
	afx_msg void OnFileSave();
	afx_msg void OnUpdateFileSave(CCmdUI* pCmdUI);
	afx_msg void OnFileSaveAs();
	afx_msg void OnUpdateFileSaveAs(CCmdUI* pCmdUI);
	afx_msg void OnFileClose();
	afx_msg void OnUpdateFileClose(CCmdUI* pCmdUI);
	afx_msg LRESULT OnChildFrm(WPARAM wParam, LPARAM lParam);
	afx_msg LRESULT OnChildClose(WPARAM wParam, LPARAM lParam);
	afx_msg LRESULT OnEPSelChange(WPARAM wParam, LPARAM lParam);	//Explorer Pane select change
	afx_msg LRESULT OnEPClose(WPARAM wParam, LPARAM lParam);		//Explorer Pane closed
	afx_msg LRESULT OnSubDocDestroy(WPARAM wParam, LPARAM lParam);
	/////Roma
	afx_msg LRESULT OnUpdatePane(WPARAM wParam, LPARAM lParam);
	
	//Yura
	afx_msg LRESULT OnReopenSaves(WPARAM wParam, LPARAM lParam);
	
	//Ribbon buttons

	//Активирует или деактивирует кнопки Ribbon в зависимости от типа документа
	void RibbonButtonOnOff(CCmdUI *pCmdUI, TYPEID nId);

	//Add element (Operation, Instrument, etc)
	afx_msg void OnExplorerElementAdd();
	afx_msg void OnUpdateExplorerElementAdd(CCmdUI* pCmdUI);
	afx_msg void OnExplorerElementDelete();	// Delete Explorer element
	afx_msg void OnUpdateExplorerElementDelete(CCmdUI* pCmdUI);

	// Cut, Copy, Paste
	int IsEdit(CWnd* pWnd);
	afx_msg void OnEditCut();
	afx_msg void OnEditCopy();
	afx_msg void OnEditPaste();
	// for Cut & Copy
	afx_msg void OnUpdateNeedSel(CCmdUI* pCmdUI);
	// for Paste
	afx_msg void OnUpdateNeedClip(CCmdUI* pCmdUI);
	

	// Категория Редактор (Чертёж, Инструмент, Траектория)
	afx_msg void OnNewNode();
	afx_msg void OnUpdateNewNode(CCmdUI *pCmdUI);
	afx_msg void OnNewCurve();
	afx_msg void OnUpdateNewCurve(CCmdUI *pCmdUI);
	afx_msg void OnNewCircleArc();
	afx_msg void OnUpdateNewCircleArc(CCmdUI *pCmdUI);
	afx_msg void OnNewFacet();
	afx_msg void OnUpdateNewFacet(CCmdUI *pCmdUI);
	afx_msg void OnClearSketch();
	afx_msg void OnUpdateClearSketch(CCmdUI *pCmdUI);
	afx_msg void OnUpdateSetCursorCheck(CCmdUI* pCmdUI, int nType);
	afx_msg void OnUpdateSketchToolTrajectoryDoc(CCmdUI *pCmdUI);	// Метод для включения кнопок Чертежа, Инструмента, Траектории
	afx_msg void OnSaveSketchToMeta();	//Сохранение всего чертежа в meta-файл
	afx_msg void OnEnableButtonSave(CCmdUI *pCmdUI);
	afx_msg void OnSetContours();

	// Категория КЭ сетка
	afx_msg void OnUpdateMesherDoc(CCmdUI *pCmdUI);	// Метод для включения всех кнопок КЭ сетки 
	afx_msg void OnClearMesh();
	//afx_msg void OnUpdateClearMesh(CCmdUI *pCmdUI);
	afx_msg void OnSetMeshType(UNINT id);
	afx_msg void OnUpdateSetMeshTypeMainButton(CCmdUI* pCmdUI);
	afx_msg void OnUpdateSetMeshType(CCmdUI *pCmdUI);
	afx_msg void OnNewMesher();
	afx_msg void OnNewMeshPoint();
	//afx_msg void OnUpdateNewMeshPoint(CCmdUI *pCmdUI);
	afx_msg void OnImportSketch();
	//afx_msg void OnUpdateImportSketch(CCmdUI *pCmdUI);
	afx_msg void OnShowContours();
	//afx_msg void OnUpdateShowContours(CCmdUI *pCmdUI);
	

	// Категория Инструмент
	afx_msg void OnNewToolNode();
	afx_msg void OnUpdateNewToolNode(CCmdUI *pCmdUI);
	
	// Категория Траектория
	afx_msg void OnPlay();
	afx_msg void OnUpdatePlay(CCmdUI *pCmdUI);


	// Категория Проект
	afx_msg void OnBuilderSave();
	afx_msg void OnUpdateBuilderDoc(CCmdUI *pCmdUI);
	
	afx_msg void OnUpdateResultsDoc(CCmdUI *pCmdUI);
	afx_msg void OnLoadResults();
	afx_msg void OnClearResults();
	

public:
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	afx_msg void OnClose();
};


