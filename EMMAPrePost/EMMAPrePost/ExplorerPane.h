#pragma once
#include "EMMADockablePane.h"
#include "../../Common/AllTypes.h"
#include "EMMAPrePost.h"
#include "RegisteredMessages.h"
#include "2DResults.h"
#include "ChildFrm.h"
#include "Gen2DOperation.h"

////////////////////////////////////////////////////////////////////////////////////////////
// CFileViewToolBar
////////////////////////////////////////////////////////////////////////////////////////////
class CExplorerToolBar : public CMFCToolBar
{
	virtual void OnUpdateCmdUI(CFrameWnd* /*pTarget*/, int bDisableIfNoHndler)	{
		CMFCToolBar::OnUpdateCmdUI((CFrameWnd*) GetOwner(), bDisableIfNoHndler);
	}
	virtual int AllowShowOnList() const { return 0; }
};

////////////////////////////////////////////////////////////////////////////////////////////
// CExplorerPane
////////////////////////////////////////////////////////////////////////////////////////////
class CExplorerPane : public CEMMADockablePane
{
	DECLARE_DYNAMIC(CExplorerPane)

public:
	CExplorerPane(size_t nID);
	virtual ~CExplorerPane();

	void FillTree();
	void SetCurSel(HTREEITEM nTreeItem);
	void AdjustLayout();
	CEMMADoc *GetCurDoc();
	void SetDocument(CEMMADoc *pDoc);
	void SetIconsToTree(CString str);
	
	//! Раскрываем все подэлементы
	void ExpandAll();	

	//////Roma
	void UpdatePane();

protected:
	
	CExplorerToolBar m_wndToolBar;
	CTreeCtrl  m_wndTree;
	
	CImageList m_imageList; //image list used by the tree
	int m_icons_name;
	CBitmap m_bitmap; //bitmap witch loads 32bits bitmap

	bool m_bDontSendWM;

protected:
	virtual void OnPressCloseButton();
protected:
	DECLARE_MESSAGE_MAP()
	
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnSize(UNINT nType, int cx, int cy);
	afx_msg void OnPaint();
	
	afx_msg void OnTvnSelchangedEPTree(NMHDR *pNMHDR, LRESULT *pResult);


public:
	afx_msg void OnAdd();
	afx_msg void OnUpdateAdd(CCmdUI *pCmdUI);
	void CExplorerPane::OnRemoveAllSubDoc();
	afx_msg void OnRemove();
	afx_msg void OnUpdateRemove(CCmdUI *pCmdUI);
	afx_msg void OnLoadResults();
	afx_msg void OnUpdateResults(CCmdUI *pCmdUI);
	
};


