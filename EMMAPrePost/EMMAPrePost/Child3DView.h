#pragma once
#include "GL3DView.h"
#include "RegisteredMessages.h"
#include "ChildFrm.h"

class CChild3DView : public CGL3DView
{
DECLARE_DYNCREATE(CChild3DView);
public:
	CChild3DView(void);
	afx_msg void OnLButtonUp(UNINT nFlags, CPoint point);
	virtual ~CChild3DView(void);
	DECLARE_MESSAGE_MAP()
};
