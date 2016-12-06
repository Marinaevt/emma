#pragma once
#include "../../Common/AllTypes.h"

const UNINT WMR_INITIALUPDATE = RegisterWindowMessage(_T("WMR_INITIALUPDATE"));
const UNINT WMR_DRAFTRECT = RegisterWindowMessage(_T("WMR_DRAFTRECT"));
const UNINT WMR_SUBDOCDESTROY = RegisterWindowMessage(_T("WMR_SUBDOCDESTROY"));
const UNINT WMR_CHILDACTIVATE = RegisterWindowMessage(_T("WMR_CHILDACTIVATE"));
const UNINT WMR_CHILDCLOSE = RegisterWindowMessage(_T("WMR_CHILDCLOSE"));
const UNINT WMR_EP_SELCHANGE = RegisterWindowMessage(_T("WMR_EP_SELCHANGE"));
const UNINT WMR_EP_CLOSE = RegisterWindowMessage(_T("WMR_EP_CLOSE"));
const UNINT WMR_ACTIVATE = RegisterWindowMessage(_T("WMR_ACTIVATE"));


// Roma
const UNINT WMR_UPDATEPANE = RegisterWindowMessage(_T("WMR_UPDATEPANE"));
const UNINT WMR_SETCURSOR = RegisterWindowMessage(_T("WMR_SETCURSOR"));

// Yura
const UNINT WMR_REOPEN_SAVES = RegisterWindowMessage(_T("WMR_REOPEN_SAVES"));