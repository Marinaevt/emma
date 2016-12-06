#include "StdAfx.h"
#include "Results.h"


IOIMPL (CResults, CRESULTS)

////////////////////////////////////////////////////////////////////////////////
CResults::CResults(void) {
	m_p2DChildViewDescriptor = nullptr;
	m_p3DChildViewDescriptor = nullptr;
	m_b3D = false;// для отображения 3d меняем значение на true
}

CResults::~CResults(void){
	DestroyViewDescriptor(m_p2DChildViewDescriptor);
	DestroyViewDescriptor(m_p3DChildViewDescriptor);
}

void CResults::DestroyViewDescriptor(CViewDescriptor *pDesctiptor) {
	if (pDesctiptor) {
		/*
		CWnd *pView = pDesctiptor->GetView();
		if (pView && pView != nullptr) {
			pView->GetTopLevelParent()->SendMessage(WMR_SUBDOCDESTROY, (WPARAM)this);
			//pView->GetTopLevelFrame()->SendMessage(WMR_SUBDOCDESTROY, (WPARAM)this);
		} //*/
		delete pDesctiptor;
		pDesctiptor = nullptr;
	}
}

void CResults::Init2DChildViewDescriptor(CRuntimeClass *pViewClass, UNINT nID) {
	if (m_p2DChildViewDescriptor) {

		if (pViewClass != m_p2DChildViewDescriptor->GetViewClass() || nID != m_p2DChildViewDescriptor->GetResID()) {
			CDlgShowError cError(ID_ERROR_RESULTS_2DVIEWCLASS_WRONG); //_T("Init2DChildViewDescriptor is wrong"));
			return;
		}

	}
	else {
		m_p2DChildViewDescriptor = new CViewDescriptor(pViewClass, nID);
	}
}

void CResults::Init3DChildViewDescriptor(CRuntimeClass *pViewClass, UNINT nID) {
	if (m_p3DChildViewDescriptor) {

		if (pViewClass != m_p3DChildViewDescriptor->GetViewClass() || nID != m_p3DChildViewDescriptor->GetResID()) {
			CDlgShowError cError(ID_ERROR_RESULTS_3DVIEWCLASS_WRONG); //_T("Init3DChildViewDescriptor is wrong"));
			return;
		}

	}
	else {
		m_p3DChildViewDescriptor = new CViewDescriptor(pViewClass, nID);
	}
}

CViewDescriptor* CResults::GetChildViewDescriptor() {
	return m_b3D ? m_p3DChildViewDescriptor : m_p2DChildViewDescriptor;
}

size_t CResults::LoadResults(void){
	std::vector<CString> files;

	//Файл задачи
	CString fpath;
	CFileDialog fileDlg(1, _T("emma_task"),0, 4|2, _T("task file (*.emma_task)|*.emma_task|"));

	if(fileDlg.DoModal() == IDOK){
		fpath = fileDlg.GetPathName();

		int nPosSlash = fpath.ReverseFind( '\\' ); 
		fpath.Truncate( nPosSlash+1 ); 
	}else{
		//CDlgShowError cErr(_T("CResults: failed"));
		return 0;
	}
	//LOGGER.Init(CString(_T("..\\..\\Logs\\Results.cpp_LoadResults.txt")));
	files = GetFolderListing(fpath, CString(_T("*.emma_task")));

	//Формируем список результатов
	size_t nLoaded=0;
	for (size_t i = 0; i < files.size(); i++){
		if (LoadStage(files[i])) {
			nLoaded++;
		}
	}

	return nLoaded;
}

CString CResults::MakeName(const CString &strFileName){
	CString strName = strFileName;
	for(int i=strName.GetLength()-1; i>=0; i--){
		//if(strName[i]==static_cast<WCHAR>('\\') || strName[i]==static_cast<WCHAR>('/')){
		if (strName[i] == _T('\\') || strName[i] == _T('/')) {
			strName =strName.Right(strName.GetLength()-1-i);
			break;
		}
	}
	for(int i=strName.GetLength()-1; i>=0; i--){
		if(strName[i] == _T('.')){
			strName = strName.Left(i);
			break;
		}
	}
	return strName;
}