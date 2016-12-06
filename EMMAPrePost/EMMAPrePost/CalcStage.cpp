#include "StdAfx.h"
#include "CalcStage.h"


IOIMPL (CCalcStage, CCALCSTAGE)

CCalcStage::CCalcStage(void): m_pTask(nullptr){
	SetIconName(IDS_ICON_RESULTDOC); 
	m_path() = _T("");
	m_IsLoaded() = false;
	RegisterMembers();
}

CCalcStage::CCalcStage(CRuntimeClass *pViewClass, size_t nID)
			:CEMMADoc(pViewClass, nID),
			 m_pTask(nullptr)
{
	SetIconName(IDS_ICON_RESULTDOC); 
	m_path() = _T("");
	m_IsLoaded() = false;
	RegisterMembers();
}

CCalcStage::~CCalcStage(void){
	DeleteTask();
}

void CCalcStage::RegisterMembers(){
	RegisterMember(&m_path);
	RegisterMember(&m_IsLoaded);
	RegisterMember(m_pTask);
}

CString CCalcStage::MakeName(const CString &strFileName){
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
			strName =strName.Left(i);
			break;
		}
	}
	return strName;
}

bool CCalcStage::Set3DFlag(bool bFlag) {
	CResults *pRes = dynamic_cast<CResults *>(m_pParentDoc);
	if (pRes) {
		pRes->m_b3D = bFlag;
		return true;
	}
	return false;
}

//пока не используется (07.02.2016)
bool CCalcStage::Get3DFlag() {
	if (dynamic_cast<CResults *>(m_pParentDoc)) {
		return (dynamic_cast<CResults *>(m_pParentDoc))->m_b3D;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////
bool CCalcStage::LoadPathToTask(const CString &strFileName){
	CStorage file;

	if(file.Open(strFileName, CStorage::modeRead | CStorage::typeBinary)){
		
		m_path() = strFileName;
		file.Close();
		return true;
	}
	return false;
}

//! После выбора результата (click) загружаем файл в m_pTask
bool CCalcStage::LoadTask(){
	if (m_path() == _T("")) return false;
	
	CStorage file;
	if(file.Open(m_path(), CStorage::modeRead | CStorage::typeBinary)){
		DeleteTask();
		m_pTask = dynamic_cast<ITask *>(IO::LoadObject(file));
		file.Close();
		if (m_pTask){
			m_pTask->SetName(MakeName(m_path()));
			m_IsLoaded() = true;
			InitialUpdateView();
			return true;
		}else{
			DeleteTask();
			return false;
		}
	}
	
	return false;
}

void CCalcStage::DeleteTask(){
	if(m_pTask != nullptr){
		delete m_pTask;
		m_pTask = nullptr;
	}
	m_IsLoaded() = false;
}

void CCalcStage::DrawGL(GLParam &parameter)
{
	if (!m_IsLoaded) LoadTask();
	if (m_pTask){
		m_pTask->DrawGL(parameter);
	}
}

void CCalcStage::DrawGL3D(GLParam &parameter)
{
	if (!m_IsLoaded) LoadTask();
	if (m_pTask) {
		m_pTask->DrawGL3D(parameter);
	}
}

CViewDescriptor *CCalcStage::GetViewDescriptor(){
	if(!m_pParentDoc)return nullptr;
	CViewDescriptor *pViewDescr = (dynamic_cast<CResults *>(m_pParentDoc))->GetChildViewDescriptor();
	return pViewDescr;
}
