#include "StdAfx.h"
#include "2DBuilder.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#define CRUTCH 0
IOIMPL (C2DBuilder, C2DBUILDER)

//потом надо будет создать собственное представление для этого документа
C2DBuilder::C2DBuilder(void): CEMMADoc(RUNTIME_CLASS(CCalcStage2DView), IDR_EMMAPrePostTYPE) {
	SetIconName(IDS_ICON_BUILDER);
}


C2DBuilder::~C2DBuilder(void){
}

bool C2DBuilder::GetBoundingBox(CRect2D &rect){

	//rect = CRect2D(0,15,0,1);
	//return true;
	CGen2DDoc *pDoc = dynamic_cast<CGen2DDoc *>(GetRoot());	//получаем ссылку на корневой документ (2D)
	CGen2DOperation *pDocTools = dynamic_cast<CGen2DOperation *>(GetParent());	//получаем ссылку на родительский документ (Операции)

	if (!pDoc && !pDocTools){
		return false;
	}

	if (pDoc){
		pDoc->GetMesher()->GetBoundingBox(rect);
	}

	if (pDocTools){
		ptrdiff_t nToolNum = pDocTools->GetToolNum();
		for (ptrdiff_t j=0; j<nToolNum; j++){

			C2DTool *pTool = pDocTools->GetTool(j);
			if(pTool){
				C2DTrajectory *pTraj = pTool->GetTrajectory();
				rect.AddRect(pTraj->GetOutline()->GetBoundingBox());
				rect.AddRect(pTool->GetOutline()->GetBoundingBox());
			}
		}
	}
	return true;
}

//! Отрисовка КЭ сетки
void C2DBuilder::DrawMesh(CGen2DDoc *pDoc, GLParam &parameter){
	//return;
	if(pDoc){
		const CMesh *pMesh = pDoc->GetMesh();

		if(pMesh->m_nElNum != 0){
			pMesh->DrawGl(parameter);
		}
	}
	
}

//! Отрисовка Инструментов
void C2DBuilder::DrawTools(CGen2DOperation *pDocTools, GLParam &parameter){
	if (pDocTools){
		ptrdiff_t nToolNum = pDocTools->GetToolNum();
		for (ptrdiff_t j=0; j<nToolNum; j++){

			C2DTool *pTool = pDocTools->GetTool(j);
			if(pTool){
				pTool->DrawGL(parameter);
			}
		}
	}
}

//! Отрисовка Траекторий с Инструментами
void C2DBuilder::DrawTrajectories(CGen2DOperation *pDocTools, GLParam &parameter){
	
	if (pDocTools){
		ptrdiff_t nToolNum = pDocTools->GetToolNum();
		for (ptrdiff_t j=0; j<nToolNum; j++){

			C2DTool *pTool = pDocTools->GetTool(j);
			if(pTool){
				C2DTrajectory *pTraj = pTool->GetTrajectory();
				pTraj->m_dCurTime = 0;
				pTraj->DrawTool(parameter);
				pTraj->DrawGL(parameter);
			}
		}
	}
}

void C2DBuilder::DrawGL(GLParam &parameter){
	CGen2DDoc *pDoc = dynamic_cast<CGen2DDoc *>(GetRoot());	//получаем ссылку на корневой документ (2D)
	CGen2DOperation *pDocTools = dynamic_cast<CGen2DOperation *>(GetParent());	//получаем ссылку на родительский документ (Операции)

	DrawMesh(pDoc, parameter);	//Отрисовка КЭ сетки

	//DrawTools(pDocTools);			//Отрисовка Инструментов (сами по себе)

	DrawTrajectories(pDocTools, parameter);	//Отрисовка Траекторий (с Инструментами)
	
}

//! Сборка всего проекта для расчётов в ядре
bool C2DBuilder::Build(){
	
	//получаем ссылку на корневой документ ("2D"/"Прокатка"/"Волочение"...)
	CGen2DDoc *pDoc = dynamic_cast<CGen2DDoc *>(GetRoot());
	//получаем ссылку на родительский документ Сборщика ("Операции")
	CGen2DOperation *pDocTools = dynamic_cast<CGen2DOperation *>(GetParent());

	if(!pDoc) return false;
	if(!pDocTools) return false;
	//LOGGER.Init(CString("..\\..\\Logs\\C2DBuilder.cpp_Build.txt"));

	//! Указатель на КЭ Сетку
	const CMesh *pMesh = pDoc->GetMesh();
	if(pMesh == nullptr){ 
		CDlgShowError Diag_err(ID_ERROR_2DBUILDER_NO_MESH);	//Показываем окно ошибки
		return false;
	}
	else if (pMesh->m_nElNum == 0) {
		CDlgShowError Diag_err(ID_ERROR_2DBUILDER_NO_MESH);	//Показываем окно ошибки
		return false;
	}
	//pMesh->WriteToLog();

	//! Указатель на список материалов
	CMaterialsList *pMaterialsList = pDoc->GetMaterialList();
	if (pMaterialsList == nullptr) {
		CDlgShowError Diag_err(ID_ERROR_2DBUILDER_MATERIAL_LOAD_FAIL); //_T("Файл материала был неправильно загружен");
		return false;
	}
	CMaterialDoc *pMaterialDoc = pMaterialsList->GetFirstMaterialDoc();	//TODO: загружаем первый материал из списка
	CString strMaterialFile; 
	if (pMaterialDoc != nullptr) {
		strMaterialFile = pMaterialDoc->GetMaterialFileName();
	}
	else {
		CDlgShowError Diag_err(ID_ERROR_2DBUILDER_MATERIAL_SELECT_FAIL); //(_T("Материал не загружен"));
		return false;
	}

	//! Узлы, материал, граничные условия
	C2DPlaneFEM* p2DFem = new C2DPlaneFEM;
	
	bool nInitMesh = p2DFem->InitMesh(pMesh);
	bool nInitMaterial = p2DFem->InitMaterial(strMaterialFile);
	//p2DFem->WriteBCToLog();

	//! Формирование файла задачи
	C2DTask Task;
	Task.SetName(pDocTools->GetName()); //В качестве имени - название операции
	Task.AddObject(p2DFem);

	//! Сохраняем инструменты и траектории
	bool bIsTools = false;
	bIsTools = PrepareTools(&Task, pDocTools);
	if (!bIsTools) {
		CDlgShowError Diag_err(ID_ERROR_2DBUILDER_INSTRUM_NOT_LOADED); //_T("Инструменты не загружены");
		return false;
	}


	// ВРЕМЕННО!!! {{
	if(pDoc->GetType() == CGEN2DDOC){
		C2DSymmetry *pSymX = new C2DSymmetry(C2DSymmetry::modeX);
		C2DSymmetry *pSymY = new C2DSymmetry(C2DSymmetry::modeY);
		Task.AddObject(pSymX);
		Task.AddObject(pSymY);
	}

	//Волочение
	if (pDoc->GetType() == C2D_DRAGGING_DOC) {
		//добавляем симметрию по Y, волочение сверху вниз
		C2DSymmetry *pSymY = new C2DSymmetry(C2DSymmetry::modeY);
		Task.AddObject(pSymY);
	}

	//Формовка
	if(pDoc->GetType() == C2D_FORMING_DOC){
		C2DSymmetry *pSymY = new C2DSymmetry(C2DSymmetry::modeY);
		Task.AddObject(pSymY);

		C2DPressure *pPressure = new C2DPressure();

		C2DMeshInterface *m = p2DFem->GetMesh();
		size_t nSize = m->GetBorderSize();
		pPressure->initBorderNodes(nSize);
		//TODO: устанавливать по-человечески
		pPressure->SetValue(3.);
		//в формовке давление идёт только на нижние горизонтальные точки (обход в одну сторону):
		for(size_t i=0; i<nSize; i++){
			Math::C2DPoint p = m->GetBorderNode(i);
			Math::C2DPoint p1 = m->GetBorderNode(i+1 < nSize ? i+1 : 0);
			/*//Для шарика {{
			Math::C2DPoint p0 = m->GetBorderNode(i-1 >= 0 ? i-1 : m->GetBorderSize()-1);
			if(p1.x-p.x > 0 && p0.y>0.0001){// }}
			/*/
			if(p1.x-p.x > 0){//*/
				pPressure->setBorderNode(i, 1);
			}
		}
		Task.AddObject(pPressure);
	}
	if (pDoc->GetType() == C2D_BALLOON_DOC) {
		C2DSymmetry *pSymY = new C2DSymmetry(C2DSymmetry::modeY);
		Task.AddObject(pSymY);
		C2DSymmetry *pSymX = new C2DSymmetry(C2DSymmetry::modeX);
		Task.AddObject(pSymX);
		C2DPressure *pPressure = new C2DPressure();
		C2DMeshInterface *m = p2DFem->GetMesh();
		size_t nSize = m->GetBorderSize();
		pPressure->initBorderNodes(nSize);
		pPressure->SetValue(3.);
		for (size_t i = 0; i<nSize; i++) {
			Math::C2DPoint p = m->GetBorderNode(i);
			Math::C2DPoint p1 = m->GetBorderNode(i + 1 < nSize ? i + 1 : 0);
			if (p1.x - p.x > 0) {//*/
#if CRUTCH
				if (p.y > 0)
#endif
					pPressure->setBorderNode(i, 1);
			}
		}
		Task.AddObject(pPressure);
	}
	// }}

	//! Сохраняем задачу в файл
	bool bIsSave = false;
	bIsSave = SaveTask(&Task);

	return bIsSave;
}

bool C2DBuilder::PrepareTools(C2DTask* Task, CGen2DOperation* pDocTools){
	
	ptrdiff_t nToolNum = pDocTools->GetToolNum();
	if (nToolNum < 1) return false;

	for (ptrdiff_t j=0; j<nToolNum; j++){

		C2DTool* pTool = pDocTools->GetTool(j);
		if(pTool){
			C2DRigid* p2DRigid = new C2DRigid;
			p2DRigid->SetShape(pTool->GetOutline());	//устанавливаем инструмент
			p2DRigid->SetStickPoint(pTool->m_StickPoint);	//устанавливаем точку привязки
			p2DRigid->SetFrictionCoeff(pTool->m_FrictionCoeff); //устанавливаем коэффициент трения
			p2DRigid->SetMotion(pTool->GetTrajectory()->GetMotion());	//устанавливаем траекторию

			Task->AddObject(p2DRigid);
		}
	}
	return true;
}

bool C2DBuilder::SaveTask(C2DTask* Task){
	CFileDialog fileDlg(0, _T("emma_task"),0, 4|2, _T("task file (*.emma_task)|*.emma_task|"));
	if(fileDlg.DoModal() == IDOK){
		CString fname = fileDlg.GetPathName();
		CString fpath = fname;
		
		int nPosSlash = fpath.ReverseFind( '\\' ); 
		fname = fname.Right(fname.GetLength() - nPosSlash - 1);
		fpath.Truncate( nPosSlash+1 ); 

		Task->SetFileName(fname);		//устанавливаем имя файла
		Task->SetFilePath(fpath);		//устанавливаем путь до директории

		CStorage file;
		//LOGGER.Init(_T("../../Logs/YD_Save.txt"));	//логгер для отлова багов сохранения проекта
		if(file.Open(fpath + fname, CStorage::modeCreate | CStorage::modeWrite | CStorage::typeBinary | CStorage::shareDenyWrite)){
			return Task->Save(file);
		}
	}else{
		return false;
	}
	return false;
}