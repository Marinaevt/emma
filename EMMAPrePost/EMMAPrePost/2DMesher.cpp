#include "StdAfx.h"
#include "2DMesher.h"

IOIMPL (C2DMesher, C2DMESHER)

C2DMesher::C2DMesher(void) : CEMMADoc(RUNTIME_CLASS(CChild2DView), IDR_EMMAPrePostTYPE){
	
	SetIconName(IDS_ICON_MESH);

	RegisterMembers();

	/*m_NumCurve = -1;
	m_NumCurveActive = -1;
	m_TypePoint = -1;
	m_TypePointActive = -1; //*/

}

C2DMesher::C2DMesher(CRuntimeClass *pViewClass, size_t nID) : CEMMADoc(pViewClass, nID){
	
	SetIconName(IDS_ICON_MESH);

	RegisterMembers();

	/*m_NumCurve = -1;
	m_NumCurveActive = -1;
	m_TypePoint = -1;
	m_TypePointActive = -1; //*/
}

void C2DMesher::RegisterMembers(){
	
	RegisterMember(&m_Outline);
	RegisterMember(&m_Mesh);	//TODO: переписать подклассы для IO
}

void C2DMesher::SetClearOutline()
{	
	GetOutline()->ClearAll();
	ClearActiveObjects();
}

C2DMesher::~C2DMesher(void){
}
//! доступ к полученной сетке для внешних классов
CMesh* C2DMesher::GetMesh()
{
	CMesh* pMesh = &m_Mesh;
	return pMesh;
}

void C2DMesher::SetMeshType(UNINT id) {
	if(GetMesh())
		GetMesh()->SetMeshType(id);
}

UNINT C2DMesher::GetMeshType() {
	//if (GetMesh()) {
		return GetMesh()->GetMeshType();
	//}
	//по-умолчанию
	return ID_MESH_TYPE_FRONT;
}

void C2DMesher::SetContours() {
	C2DOutline *pOutline = GetOutline();

	if (!pOutline) {
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_CONTOURS); //_T("Contours not found."));	//Показываем окно ошибки
		return;
	}

	if(pOutline->FindContours(GetMesh()->m_dRibLen))
		GetView()->SendMessage(WMR_INITIALUPDATE); // обновляем представление
	else
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_CONTOURS); //_T("Contours not found."));	//Показываем окно ошибки
}


//Установка чертежа (после нажатия ОК в диалоге)
void C2DMesher::SetOutline(C2DOutline *pOutline)
{
	//LOGGER.Init(CString("..\\..\\Logs\\2DMesher_SetOutline.txt"));
	if (pOutline == nullptr) return;

	SetClearOutline();
	
	VERIFY(m_Outline.Copy(pOutline));
	
	SetContours();
	
	//bool bCondition = m_Outline.PrerareToMeshing();
	//if (!bCondition) return;	//TODO: показать диалог ошибки о некорректном чертеже\контуре

	//InitialUpdateView();	//Установка вида по центру
	DBL dLen = GetInitLen();	//Получение начальной длины разбиения

	GetMesh()->DeleteMesh();
	bool bGen = GetMesh()->GenerateMesh(dLen, &m_Outline);
	GetMesh()->ConstructMap();

	InitialUpdateView();	//Обновление вида после генерации и оптимизации КЭ сетки
	InitialUpdatePane();	//Обновление панели после генерации
	
}

DBL C2DMesher::GetInitLen(){
	
	DBL dLenMin = GetOutline()->GetCurve(0)->CalcLength(),
		dLenMax = dLenMin;	//в начале устанавливаем равными
	DBL dTmp = 0.0;

	//Находим меньшую из сторон
	for (size_t i=1; i<GetOutline()->GetCurveCount(); i++){
		
		dTmp = GetOutline()->GetCurve(i)->CalcLength();
		
		if (dTmp < dLenMin){
			dLenMin = dTmp;
		}

		if (dTmp > dLenMax){
			dLenMax = dTmp;
		}
	}

	DBL nRet = 0.0;
	if (dLenMin*10 >= dLenMax){	//отличия почти нет
		nRet = dLenMin;	//берём минимальную сторону
	}else{	//отличие на порядок и более
		nRet = dLenMax / dLenMin;
	}
	
	return nRet;
}

void C2DMesher::DrawGL(GLParam &parameter)
{
	glColor3d(0,0,0);
	glPointSize(1);
	GetOutline()->DrawGL(parameter); //рисуем чертеж
	GetMesh()->DrawGl(parameter);  //рисуем сетку
}

//! Добавление свойства в таблицу справа
CMFCPropertyGridProperty* C2DMesher::AddProp(CStringW name, double eval, size_t pos){

	CMFCPropertyGridProperty *pProp = new CMFCPropertyGridProperty(name, eval, name);	//The group name, DATA, Description
	pProp->SetData(pos);	//устанавливает свойство в нужную позицию
	if (pos == 3) pProp->AllowEdit(0);	//кол-во эл-ов
	if (pos == 4) pProp->AllowEdit(0);	//кол-во узлов
	return pProp;
}

//! Заполняет собой таблицу свойств
void C2DMesher::FillPropList(CMFCPropertyGridCtrl *pGrid)
{
	if (pGrid == 0) return;	

	pGrid->RemoveAll();
	CMesh* pMesh = GetMesh();
	if (pMesh->GetElNum() > 0) {
		pGrid->AddProperty(AddProp(_T("RibLength"), GetMesh()->m_dRibLen, 1));
		pGrid->AddProperty(AddProp(_T("PointRadius"), GetMesh()->m_dMeshRad, 2));
		pGrid->AddProperty(AddProp(_T("Elements"), GetMesh()->m_nElNum, 3));
		pGrid->AddProperty(AddProp(_T("Points"), GetMesh()->m_nPointsNum, 4));
		return;
	}
	
	pGrid->Invalidate();	//ненужно обновлять

}


void C2DMesher::UpdateProp(double *pEval, UNLONG *pId)
{
	if( (pEval == 0) || (pId == 0)) return;

	switch(*pId){
		
		//case 1: Кол-во элементов (не редактируется)
		
		//case 2: Кол-во точек (не редактируется)

		case 1:
			
			GetMesh()->m_dRibLen = *pEval;
			if (GetMesh()->m_dRibLen <= 0.0) {
				GetMesh()->m_dRibLen = 1.0;
			}
			if (GetOutline()->GetCurveCount() > 1 && GetOutline()->GetContourCount() > 0){
				GetMesh()->DeleteMesh();
				bool bGen = GetMesh()->GenerateMesh(m_Mesh.m_dRibLen, &m_Outline);
				GetMesh()->ConstructMap();
			}
			break;

		case 2:
			GetMesh()->m_dMeshRad = *pEval;
			if (GetMesh()->m_dMeshRad < 0.0) {
				GetMesh()->m_dMeshRad = 1.0;
			}
			if(GetMesh()->m_IsSetPoint){
				GetMesh()->DeleteMesh();
				bool bGen = GetMesh()->GenerateMesh(m_Mesh.m_dRibLen, &m_Outline);
				GetMesh()->ConstructMap();
			}
			break;
	}

	InitialUpdateView();
	//InitialUpdatePane();
}

void C2DMesher::UpdatePropList(CMFCPropertyGridCtrl * pGrid)
{
	//кол-во строк в таблице свойств (для проверки)
	int ncount = pGrid->GetPropertyCount();
	
	//конвертируем в MS дурацкий формат для таблицы
	COleVariant varElements((ULONGLONG) GetMesh()->GetElNum());
	COleVariant varNodes((ULONGLONG)GetMesh()->GetNodesNum());

	//устанавливаем значения для кол-ва элементов и узлов
	/*
	CMFCPropertyGridProperty* pProperty = pGrid->GetProperty(0);
	pProperty->SetValue(&varElements);
	CMFCPropertyGridProperty* pProperty2 = pGrid->GetProperty(1);
	pProperty2->SetValue(&varNodes);
	//*/
	CMFCPropertyGridProperty* pProperty1 = pGrid->FindItemByData(3);
	pProperty1->AllowEdit(0);
	CMFCPropertyGridProperty* pProperty2 = pGrid->FindItemByData(4);
	pProperty2->AllowEdit(0);

	pProperty1->SetValue(&varElements);
	pProperty2->SetValue(&varNodes);

}

void C2DMesher::AddMeshPoint(Math::C2DPoint ppoint){
	
	CMesh* pMesh = GetMesh();
	//Устанавливаем точку в любом случае
	pMesh->m_MeshPoint = ppoint;
	pMesh->m_IsSetPoint = true;

	if (pMesh->m_nElNum > 0){

		pMesh->DeleteMesh();
		pMesh->GenerateMesh(m_Mesh.m_dRibLen, &m_Outline);
		pMesh->ConstructMap();
	}
	InitialUpdateView();
}

bool C2DMesher::GetBoundingBox(CRect2D &rect){
	if (m_Outline.GetBoundingBox(rect)){
		return true;
	}
	return false;
}