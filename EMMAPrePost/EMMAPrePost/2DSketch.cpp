#include "StdAfx.h"
#include "2DSketch.h"


IOIMPL (C2DSketch, C2DSKETCH)

C2DSketch::C2DSketch(CString strName) : CEMMADoc(RUNTIME_CLASS(CChild2DView), IDR_EMMAPrePostTYPE, strName){
	RegisterMember(&m_Outline);
	SetIconName(IDS_ICON_SKETCH);
}

C2DSketch::C2DSketch(int icon, CString strName) : CEMMADoc(RUNTIME_CLASS(CChild2DView), IDR_EMMAPrePostTYPE, strName) {
	RegisterMember(&m_Outline);
	SetIconName(icon);
}

C2DSketch::C2DSketch(CRuntimeClass *pViewClass, size_t nID, CString strName) : CEMMADoc(pViewClass, nID, strName){
	RegisterMember(&m_Outline);
	SetIconName(IDS_ICON_SKETCH);
}

C2DSketch::C2DSketch(CRuntimeClass *pViewClass, size_t nID, int icon, CString strName) : CEMMADoc(pViewClass, nID, strName) {
	RegisterMember(&m_Outline);
	SetIconName(icon);
}

C2DSketch::~C2DSketch(void){
}

void C2DSketch::SetOutline(C2DOutline *pOutline){
	//m_Outline.ClearAll();
	VERIFY(m_Outline.Copy(pOutline));
	//m_OutlineAdapter.SetObject(m_Outline);

/////////////////////Roma
	m_CurvesActive.clear();
	m_NodesActive.clear();
///////////////////////////////

	InitialUpdateView();
	InitialUpdatePane();
}

void C2DSketch::SetContours() {
	C2DOutline *pOutline = GetOutline();

	if (!pOutline) {
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_CONTOURS); //_T("Contours not found."));	//Показываем окно ошибки
		return;
	}

	if (pOutline->FindContours(1.0))// длина ребра сетки
		GetView()->SendMessage(WMR_INITIALUPDATE); // обновляем представление
	else
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_CONTOURS); //_T("Contours not found."));	//Показываем окно ошибки
}


void C2DSketch::SetClearOutline(){
	
	GetOutline()->ClearAll();

	ClearActiveObjects();
	
	InitialUpdateView();
	InitialUpdatePane();
}


/*void C2DSketch::InitialClearView(){

	CRect2D rect;
	rect.left = 0;
	rect.right = 50;
	rect.bottom = 25;
	rect.top = 0;

	if(GetView()) GetView()->SendMessage(WMR_DRAFTRECT, (WPARAM)&rect, (LPARAM)1);
} //*/

//! Добавление свойства в таблицу справа
CMFCPropertyGridProperty* C2DSketch::AddProp(CStringW name, double eval, size_t pos){

	CMFCPropertyGridProperty *pProp = new CMFCPropertyGridProperty(name, eval, name);	//The group name, DATA, Description
	pProp->SetData(pos);	//устанавливает свойство в нужную позицию
	if (pos == 5) pProp->AllowEdit(0);	//точка внутри или снаружи
	return pProp;
}

//! Заполняет таблицу свойств
void C2DSketch::FillPropList(CMFCPropertyGridCtrl *pGrid){
	if(pGrid == 0) return;

	pGrid->RemoveAll();
	if ((m_CurvesActive.size() == 1) && (m_NodesActive.size() == 0)){
		
		C2DCurve * t_curve;
		t_curve = m_Outline.GetCurve(m_CurvesActive[0]);
			
		pGrid->AddProperty(AddProp(_T("X1"), t_curve->GetStartNode()->GetPoint().x,1));
		pGrid->AddProperty(AddProp(_T("Y1"), t_curve->GetStartNode()->GetPoint().y,2));
		pGrid->AddProperty(AddProp(_T("X2"), t_curve->GetEndNode()->GetPoint().x,3));
		pGrid->AddProperty(AddProp(_T("Y2"), t_curve->GetEndNode()->GetPoint().y,4));
		pGrid->AddProperty(AddProp(_T("L"), t_curve->CalcLength(),5));

		if (m_Outline.GetCurve(m_CurvesActive[0])->GetType() == C2DCIRCARC) {
			C2DCircleArc *pCurve;
			pCurve = dynamic_cast<C2DCircleArc*>(m_Outline.GetCurve(m_CurvesActive[0]));
			pGrid->AddProperty(AddProp(_T("R"),pCurve->GetRadius(),6));
			
		}

		return;
	}
	
	//Если кривая не выделена, значит нужен узел
	if (m_NodesActive.size() == 1){
	
		C2DNode * t_node;
		t_node = m_Outline.GetNode(m_NodesActive[0]);

		pGrid->AddProperty(AddProp(_T("X"), t_node->GetPoint().x,7));
		pGrid->AddProperty(AddProp(_T("Y"), t_node->GetPoint().y,8));

		return;
		}
		
	pGrid->Invalidate();
	
}

//! Обновляет таблицу свойств
void C2DSketch::UpdateProp(double *pEval, UNLONG *pId){

	if ((pEval == 0) || (pId == 0)) {
		return;
	}
	Math::C2DPoint tPoint;
	C2DCurve *pCur = nullptr;
	C2DCircleArc *pArc = nullptr;
	bool isCurveActive = false;

////////////// ПЕРЕПИСАТЬ ВСЕ В ЦИКЛАХ!!!

	if (m_CurvesActive.size() > 0){
		pCur = m_Outline.GetCurve(m_CurvesActive[0]);
		isCurveActive = true;
	}

		switch(*pId){
		case 1:
			if (isCurveActive){
				tPoint = pCur->GetStartNode()->GetPoint();
				tPoint.x = *pEval;
				pCur->GetStartNode()->SetPoint(tPoint);
			}
			break;

		case 2:
			if (isCurveActive) {
				tPoint = pCur->GetStartNode()->GetPoint();
				tPoint.y = *pEval;
				pCur->GetStartNode()->SetPoint(tPoint);
			}
			break;

		case 3:
			if (isCurveActive) {
				tPoint = pCur->GetEndNode()->GetPoint();
				tPoint.x = *pEval;
				pCur->GetEndNode()->SetPoint(tPoint);
			}
			break;

		case 4:
			if (isCurveActive) {
				tPoint = pCur->GetEndNode()->GetPoint();
				tPoint.y = *pEval;
				pCur->GetEndNode()->SetPoint(tPoint);
			}
			break;
	
		//5 - длина, не изменяем

		case 6:
			if (isCurveActive) {
				pArc = dynamic_cast<C2DCircleArc*>(pCur);
				pArc->SetRadius(*pEval);
			}
			break;

		//Точки в узлах
		case 7:
			tPoint = m_Outline.GetNode(m_NodesActive[0])->GetPoint();
			tPoint.x = *pEval;
			m_Outline.GetNode(m_NodesActive[0])->SetPoint(tPoint);
			break;

		case 8:
			tPoint = m_Outline.GetNode(m_NodesActive[0])->GetPoint();
			tPoint.y = *pEval;
			m_Outline.GetNode(m_NodesActive[0])->SetPoint(tPoint);
			break;
		}

	for(size_t i = 0; i < m_Outline.GetCurveCount() - 1; i++){
		if(m_Outline.GetCurve(i)->GetType() == C2DCIRCARC){
			C2DCircleArc *t_pArc = dynamic_cast<C2DCircleArc*>(m_Outline.GetCurve(i));
			t_pArc->CheckRadius(100);
		}
	}

	InitialUpdateView();
	//InitialUpdatePane();
}

void C2DSketch::UpdatePropList(CMFCPropertyGridCtrl * pGrid)
{
	//кол-во строк в таблице свойств (для проверки)
	int ncount = pGrid->GetPropertyCount();

	//конвертируем в MS дурацкий формат для таблицы
	C2DCurve * t_curve;
	t_curve = m_Outline.GetCurve(m_CurvesActive[0]);
	COleVariant varLength(t_curve->CalcLength());
	
	//устанавливаем значения для кол-ва эементов и узлов
	CMFCPropertyGridProperty* pProperty = pGrid->GetProperty(4);
	CMFCPropertyGridProperty* pProperty2 = pGrid->FindItemByData(5);
	pProperty->SetValue(&varLength);
}

bool C2DSketch::GetBoundingBox(CRect2D &rect){
	if (m_Outline.GetBoundingBox(rect)){
		return true;
	}
	return false;
}