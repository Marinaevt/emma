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

	if (pOutline->FindContours(1.0)) {// длина ребра сетки
		GetView()->SendMessage(WMR_INITIALUPDATE); // обновляем представление
	}
	else {
		CDlgShowError Diag_err(ID_ERROR_2DMESHERPANE_NO_CONTOURS); //_T("Contours not found."));	//Показываем окно ошибки
	}
}


void C2DSketch::SetClearOutline(){
	
	GetOutline()->ClearAll();

	ClearActiveObjects();
	
	InitialUpdateView();
	InitialUpdatePane();
}


//! Добавление свойства в таблицу справа
CMFCPropertyGridProperty* C2DSketch::AddProp(CStringW name, double eval, size_t pos){
	
	CMFCPropertyGridProperty *pProp = new CMFCPropertyGridProperty(name, eval, name);	//The group name, DATA, Description
	pProp->SetData(pos);	//устанавливает свойство в нужную позицию
	if (pos == 6) {
		pProp->Enable(0);	//выключает поле свойства, не редактируется
	}
	return pProp;
}

//! Заполняет таблицу свойств
void C2DSketch::FillPropList(CMFCPropertyGridCtrl *pGrid){
	
	if (pGrid == 0) {
		return;
	}

	pGrid->RemoveAll();
	if ((m_CurvesActive.size() == 1) && (m_NodesActive.size() == 0)){
		
		C2DCurve * t_curve;
		t_curve = m_Outline.GetCurve(m_CurvesActive[0]);
			
		pGrid->AddProperty(AddProp(_T("X1"), t_curve->GetStartNode()->GetPoint().x,1));
		pGrid->AddProperty(AddProp(_T("Y1"), t_curve->GetStartNode()->GetPoint().y,2));
		pGrid->AddProperty(AddProp(_T("X2"), t_curve->GetEndNode()->GetPoint().x,3));
		pGrid->AddProperty(AddProp(_T("Y2"), t_curve->GetEndNode()->GetPoint().y,4));

		if (m_Outline.GetCurve(m_CurvesActive[0])->GetType() == C2DCIRCARC) {
			C2DCircleArc *pCurve;
			pCurve = dynamic_cast<C2DCircleArc*>(m_Outline.GetCurve(m_CurvesActive[0]));
			pGrid->AddProperty(AddProp(_T("R"),pCurve->GetRadius(),5));
			
		}

		pGrid->AddProperty(AddProp(_T("L"), t_curve->CalcLength(), 6));

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

//! Обновляем все поля таблицы свойств
void C2DSketch::UpdatePropList(CMFCPropertyGridCtrl * pGrid)
{
	if (pGrid == nullptr) {
		CDlgShowError err(_T("CMFCPropertyGridCtrl is null"));
		return;
	}
	
	//кол-во строк в таблице свойств (для проверки)
	int nCount = pGrid->GetPropertyCount();
	if (nCount == 0) {
		CDlgShowError err(_T("CMFCPropertyGridCtrl is empty"));
		return;
	}

	Math::C2DPoint tPoint;
	C2DCurve *pCur = nullptr;
	C2DCircleArc *pArc = nullptr;
	bool isCurveActive = false;

	if (m_CurvesActive.size() > 0) {
		pCur = m_Outline.GetCurve(m_CurvesActive[0]);	//берём только первую из выделенных
		isCurveActive = true;
	}

	if (m_NodesActive.size() > 0) {	//захватили узел, но выделение кривой не снялось (близко)
		isCurveActive = false;
	}

	//обновляем все поля

	//если выбрали кривую, то определяем узлы, длину (и радиус для дуг)
	if (isCurveActive) {
		CMFCPropertyGridProperty* pProperty1 = pGrid->FindItemByData(1);	//X1
		if (pProperty1) {
			tPoint = pCur->GetStartNode()->GetPoint();
			tPoint.x = pProperty1->GetValue().dblVal;
			pCur->GetStartNode()->SetPoint(tPoint);
		}

		CMFCPropertyGridProperty* pProperty2 = pGrid->FindItemByData(2);	//Y1
		if (pProperty2) {
			tPoint = pCur->GetStartNode()->GetPoint();
			tPoint.y = pProperty2->GetValue().dblVal;
			pCur->GetStartNode()->SetPoint(tPoint);
		}

		CMFCPropertyGridProperty* pProperty3 = pGrid->FindItemByData(3);	//X2
		if (pProperty3) {
			tPoint = pCur->GetEndNode()->GetPoint();
			tPoint.x = pProperty3->GetValue().dblVal;
			pCur->GetEndNode()->SetPoint(tPoint);
		}

		CMFCPropertyGridProperty* pProperty4 = pGrid->FindItemByData(4);	//Y2
		if (pProperty4) {
			tPoint = pCur->GetEndNode()->GetPoint();
			tPoint.y = pProperty4->GetValue().dblVal;
			pCur->GetEndNode()->SetPoint(tPoint);
		}

		//5 - радиус
		CMFCPropertyGridProperty* pProperty5 = pGrid->FindItemByData(5);	//Radius
		if (pProperty5) {
			pArc = dynamic_cast<C2DCircleArc*>(pCur);
			pArc->SetRadius(pProperty5->GetValue().dblVal);
		}

		//6 - длина, не изменяем, но обновляем
		CMFCPropertyGridProperty* pProperty6 = pGrid->FindItemByData(6);	//Length
		if (pProperty6) {
			COleVariant varLength(pCur->CalcLength());
			pProperty6->SetValue(&varLength);
		}
	}
	else {
		//Координаты выделенного узла
		CMFCPropertyGridProperty* pProperty7 = pGrid->FindItemByData(7);	//X
		if (pProperty7) {
			tPoint = m_Outline.GetNode(m_NodesActive[0])->GetPoint();
			tPoint.x = pProperty7->GetValue().dblVal;
			m_Outline.GetNode(m_NodesActive[0])->SetPoint(tPoint);
		}

		CMFCPropertyGridProperty* pProperty8 = pGrid->FindItemByData(8);	//Y
		if (pProperty8) {
			tPoint = m_Outline.GetNode(m_NodesActive[0])->GetPoint();
			tPoint.y = pProperty8->GetValue().dblVal;
			m_Outline.GetNode(m_NodesActive[0])->SetPoint(tPoint);
		}
	}

	//Проверяем радиусы у дуг
	size_t nCurves = m_Outline.GetCurveCount();
	for (size_t i = 0; i < nCurves; i++) {
		if (m_Outline.GetCurve(i)->GetType() == C2DCIRCARC) {
			C2DCircleArc *t_pArc = dynamic_cast<C2DCircleArc*>(m_Outline.GetCurve(i));
			bool bCheck = t_pArc->CheckRadius(100);
			//TODO: ??? in UpdatePropList() ???
		}
	}

	InitialUpdateView();
}

bool C2DSketch::GetBoundingBox(CRect2D &rect){
	if (m_Outline.GetBoundingBox(rect)){
		return true;
	}
	return false;
}