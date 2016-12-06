#include "StdAfx.h"
#include "2DTrajectory.h"
#include "2DTool.h"	//исключение

IOIMPL (C2DTrajectory, C2DTRAJECTORY)

C2DTrajectory::C2DTrajectory(void) : CEMMADoc(RUNTIME_CLASS(CChild2DView), IDR_EMMAPrePostTYPE){
	RegisterMember(&m_Motion);

	SetClearOutline();
	Math::C2DPoint *p = new Math::C2DPoint(0,0);
	m_Motion.CreateNode(C2DNODE, *p);
	
//	m_dCurTime = 0;
//	m_nState = 10000;
	m_fPlay = false;
	SetIconName(IDS_ICON_TRAJECTORY);
}

C2DTrajectory::~C2DTrajectory(void){
}

//! доступ к полученной траектории для внешних классов
CIO2DMotion* C2DTrajectory::GetMotion(){
	return &m_Motion;
}

void C2DTrajectory::SetOutline(C2DOutline *pOutline){
	
	SetClearOutline();
	VERIFY(m_Motion.Copy(pOutline)); //Так скорее всего нельзя делать

	m_CurvesActive.clear();
	m_NodesActive.clear();

	InitialUpdateView();
	InitialUpdatePane();
}

void C2DTrajectory::SetClearOutline(){

	GetOutline()->ClearAll();
	
	ClearActiveObjects();
	InitialClearView();

	GetOutline()->CreateContour(C2DCONTOUR);

	C2DPieceLinearFunction *t_const = new C2DPieceLinearFunction();
	t_const->m_args().push_back(0.0);
	t_const->m_values().push_back(1.0);

	C2DConstantFunction *t_const1 = new C2DConstantFunction();
	t_const1->m_value() = 1.0;

	CIO2DMotion *t_motion;
	t_motion = dynamic_cast<CIO2DMotion*>(GetOutline());
	
	t_motion->m_vels().push_back(t_const);
	t_motion->m_vels().push_back(t_const1);
	InitialUpdateView();
	InitialUpdatePane();
}

void C2DTrajectory::InitialClearView(){
	
	CRect2D rect = m_Outline.GetBoundingBox();
	if(GetView()) GetView()->SendMessage(WMR_DRAFTRECT, (WPARAM)&rect, (LPARAM)1);
}

//! Добавление свойства в таблицу справа
CMFCPropertyGridProperty* C2DTrajectory::AddProp(CStringW name, double eval, size_t pos){

	CMFCPropertyGridProperty *pProp = new CMFCPropertyGridProperty(name, eval, name);
	pProp->SetData(pos);
	if (pos == 5) pProp->AllowEdit(0);	//длина дуги
	if (pos == 10) pProp->AllowEdit(0);	//время в
	if (pos == 11) pProp->AllowEdit(0);	//время из
	return pProp;
}

//! Заполнение параметров (правая панель)
void C2DTrajectory::FillPropList(CMFCPropertyGridCtrl *pGrid){
	
	if(pGrid == nullptr) return;

	//! Если активная кривая
	if ((m_CurvesActive.size() == 1) && (m_NodesActive.size() == 0)){
		pGrid->RemoveAll();
		C2DCurve * t_curve;
		t_curve = m_Motion.GetCurve(m_CurvesActive[0]);
		

		pGrid->AddProperty(AddProp(_T("X1"), t_curve->GetStartNode()->GetPoint().x,1));
		pGrid->AddProperty(AddProp(_T("Y1"), t_curve->GetStartNode()->GetPoint().y,2));
		pGrid->AddProperty(AddProp(_T("X2"), t_curve->GetEndNode()->GetPoint().x,3));
		pGrid->AddProperty(AddProp(_T("Y2"), t_curve->GetEndNode()->GetPoint().y,4));
		pGrid->AddProperty(AddProp(_T("L"), t_curve->CalcLength(),5));

		if (m_Motion.GetCurve(m_CurvesActive[0])->GetType() == C2DCIRCARC) {
			C2DCircleArc *pCurve;
			pCurve = dynamic_cast<C2DCircleArc*>(m_Motion.GetCurve(m_CurvesActive[0]));
			pGrid->AddProperty(AddProp(_T("R"),pCurve->GetRadius(),6));
			
		}

		//скорость изменения времени (только PrePost и только для анимации передвижения)
		pGrid->AddProperty(AddProp(_T("SpeedOfTime"), m_dTimeSpeed,15));

		pGrid->Invalidate();
		return;
	}

	//! Если активная точка
	if ((m_NodesActive.size() == 1)){
		pGrid->RemoveAll();
		C2DNode * t_node;
		t_node = m_Motion.GetNode(m_NodesActive[0]);
		C2DPieceLinearFunction* vel = dynamic_cast<C2DPieceLinearFunction *>(m_Motion.m_vels()[0]);
		C2DFunction* angVel = dynamic_cast<C2DFunction *>(m_Motion.m_vels()[1]);
		C2DMotionParams mot_param;
		
		pGrid->AddProperty(AddProp(_T("X"), t_node->GetPoint().x,7));
		pGrid->AddProperty(AddProp(_T("Y"), t_node->GetPoint().y,8));
		
		if (m_Motion.m_params().size() != 0){
			mot_param = m_Motion.m_params()[m_NodesActive[0]];
			pGrid->AddProperty(AddProp(_T("Время (в)"),mot_param.GetEntrTime(),10));
			pGrid->AddProperty(AddProp(_T("Время (из)"),mot_param.GetExitTime(),11));
			pGrid->AddProperty(AddProp(_T("Пауза"),mot_param.GetPause(),12));
		}

		//линейная скорость
		pGrid->AddProperty(AddProp(_T("V"), vel->m_values().at(m_NodesActive[0]),13));
		//угловая скорость
		pGrid->AddProperty(AddProp(_T("W"), angVel->GetValue(m_NodesActive[0]),14));
		//скорость изменения времени (только PrePost и только для анимации передвижения)
		pGrid->AddProperty(AddProp(_T("SpeedOfTime"), m_dTimeSpeed,15));
		
		pGrid->Invalidate();
		return;

	}
		
	pGrid->RemoveAll();
	pGrid->Invalidate();

}

//! Обновление параметров (правая панель)
void C2DTrajectory::UpdateProp(double *pEval, UNLONG *pId){
	
	if( (pEval == nullptr) || (pId == nullptr)) return;
	Math::C2DPoint tPoint;
	C2DCircleArc *pArc=nullptr;
	C2DCurve *pCur=nullptr;
	bool isCurveActive = false;
////////////// ПЕРЕПИСАТЬ ВСЕ В ЦИКЛАХ!!!
		
	if (m_CurvesActive.size() > 0){
		pCur = m_Motion.GetCurve(m_CurvesActive[0]);
		isCurveActive = true;
	}

	int p;
	C2DPieceLinearFunction* vel = dynamic_cast<C2DPieceLinearFunction *>(m_Motion.m_vels()[0]);	// линейная скорость
	C2DPieceLinearFunction* angVel = dynamic_cast<C2DPieceLinearFunction *>(m_Motion.m_vels()[1]);	//угловая скорость

		switch(*pId){
		
		//! Для кривых
		case 1:
			if (isCurveActive) {
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
		
		//case 5 возвращает длину сразу

		case 6:
			if (isCurveActive) {
				pArc = dynamic_cast<C2DCircleArc*>(pCur);
				pArc->SetRadius(*pEval);
			}
			break;

		//! Для точек
		case 7:
			tPoint = m_Motion.GetNode(m_NodesActive[0])->GetPoint();
			tPoint.x = *pEval;
			m_Motion.GetNode(m_NodesActive[0])->SetPoint(tPoint);
			break;

		case 8:
			tPoint = m_Motion.GetNode(m_NodesActive[0])->GetPoint();
			tPoint.y = *pEval;
			m_Motion.GetNode(m_NodesActive[0])->SetPoint(tPoint);
			break;
		
		//TODO: 9 почему-то пропустили

		case 10:	//время в
			m_Motion.m_params()[m_NodesActive[0]].SetEntrTime(*pEval);
			break;

		case 11:	//время из
			m_Motion.m_params()[m_NodesActive[0]].SetExitTime(*pEval);
			break;

		case 12:	//пауза
			p = m_NodesActive[0];
			m_Motion.m_params()[p].SetPause(*pEval);
			break;

		case 13:	//линейная скорость [0]
			
			vel->m_values().at(m_NodesActive[0]) = (*pEval);
			
			//vel->m_args().SetAt(m_NodesActive[0],vel->m_args().GetAt(m_NodesActive[0]-1)+(m_Motion.GetNode(m_NodesActive[0]-1)->GetPoint() - m_Motion.GetNode(m_NodesActive[0])->GetPoint()).Norm()/(vel->m_values()[m_NodesActive[0]-1]+(vel->m_values()[m_NodesActive[0]]-vel->m_values()[m_NodesActive[0]-1])/2));
			//vel->m_args().SetAt(m_NodesActive[0],vel->m_args().GetAt(m_NodesActive[0])+2*(m_Motion.GetNode(m_NodesActive[0])->GetPoint()-m_Motion.GetNode(m_NodesActive[0]-1)->GetPoint()).Norm()/(vel->m_values()[m_NodesActive[0]-1]+vel->m_values()[m_NodesActive[0]]));
			//m_Motion.m_vels()[0].

			//m_Motion.CalcTime();
			//m_Motion.m_params()[m_NodesActive[0]].SetPause(*pEval);
			break;

		case 14:	// угловая скорость [1]
			
			//TODO: !!!!!! при изменении значения - валится
			//angVel->m_values().at(m_NodesActive[0]) = (*pEval);
			//m_Motion.m_params()[m_NodesActive[0]].SetPause(*pEval);
			break;

		case 15:	// новая скорость изменения времени (dt) 
			
			SetSpeedOfTime(*pEval);
			break;
		
		}

	m_Motion.CalcTime(); //перерасчёт времени после изменений

	for(size_t i = 0; i < m_Motion.GetCurveCount() - 1; i++){
		if(m_Motion.GetCurve(i)->GetType() == C2DCIRCARC){
			C2DCircleArc *t_pArc = dynamic_cast<C2DCircleArc*>(m_Motion.GetCurve(i));
			t_pArc->CheckRadius(100);
		}
	}

	InitialUpdateView();
	InitialUpdatePane();
}

void C2DTrajectory::SetContours() {
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

//! Отрисовка Траектории
void C2DTrajectory::DrawGL(GLParam &parameter){
	
	I2DSketch::DrawGL(parameter);
	if (m_fPlay)	DrawTool(parameter);		
}

//! Установка Инструмента в начало Траектории и отрисовка
void C2DTrajectory::DrawTool(GLParam &parameter){
	
	C2DPosition posit;
	if (!m_Motion.GetPos(m_dCurTime,posit)){
		return;
	};

	DBL dX, dY;
	m_Motion.GetVel(0.5, dX, dY);

	C2DTool *tool = dynamic_cast<C2DTool*> (m_pParentDoc);
	C2DTool tmp;
	if (tool){
		if(tmp.C2DTool::Copy(tool)){
			tmp.ChangeNodes(&posit.m_pos);
			tmp.DrawGL(parameter);
		}
	}
}

bool C2DTrajectory::GetBoundingBox(CRect2D &rect){
	if (m_Motion.GetBoundingBox(rect)){
		return true;
	} 
	return false;
}
