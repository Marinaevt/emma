#include "StdAfx.h"
#include "I2DSketch.h"

I2DSketch::I2DSketch(void)
{
	m_NumCurve = -1;	//Номер кривой
	m_NumNode = -1;		//Номер узла

	m_CurvesActive.clear();
	m_NodesActive.clear();

	m_fAdd = false;
	m_fCtrl = false;
	m_dCurTime = 0;		//Текущее время (для PrePost только!)
	m_dTimeSpeed = 0.01;	//Скорость изменения Текущего времени m_dCurTime (для PrePost только!)
	m_fPlay = false;
}

I2DSketch::~I2DSketch(void)
{
}

void I2DSketch::SetSpeedOfTime(DBL dNewSpeedOfTime){
	m_dTimeSpeed = dNewSpeedOfTime;
}

/////Roma
void I2DSketch::DrawGL(GLParam &parameter){
	
	glColor3d(0,0,0);
	//glPointSize(1);
	GetOutline()->DrawGL(parameter);

	C2DCurve *pCur;
	// Отрисовка активных кривых.
	for (size_t i = 0; i < m_CurvesActive.size(); i++) {
		glColor3d(1,1,0);
		glLineWidth(3);
		pCur = GetOutline()->GetCurve(m_CurvesActive[i]);
		pCur->DrawGL(parameter);

		glColor3d(1,0,1);
		glPointSize(5);
		pCur->GetStartNode()->DrawGL(parameter);
		pCur->GetEndNode()->DrawGL(parameter);
	}
	
	// Отрисовка активных узлов.
	for (size_t i = 0; i < m_NodesActive.size(); i++) {
		glColor3d(0,1,1);
		glPointSize(7);
		GetOutline()->GetNode(m_NodesActive[i])->DrawGL(parameter);
	}
	
	if ((m_NumCurve != -1) && (m_NumNode == -1)){
		glColor3d(0,1,0);
		glLineWidth(3);
		pCur = GetOutline()->GetCurve(m_NumCurve);
		pCur->DrawGL(parameter);
		glColor3d(1,0,0);
		glPointSize(5);
		pCur->GetStartNode()->DrawGL(parameter);
		pCur->GetEndNode()->DrawGL(parameter);
	}

	if (m_NumNode != -1){
		glColor3d(0,0,1);
		glPointSize(7);
		GetOutline()->GetNode(m_NumNode)->DrawGL(parameter);
	}

	//Отрисовка всего остального.
	glColor3d(1,0,0);
	glPointSize(3);
	glLineWidth(2);
	m_Outline.DrawGL(parameter);
			
}

void I2DSketch::SetActiveObject(Math::C2DPoint *ppoint)
{	
	m_Outline.ClearAll();	//очищаем предыдущие вспомогательные рисунки
	Math::C2DPoint point_on_curve;
	
	DBL dL = 0.0; //расстояние от текущей точки до ближайшего узла
	DBL dEps = 0.0; //Расстояние до кривой или точки 
	
	//Нашли номер ближайшего узла и расстояние до него
	m_NumNode = GetOutline()->GetClosestNode(*ppoint, dL);

	//если у нас только узлы...
	if (GetOutline()->GetCurveCount() < 1) {
		return;
	}
	
	//если есть кривые, то смотрим, что ближе
	dEps = GetOutline()->GetSmallestCurveLenght() / 4.0;	 //(делим на 4, чтобы наводить на кривую, иначе только узлы)

	//узел ближе, чем кривая?
	if (dL < dEps) {
		//Мы уже запомнили его номер в m_NumNode, после клика узел станет активным
		return;
	}else {
		//Сбрасываем выделение
		m_NumNode = -1;
	}

	int nTypePoint = GetOutline()->GetClosestPoint(*ppoint, point_on_curve, m_NumCurve);

	//при ошибке не выделяем ничего
	if (nTypePoint == -1 ) {
		 m_NumCurve = -1;
		 return;
	}
	
	//сбрасываем выделение, если слишком далеко (больше длины наименьшей кривой)
	if ((*ppoint - point_on_curve).Norm() > dEps*4.0){
		m_NumCurve = -1;
		return;
	}

	// !!! У ВСЕХ узлов уже была возможность "выделиться" выше, доп. проверки (ниже) не нужны
	
}

//очищаем все активные элементы
void I2DSketch::ClearActiveObjects(){
	
	m_Outline.ClearAll();
	m_CurvesActive.clear();
	m_NodesActive.clear();
	m_NumCurve = -1;
	m_NumNode = -1;
}

//Создаём вспомогательную линию для подсказки
void I2DSketch::SetTempLine(Math::C2DPoint *ppoint) {

	m_Outline.ClearAll();

	if (GetOutline()->GetNodeCount() == 0) {
		return;
	}
	
	int nSecondNode = GetOutline()->GetNodeCount() - 1;	//если нет активных узлов - берём последний
	if (m_NodesActive.size() > 0) {

		if (m_NodesActive[m_NodesActive.size() - 1] != m_NumNode) {
			m_NumNode = m_NodesActive[m_NodesActive.size() - 1];
		}

		nSecondNode = m_NumNode;
	}

	DBL dL = 0.0, //расстояние от текущей точки до ближайшего узла
		dEps = GetOutline()->GetSmallestCurveLenght() / 10.0;	//Расстояние до кривой или точки (делим на 4, чтобы наводить на кривую, иначе только узлы)

	//Нашли номер ближайшего узла и расстояние до него
	int nClosestNode = GetOutline()->GetClosestNode(*ppoint, dL);
	if (dL < dEps) {
		//если расстояние до ближайшего узла "мало", то выбираем этот узел
		*ppoint = GetOutline()->GetNode(nClosestNode)->GetPoint();
	}

	int n1 = m_Outline.CreateNode(C2DNODE, *ppoint),
		n2 = m_Outline.CreateNode(C2DNODE, GetOutline()->GetNode(nSecondNode)->GetPoint());
	int c1 = m_Outline.CreateCurve(C2DCURVE, n1, n2);

}

//Создаём вспомогательную дугу для подсказки
void I2DSketch::SetTempArc(Math::C2DPoint *ppoint) {

	m_Outline.ClearAll();

	if (GetOutline()->GetNodeCount() == 0) {
		return;
	}

	int nSecondNode = GetOutline()->GetNodeCount() - 1;	//если нет активных узлов - берём последний
	if (m_NodesActive.size() > 0) {

		if (m_NodesActive[m_NodesActive.size() - 1] != m_NumNode) {
			m_NumNode = m_NodesActive[m_NodesActive.size() - 1];
		}

		nSecondNode = m_NumNode;
	}

	DBL dL = 0.0, //расстояние от текущей точки до ближайшего узла
		dEps = GetOutline()->GetSmallestCurveLenght() / 10.0;	//Расстояние до кривой или точки (делим на 4, чтобы наводить на кривую, иначе только узлы)

	//Нашли номер ближайшего узла и расстояние до него
	int nClosestNode = GetOutline()->GetClosestNode(*ppoint, dL);
	if (dL < dEps) {
		//если расстояние до ближайшего узла "мало", то выбираем этот узел
		*ppoint = GetOutline()->GetNode(nClosestNode)->GetPoint();
	}

	int n1 = m_Outline.CreateNode(C2DNODE, *ppoint),
		n2 = m_Outline.CreateNode(C2DNODE, GetOutline()->GetNode(nSecondNode)->GetPoint());
	DBL dRad = m_Outline.GetNode(n1)->GetPoint().Len(m_Outline.GetNode(n2)->GetPoint()) / 2.0 + EPS;	//диаметр пополам
	C2DCircleArc *pCircleArc = new C2DCircleArc(n2, n1, dRad);

	m_Outline.AddCurve(pCircleArc);
}

int I2DSketch::FindFacetNode(Math::C2DPoint *pPoint, int &nCurve1, int &nCurve2){
	double dDistance = 0.0; // Переменная для расстояния до ближайшего узла
	int nNode = GetOutline()->GetClosestNode(*pPoint, dDistance); //Ближайший узел

	if(nNode == -1 || dDistance > 1) return -1; // Если узел не найден или слишком далеко - фаску построить нельзя

	nCurve1 = nCurve2 = -1;
	
	for(size_t i = 0; i < GetOutline()->GetCurveCount(); i++){
		if((GetOutline()->GetCurve(i)->GetStart() == nNode) || (GetOutline()->GetCurve(i)->GetEnd() == nNode)){
			if(nCurve1 == -1){
				nCurve1 = i;
			}else if(nCurve2 == -1){
				nCurve2 = i;
			}else{//Если в узле сходится более трех кривых, прекращаем поиск
				nCurve2 = -1;
				break;
			}
		}
	}

	if(nCurve2 == -1) nNode = -1; // Если в узле сходится более или менее 2-х кривых - фаску построить нельзя
	return nNode;
}

void I2DSketch::CopyInM_OutlineLine(C2DFacet_2lines *out){
	
	int n1 = m_Outline.CreateNode(C2DNODE, out->GetCurve(0)->GetStartNode()->GetPoint()),
		n2 = m_Outline.CreateNode(C2DNODE, out->GetCurve(0)->GetEndNode()->GetPoint());
	int c1 = m_Outline.CreateCurve(C2DCURVE,n1,n2);

	n1 = m_Outline.CreateNode(C2DNODE,out->GetCurve(1)->GetStartNode()->GetPoint()),
	n2 = m_Outline.CreateNode(C2DNODE,out->GetCurve(1)->GetEndNode()->GetPoint());
	c1 = m_Outline.CreateCurve(C2DCURVE,n1,n2);

	n1 = m_Outline.CreateNode(C2DNODE,out->GetCurve(2)->GetStartNode()->GetPoint()),
	n2 = m_Outline.CreateNode(C2DNODE,out->GetCurve(2)->GetEndNode()->GetPoint());
	c1 = m_Outline.CreateCurve(C2DCIRCARC,n1,n2);

	C2DCircleArc *pCur1 = dynamic_cast<C2DCircleArc *>(m_Outline.GetCurve(c1));
	pCur1->SetRadius(dynamic_cast<C2DCircleArc *>(out->GetCurve(2))->GetRadius());
	pCur1->SetCW(dynamic_cast<C2DCircleArc *>(out->GetCurve(2))->GetCW());	
}

void I2DSketch::CopyInM_OutlineLineArc(C2DFacet_line_arc *out){
	int n1 = m_Outline.CreateNode(C2DNODE, out->GetCurve(0)->GetStartNode()->GetPoint()),
		n2 = m_Outline.CreateNode(C2DNODE, out->GetCurve(0)->GetEndNode()->GetPoint());
	int c1 = m_Outline.CreateCurve(C2DCURVE,n1,n2);
	
	//Первая дуга
	n1 = m_Outline.CreateNode(C2DNODE,out->GetCurve(1)->GetStartNode()->GetPoint()),
	n2 = m_Outline.CreateNode(C2DNODE,out->GetCurve(1)->GetEndNode()->GetPoint());
	int c2 = m_Outline.CreateCurve(C2DCIRCARC,n1,n2);

	C2DCircleArc *pCur1 = dynamic_cast<C2DCircleArc *>(m_Outline.GetCurve(c2));
	pCur1->SetRadius( dynamic_cast<C2DCircleArc *>(out->GetCurve(1))->GetRadius());
	pCur1->SetCW( dynamic_cast<C2DCircleArc *>(out->GetCurve(1))->GetCW());

	//Вторая дуга
	n1 = m_Outline.CreateNode(C2DNODE,out->GetCurve(2)->GetStartNode()->GetPoint()),
	n2 = m_Outline.CreateNode(C2DNODE,out->GetCurve(2)->GetEndNode()->GetPoint());
	int c3 = m_Outline.CreateCurve(C2DCIRCARC,n1,n2);

	C2DCircleArc *pCur2 = dynamic_cast<C2DCircleArc *>(m_Outline.GetCurve(c3));
	pCur2->SetRadius( dynamic_cast<C2DCircleArc *>(out->GetCurve(2))->GetRadius());
	pCur2->SetCW( dynamic_cast<C2DCircleArc *>(out->GetCurve(2))->GetCW());
}

void I2DSketch::CopyInM_OutlineArcs(C2DFacet_2arcs *out){
	int n1 = m_Outline.CreateNode(C2DNODE, out->GetCurve(0)->GetStartNode()->GetPoint()),
		n2 = m_Outline.CreateNode(C2DNODE, out->GetCurve(0)->GetEndNode()->GetPoint());
	int c1 = m_Outline.CreateCurve(C2DCIRCARC,n1,n2);
	C2DCircleArc *pCur1 = dynamic_cast<C2DCircleArc *>(m_Outline.GetCurve(c1));
	pCur1->SetRadius(dynamic_cast<C2DCircleArc *>(out->GetCurve(0))->GetRadius());
	pCur1->SetCW(dynamic_cast<C2DCircleArc *>(out->GetCurve(0))->GetCW());

	n1 = m_Outline.CreateNode(C2DNODE,out->GetCurve(1)->GetStartNode()->GetPoint()),
	n2 = m_Outline.CreateNode(C2DNODE,out->GetCurve(1)->GetEndNode()->GetPoint());
	c1 = m_Outline.CreateCurve(C2DCIRCARC,n1,n2);
	C2DCircleArc *pCur2 = dynamic_cast<C2DCircleArc *>(m_Outline.GetCurve(c1));
	pCur2->SetRadius(dynamic_cast<C2DCircleArc *>(out->GetCurve(1))->GetRadius());
	pCur2->SetCW(dynamic_cast<C2DCircleArc *>(out->GetCurve(1))->GetCW());

	n1 = m_Outline.CreateNode(C2DNODE,out->GetCurve(2)->GetStartNode()->GetPoint()),
	n2 = m_Outline.CreateNode(C2DNODE,out->GetCurve(2)->GetEndNode()->GetPoint());
	c1 = m_Outline.CreateCurve(C2DCIRCARC,n1,n2);
	C2DCircleArc *pCur3 = dynamic_cast<C2DCircleArc *>(m_Outline.GetCurve(c1));
	pCur3->SetRadius(dynamic_cast<C2DCircleArc *>(out->GetCurve(2))->GetRadius());
	pCur3->SetCW(dynamic_cast<C2DCircleArc *>(out->GetCurve(2))->GetCW());
}

void I2DSketch::CopyInOutlineLine(int i){
	
	C2DCurve *pCur = m_Outline.GetCurve(i);
	int n1 = GetOutline()->CreateNode(C2DNODE,pCur->GetStartNode()->GetPoint()),
		n2 = GetOutline()->CreateNode(C2DNODE,pCur->GetEndNode()->GetPoint());

	//C2DCurve *pLine = pCur->Clone();
	int c1 = GetOutline()->CreateCurve(C2DCURVE,n1,n2);
	//TODO: надо уточнять номер контура
	bool b = GetOutline()->AddCurveToContour(c1,0);
	//GetOutline()->GetCurve(c1)->SetNodes(n1,n2);
}

//приходится пересоздавать, т.к. номера кривых и узлов не совпадают
void I2DSketch::CopyInOutlineArc(int i){
	
	C2DCircleArc *pCur = dynamic_cast<C2DCircleArc *>(m_Outline.GetCurve(i));
	int n1 = GetOutline()->CreateNode(C2DNODE,pCur->GetStartNode()->GetPoint()),
		n2 = GetOutline()->CreateNode(C2DNODE,pCur->GetEndNode()->GetPoint());

	//C2DCircleArc *pArc = ((C2DCircleArc *) pCur)->Clone();
	//int c1 = GetOutline()->AddCurve(pArc);
	int c1 = GetOutline()->CreateCurve(C2DCIRCARC,n1,n2);
	//TODO: надо уточнять номер контура
	bool b = GetOutline()->AddCurveToContour(c1,0);
	
	C2DCircleArc *pCurOutline = dynamic_cast<C2DCircleArc *>(GetOutline()->GetCurve(c1));
	pCurOutline->SetRadius(pCur->GetRadius());
	pCurOutline->SetCW(pCur->GetCW());	
}

void I2DSketch::CopyInOutline(int nC1,int nC2, int nP, short t1, short t2){
	
	bool i1 = false,
		i2 = false,
		i3 = false;

	C2DCurve *pCur1 = GetOutline()->GetCurve(nC1),
			*pCur2 = GetOutline()->GetCurve(nC2);

	i1 = GetOutline()->DelCurve(nC1);
	i2 = GetOutline()->DelCurve(nC2-1);	//потому что при удалении предыдущего номера меняются
	i3 = GetOutline()->DelNode(nP);
	//*/
	if(t1 == C2DCIRCARC && t2 == C2DCURVE){
		CopyInOutlineArc(1);
		CopyInOutlineLine(0);
	}else if (t1 == C2DCURVE && t2 == C2DCIRCARC){
		CopyInOutlineLine(0);
		CopyInOutlineArc(1);
	}else if (t1 == C2DCURVE && t2 == C2DCURVE){
		//для отрезка-отрезка
		CopyInOutlineLine(0);
		CopyInOutlineLine(1);
	}else{
		CopyInOutlineArc(0);
		CopyInOutlineArc(1);
	}
	CopyInOutlineArc(2);
	
	//LOGGER.Init(CString(_T("..\\..\\Logs\\I2DSketch.cpp_CopyInM_Outline_end.txt")));
	//GetOutline()->WriteToLog();
}

// Вспомогательная дуговая фаска (сопряжение)
void I2DSketch::SetTempFacet(Math::C2DPoint *pPoint, bool set){
	
	m_Outline.ClearAll();  // Очищаем временный Чертёж
	int nCurve1, nCurve2; // Номера сопрягаемых кривых
	int nNode = FindFacetNode(pPoint, nCurve1, nCurve2); //Ищем узел для построения фаски
	if(nNode == -1) return;//Если узел не найден - выходим

	//Yura////////////////////////////////////////
	C2DCurve* pCur1 = GetOutline()->GetCurve(nCurve1);
	C2DCurve* pCur2 = GetOutline()->GetCurve(nCurve2);

	Math::C2DPoint s1=pCur1->GetStartNode()->GetPoint(),
		e1=pCur1->GetEndNode()->GetPoint(),
		s2=pCur2->GetStartNode()->GetPoint(),
		e2=pCur2->GetEndNode()->GetPoint(),
		ugl=GetOutline()->GetNode(nNode)->GetPoint();	// отмеченная точка
	short type1 = pCur1->GetType(),
		type2 = pCur2->GetType();
	bool retval = false;
	
	//Вычисление r, потом будем его принимать как параметр
	DBL r = 1.0; //(GetOutline()->GetCurve(nCurve1)->CalcLength()+GetOutline()->GetCurve(nCurve2)->CalcLength())/8;

	//Типы кривых - отрезки
	if ((type1 == C2DCURVE) && (type2 == C2DCURVE) ){
		C2DFacet_2lines Facet(ugl,r);

		Facet.SetTempFacetPointsOrder(s1,e1,s2,e2);
		Facet.CreateData();
		Facet.Facet_ugl(&(Facet.m_pVect1),&(Facet.m_pVect2),s1,e1,s2,e2);
		retval = Facet.Create();
		
		if (retval) {
			CopyInM_OutlineLine(&Facet);
		}

	//Типы кривых - дуга и отрезок
	}else if ( ((type1 == C2DCIRCARC) && (type2 == C2DCURVE)) || ((type1 == C2DCURVE) && (type2 == C2DCIRCARC)) ){
		
		C2DFacet_line_arc Facet(ugl,r);
		
		Facet.SetTempFacetPointsOrder(s1,e1,s2,e2);
		
		if (type1 == C2DCIRCARC){
			Facet.CreateData(type1,type2,(dynamic_cast<C2DCircleArc *>(pCur1))->GetCW(),(dynamic_cast<C2DCircleArc *>(pCur1))->GetRadius());
		}else{
			Facet.CreateData(type1,type2,(dynamic_cast<C2DCircleArc *>(pCur2))->GetCW(),(dynamic_cast<C2DCircleArc *>(pCur2))->GetRadius());
		}
		
		//Правка дуг при заполнении Facet
		Facet.CorrectArcs(type1,type2,nCurve1,nCurve2, GetOutline());
		
		//Нахождение направляющих векторов
		Facet.Facet_ugl(&(Facet.m_pVect1),&(Facet.m_pVect2),s1,e1,s2,e2);

		//Создание фаски
		retval = Facet.Create();

		if (retval) {
			CopyInM_OutlineLineArc(&Facet);		//Копируем в Outline как подсказку
		}
	
	//Типы кривых - дуги
	}else if ( (type1 == C2DCIRCARC) && (type2 == C2DCIRCARC) ){
		C2DFacet_2arcs Facet(ugl,r);
		
		Facet.SetTempFacetPointsOrder(s1,e1,s2,e2);
		Facet.CreateData( dynamic_cast<C2DCircleArc *>(pCur1)->GetCW(), dynamic_cast<C2DCircleArc *>(pCur1)->GetRadius(),	
			dynamic_cast<C2DCircleArc *>(pCur2)->GetCW(), dynamic_cast<C2DCircleArc *>(pCur2)->GetRadius());

		//{{Правка дуг при заполнении Facet
		Facet.CorrectArcs(nCurve1,nCurve2, GetOutline());

		//Нахождение направляющих векторов
		Facet.Facet_ugl(&(Facet.m_pVect1),&(Facet.m_pVect2),s1,e1,s2,e2);
		
		retval = Facet.Create();

		if (retval) {
			CopyInM_OutlineArcs(&Facet);	//Копируем в Outline как подсказку
		}
	}

	//Если нажали на кнопку мыши, то фиксируем результат 
	if (retval==true && set==true){
		CopyInOutline(nCurve1,nCurve2, nNode, type1, type2);
	}

	//Yura////////////////////////////////////////
}

void I2DSketch::OnClick(int *cur){
	
	if ((GetOutline()->GetType() == C2DMOTION) && (*cur != 0)) {
		return;
	}

	if (!m_fCtrl){

		m_CurvesActive.clear();
		m_NodesActive.clear();

		if (m_NumCurve != -1) {
			m_CurvesActive.push_back(m_NumCurve);
		}
		if (m_NumNode != -1) {
			m_NodesActive.push_back(m_NumNode);
		}
	
	}else{
		
		if (m_NumCurve != -1){
			for (size_t i = 0; i < m_CurvesActive.size(); i++)
				if (m_NumCurve == m_CurvesActive[i]){
					m_CurvesActive.erase(m_CurvesActive.begin() + i);
					return;
				}
			m_CurvesActive.push_back(m_NumCurve);
		}
					
		if (m_NumNode != -1){
			for (size_t i = 0; i < m_NodesActive.size(); i++)
				if (m_NumNode == m_NodesActive[i]){
					m_NodesActive.erase(m_NodesActive.begin() + i);
					return;
				}
			m_NodesActive.push_back(m_NumNode);
		}
				

	}

}

void I2DSketch::AddNode(Math::C2DPoint *ppoint){
	
	//m_CurvesActive.clear();
	//m_NodesActive.clear();

	//C2DNode *t_node = new C2DNode(*ppoint);

	int nNode = GetOutline()->CreateNode(C2DNODE, *ppoint);
	m_NumNode = nNode;	//переключили выделение на новый узел
	
	if (GetOutline()->GetType() == C2DMOTION){
		CIO2DMotion* t_motion;
		t_motion = dynamic_cast<CIO2DMotion*>(GetOutline());
		t_motion->InitParams(GetOutline()->GetNodeCount());
	}

}


int I2DSketch::AddFirstNode(Math::C2DPoint *ppoint) {

	//C2DNode *t_node = new C2DNode(*ppoint);

	int nNode = GetOutline()->CreateNode(C2DNODE,*ppoint);

	if (GetOutline()->GetType() == C2DMOTION) {
		CIO2DMotion* t_motion;
		t_motion = dynamic_cast<CIO2DMotion*>(GetOutline());
		t_motion->InitParams(GetOutline()->GetNodeCount());
	}
	return nNode;
}
//*/

void I2DSketch::AddCurve(Math::C2DPoint *ppoint){
	
	//! Добавляем кривые в последний созданный контур
	C2DOutline* pOut = GetOutline();
	if (pOut->GetContourCount() < 1) {
		size_t Res = pOut->CreateContour(C2DCONTOUR);
	}
	C2DContour* t_contour = pOut->GetContour(pOut->GetContourCount() - 1);

	int nFirstNode = 0;
	//если забыли поставить первую точку
	if (pOut->GetNodeCount() == 0 ){

		Math::C2DPoint* pFirstNode = new Math::C2DPoint(ppoint->m_X, ppoint->m_Y);
		nFirstNode = AddFirstNode(pFirstNode);
		m_NumNode = nFirstNode;
		m_NodesActive.push_back(m_NumNode);
		return;
	}else{
		DBL dL = 0.0, //расстояние от текущей точки до ближайшего узла
			dEps = pOut->GetSmallestCurveLenght() / 10.0;	//Расстояние до кривой или точки (делим на 4, чтобы наводить на кривую, иначе только узлы)

		//Нашли номер ближайшего узла и расстояние до него
		int nClosestNode = pOut->GetClosestNode(*ppoint, dL);
		if (dL < dEps) {
			*ppoint = pOut->GetNode(nClosestNode)->GetPoint();
		}
		nFirstNode = AddFirstNode(ppoint);
	}

	int nSecondNode = 0;
	if (m_NodesActive.size() > 0) {
		nSecondNode = m_NodesActive[m_NodesActive.size() - 1];
	}
	else {
		//на данном этапе мы точно создали вторую точку
		nSecondNode = pOut->GetNodeCount() - 2;
	}

	//если выбрали ту же точку, с которой начали, то ничего не делаем
	if (nSecondNode == nFirstNode) {
		return;
	}

	int nAddNum = pOut->CreateCurve(C2DCURVE, nSecondNode, nFirstNode);
	m_CurvesActive.push_back(nAddNum);
	m_NumNode = nFirstNode;
	m_NodesActive.push_back(m_NumNode);

	bool bAdd = t_contour->AddCurve(pOut->GetCurveCount()-1);
	bool bAdd2 = false;
	if (!bAdd) {
		//Если не удалось записать кривую в контур, скорее всего из-за разветвлений (из узла несколько кривых).
		//Поэтому создаём новый контур и пишем в него.
		//t_contour = new C2DContour();
		size_t Res = pOut->CreateContour(C2DCONTOUR); //AddContour(t_contour);
		bAdd2 = pOut->GetContour(Res)->AddCurve(pOut->GetCurveCount() - 1);
	}
	
	//! Если мы рисуем траекторию
	if (pOut->GetType() == C2DMOTION) {

		CIO2DMotion * pOutline = dynamic_cast<CIO2DMotion*>(pOut);
		pOutline->CalcTime();
	}
}

void I2DSketch::AddCircleArc(Math::C2DPoint *ppoint){
	
	//! Добавляем кривые в последний созданный контур
	C2DOutline* pOut = GetOutline();
	if (pOut->GetContourCount() < 1) {
		size_t Res = pOut->CreateContour(C2DCONTOUR);
	}
	C2DContour* t_contour = pOut->GetContour(pOut->GetContourCount() - 1);

	int nFirstNode = 0;
	//если забыли поставить первую точку
	if (pOut->GetNodeCount() == 0 ){
		Math::C2DPoint* pFirstNode = new Math::C2DPoint(ppoint->m_X, ppoint->m_Y);
		nFirstNode = AddFirstNode(pFirstNode);
		m_NumNode = nFirstNode;
		m_NodesActive.push_back(m_NumNode);
		return;
	}else{
		DBL dL = 0.0, //расстояние от текущей точки до ближайшего узла
			dEps = pOut->GetSmallestCurveLenght() / 10.0;	//Расстояние до кривой или точки (делим на 4, чтобы наводить на кривую, иначе только узлы)

															//Нашли номер ближайшего узла и расстояние до него
		int nClosestNode = pOut->GetClosestNode(*ppoint, dL);
		if (dL < dEps) {
			*ppoint = pOut->GetNode(nClosestNode)->GetPoint();
		}
		nFirstNode = AddFirstNode(ppoint);
	}

	int nSecondNode = 0;
	if (m_NodesActive.size() > 0) {
		nSecondNode = m_NodesActive[m_NodesActive.size() - 1];
	}
	else {
		//на данном этапе мы точно создали вторую точку
		nSecondNode = pOut->GetNodeCount() - 2;
	}

	//если выбрали ту же точку, с которой начали, то ничего не делаем
	if (nSecondNode == nFirstNode) {
		return;
	}

	DBL dRad = pOut->GetNode(nFirstNode)->GetPoint().Len(pOut->GetNode(nSecondNode)->GetPoint()) / 2.0 + EPS;	//диаметр пополам
	C2DCircleArc* t_arc = new C2DCircleArc(nSecondNode, nFirstNode, dRad);
	
	int nAddNum = pOut->AddCurve(t_arc);
	m_CurvesActive.push_back(nAddNum);
	m_NumNode = nFirstNode;
	m_NodesActive.push_back(m_NumNode);

	bool bAddArc = t_contour->AddCurve(GetOutline()->GetCurveCount() - 1);
	bool bAddArc2 = false;
	if (!bAddArc) {
		//Если не удалось записать кривую в контур, скорее всего из-за разветвлений (из узла несколько кривых).
		//Поэтому создаём новый контур и пишем в него.
		t_contour = new C2DContour();
		size_t Res = pOut->AddContour(t_contour);
		bAddArc2 = t_contour->AddCurve(pOut->GetCurveCount() - 1);
	}


	//! Если мы рисуем траекторию
	if (GetOutline()->GetType() == C2DMOTION) { 

		CIO2DMotion * pOutline = dynamic_cast<CIO2DMotion*>(GetOutline());
		pOutline->CalcTime();
	}

}

//! Добавление элемента Траектории на чертёж
void I2DSketch::AddItem(Math::C2DPoint *ppoint, int cur){

	//При событии OnLButtonUp
	//Зависит от m_nCurType - курсор;
	if (cur == 100) AddNode(ppoint);
	if (cur == 101) AddCurve(ppoint);
	if (cur == 102) AddCircleArc(ppoint);
}

void I2DSketch::ChangeTime(GLParam &parameter){
	
	//(анимация перемещения инструмента по траектории PrePost)
	m_dCurTime = m_dCurTime + m_dTimeSpeed;	//Меняем время на величину dt 
	DrawGL(parameter);

}
