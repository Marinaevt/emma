#include "StdAfx.h"
#include "2DRigid.h"


IOIMPL (C2DRigid, T2DRIGID)

C2DRigid::C2DRigid()
{
	RegisterMember(&m_shape);
	RegisterMember(&m_motion);
	RegisterMember(&m_stickpoint);
	RegisterMember(&m_FrictionCoeff);
	m_oldpos = Math::C2DPoint::Zero;
	//m_IsStuck = false;
	//m_IsFirstStuck =false;
	ResetTime();
}

void C2DRigid::ResetTime()
{
	m_time = 0;
}

bool C2DRigid::Init()
{
	//LOGGER.Init(CString(_T("..\\..\\Logs\\C2DRigid.cpp_Init.txt")));
	Math::C2DPoint dpos;
	
	// Находим начальную точку траектории, это будет первоначальная "старая" позиция
	m_oldpos = m_motion.GetNode(0)->GetPoint();
	
	// Смещение Инструмента
	dpos = m_oldpos - m_stickpoint();

	// Меняем положение Инструмента
	for (size_t i = 0; i < m_shape.GetNodeCount(); i++){
		m_shape.GetNode(i)->SetPoint(m_shape.GetNode(i)->GetPoint() += dpos);
	}
	
	//m_shape.WriteToLog();
	return true;
}

void C2DRigid::Move(double dt)
{
	
	Math::C2DPoint pos, tau, dpos;
	
	DBL l = dynamic_cast<C2DFunction*>(m_motion.m_vels()[0])->GetIntegral(0, m_time + dt);
	m_motion.GetPoint(0, l, pos, tau);

	//Изменение позиции Инструмента
	dpos = pos - m_oldpos;

	for (size_t i = 0; i < m_shape.GetNodeCount(); i++){
		tau = m_shape.GetNode(i)->GetPoint();
		m_shape.GetNode(i)->SetPoint(tau + dpos);
	}

	m_oldpos = pos;
	m_time += dt;
	
}

void C2DRigid::Calc(double dt)
{
}

//! Точка внутри инструмента или нет (точность len)
/*
bool C2DRigid::IsInside(const Math::C2DPoint& point, Math::C2DPoint& closep, DBL len)
{
	
	//точность пока не учитывается, ищется ближайший
	
	short inside = m_shape.IsInside(point);
	if (inside == 0 || inside == 1) {
		
		DBL dist;
		int nNode = m_shape.GetClosestNode(point,dist);	//получаем ближайший узел
		if(nNode == -1) return false;

		closep = m_shape.GetNode(nNode)->GetPoint();
		Math::C2DPoint minim;
		
		//находим все кривые с этим узлом
		for(size_t i=0; i<m_shape.GetCurveCount(); i++){

			C2DCurve *pCur = m_shape.GetCurve(i);
			if (pCur->GetStart() == nNode || pCur->GetEnd() == nNode){
				
				//находим ближайшую точку на кривой и сравниваем с предыдущей
				int p = pCur->GetClosestPoint(point,minim);
				if (p == -1) return false;
				if (dist > point.Len(minim)) { 
					dist = point.Len(minim);
					closep = minim;
				}

				//
				// если точность достаточна
				//if(dist <= len){

				//}
			}
		}
		
		return true;
	}
	return false;
}

*/

bool C2DRigid::GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc)
{
	/*// не используем этот метод (раньше - да)
	if (!pMesh) return false;

	// Если с прошлого раза мы были присоединены, то это уже не первая стыковка
	// (значит трение или отцепились)
	if (m_IsStuck){
		m_IsFirstStuck = false;
	}

	m_IsStuck = false;	//заново проходим по граничным узлам (вдруг отцепились)

	//Устанавливаем ГУ только точкам внутри\на границе инуструмента
	if (m_shape.IsInside(pMesh->GetBorderNode(nBoundaryNode)) > 0) {
			
		m_IsStuck = true; //Пристыковались

		C2DPosition tmp;
		m_motion.GetPos(m_time, tmp);

		// <-- для осесимметричной задачи (х == 0, y == 0 - заделка)
		if (m_IsStuck && m_IsFirstStuck 
			|| fabs(m_FrictionCoeff - 1.0) < EPS
			|| (fabs(pMesh->GetBorderNode(nBoundaryNode).x) < EPS 
			//	&&	fabs(pMesh->GetBorderNode(nBoundaryNode).y) < EPS
			)){ 
			//Заделка
			bc.setKinematic(tmp.m_vel.m_x, tmp.m_vel.m_y);

		}else if (m_IsStuck && !m_IsFirstStuck){
	
			DBL dSigmaNormal_x = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::sigma_x);
			DBL dSigmaNormal_y = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::sigma_y);
			DBL dSigmaInt = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::int_s);
			DBL dSquare = pMesh->GetCircleSquare(nBoundaryNode);

			double dRes = Friction(dSigmaNormal_x,dSigmaNormal_y,dSigmaInt,dSquare);
			//dRes /= (1.000000001 - m_FrictionCoeff);

			//Случай для одноосного сжатия
			bc.setSymX(tmp.m_vel.m_y, dRes);
		}
	}

	//Если контакта не было (после перебора), то сбрасываем значения и отбой
	// отцепились
	if (!m_IsStuck){
		m_IsFirstStuck = true;	//следующая стыковка будет "первой"
		return false;
	}

	return true;
	//*/
	return false;	//тут остановки быть не должно
}


//используем этот метод
bool C2DRigid::GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom> *bc) {

	if (!pMesh || pMesh == nullptr) {
		CDlgShowError cError(ID_ERROR_2DRIGID_MESHINTERFACE_NULL); //_T("Указатель на C2DMeshInterface null"));
		return false;
	}

	size_t number_of_element = 0;
	size_t nSize = pMesh->m_bordernodes.GetSize();// количество граничных узлов

	//Устанавливаем ГУ только точкам внутри\на границе инуструмента
	//Цикл по граничным узлам
	for (size_t nBoundaryNode = 0; nBoundaryNode < nSize; nBoundaryNode++)
	{
		//C2DCurve* closestCurve = m_shape.GetClosestCurve((pMesh->GetBorderNode(nBoundaryNode)));// ближайшая кривая к узлу
		size_t nBoundaryNode2 = nBoundaryNode > 0 ? nBoundaryNode - 1 : nSize - 1;	//соседний узел
		short inside = m_shape.IsInside(pMesh->GetBorderNode(nBoundaryNode));
		
		if (inside != -1)
		{

			C2DPosition tmp;
			C2DNode node;
			m_motion.GetPos(m_time, tmp);
			number_of_element = nBoundaryNode;
			DBL dist;
			//DBL dist_1;
			int nNode = m_shape.GetClosestNode(pMesh->GetBorderNode(nBoundaryNode), dist);	//получаем ближайший узел
			if (nNode == -1) return false;
			DBL	closep = m_shape.GetNode(nNode)->GetPoint();
			Math::C2DPoint minim, clstnd;
			//находим все кривые с этим узлом
			for (size_t i = 0; i < m_shape.GetCurveCount(); i++) 
			{
				C2DCurve *pCur = m_shape.GetCurve(i);
				if (pCur->GetStart() == nNode || pCur->GetEnd() == nNode)
				{

					//находим ближайшую точку на кривой и сравниваем с предыдущей
					int p = pCur->GetClosestPoint(pMesh->GetBorderNode(nBoundaryNode), minim);
					//int p_1 = pCur->GetClosestPoint(pMesh->GetBorderNode(nBoundaryNode2), minim);
					if (p == -1) return false;
					
					//	DBL m = pMesh->GetBorderNode(nBoundaryNode).Len(minim);
					//n1 = m_shape.GetNode((pCur->GetStart() == nNode ? pCur->GetEnd() : pCur->GetStart()))->GetPoint();

					if (dist > pMesh->GetBorderNode(nBoundaryNode).Len(minim)) 
					{
						dist = pMesh->GetBorderNode(nBoundaryNode).Len(minim);// получаем расстояние от точки Заготовки до Инструмента
						//	dist_1 = pMesh->GetBorderNode(nBoundaryNode2).Len(minim);

						clstnd = minim;
						closep = dist;
					}
				}
			}


			DBL testangle = m_shape.GetNode(nNode)->GetPoint().x - clstnd.x ? atan((m_shape.GetNode(nNode)->GetPoint().y - clstnd.y) / (m_shape.GetNode(nNode)->GetPoint().x - clstnd.x)) : 1.5708;




			
			DBL node_2_x = pMesh->GetBorderNode(nBoundaryNode).x;
			DBL node_2_y = pMesh->GetBorderNode(nBoundaryNode).y;
			DBL node_1_x = pMesh->GetBorderNode(nBoundaryNode2).x;
			DBL node_1_y = pMesh->GetBorderNode(nBoundaryNode2).y;

		

			
			//pMesh->GetBorderNode(nBoundaryNode).Splitting();

			DBL phy;
			//if ((nBoundaryNode || nBoundaryNode2) && nBoundaryNode > nBoundaryNode2)
			//{   // находим угол между точкой и горизонтом
				DBL cos_phy = abs(((node_2_x - node_1_x) /
					(sqrt((node_2_x - node_1_x)*(node_2_x - node_1_x) + (node_2_y - node_1_y)*(node_2_y - node_1_y)))));
				DBL sin_phy = abs(((node_2_y - node_1_y) /
					(sqrt((node_2_x - node_1_x)*(node_2_x - node_1_x) + (node_2_y - node_1_y)*(node_2_y - node_1_y)))));
				DBL angle_1 = asin(sin_phy);
				DBL angle_2 = acos(cos_phy);

				

				//tmp.m_pos.m_x = node_2_x - dist;
				//	tmp.m_pos.m_y = node_2_y - dist*sin_phy;
				//tmp.m_pos.m_x = pMesh->m_nodes[nBoundaryNode].m_x - dist;
				//tmp.m_pos.m_y = pMesh->m_nodes[nBoundaryNode].m_y + dist*sin_phy;

				//Узел заготовки движется с инструментом / прилипание
				if ((fabs(m_FrictionCoeff - 1.0) < EPS || fabs(node_2_x) < EPS)) //fabs(pMesh->GetBorderNode(nBoundaryNode).x) < EPS)
				{
					ALOGI("DD", AllToString(m_time) + CString(_T(" | ")) // было g_Step 
						+ AllToString(pMesh->GetBorderNode(nBoundaryNode).x) + CString(_T(" | "))
						+ AllToString(pMesh->GetBorderNode(nBoundaryNode).y) + CString(_T(" | "))
						+ AllToString(tmp.m_pos.m_x) + CString(_T(" | "))
						+ AllToString(tmp.m_pos.m_y) + CString(_T(" | "))
						+ AllToString(tmp.m_vel.m_x) + CString(_T(" | "))
						+ AllToString(tmp.m_vel.m_y) + CString(_T(" | "))
						+ AllToString(tmp.m_vel) + CString(_T(" | "))
						+ AllToString("Kinematic") + CString(_T(" | "))
						+ AllToString(number_of_element) + CString(_T(" | "))
						+ AllToString(dist) + CString(_T(" | "))
						+ AllToString(closep)
						);

					bc->at(nBoundaryNode).setKinematic(tmp.m_vel.m_x, tmp.m_vel.m_y, angle_1);
					
					//Узел заготовки скользит (в данном случае с учетом силы трения)
				}
				else// if (0 < m_shape.IsInside(pMesh->GetBorderNode(nBoundaryNode))
					///&& (0 < m_FrictionCoeff < 0.9)	
				{

					DBL dSigmaNormal_x = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::sigma_x);
					DBL dSigmaNormal_y = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::sigma_y);
					DBL dSigmaInt = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::int_s);
					DBL dEpsilonInt = pMesh->GetNField(pMesh->m_bordernodes[nBoundaryNode], eFields::int_d);
					DBL dSquare = pMesh->GetCircleSquare(nBoundaryNode);

					double dRes = Friction(dSigmaNormal_x, dSigmaNormal_y, dSigmaInt, dSquare);
					DBL tmp1 = m_time;//g_Step.dt;
					DBL phy_degree = (phy)*57.3;

					//DBL	velX = tmp.m_vel.m_x//*cos_phy
					//	     + tmp.m_vel.m_y//*sin_phy
					 //tmp.m_pos.m_y/tmp1
					//+ tmp.m_pos.m_x / tmp1
					
					//DBL velX = sqrt(tmp.m_vel.m_x*tmp.m_vel.m_x + tmp.m_vel.m_y*tmp.m_vel.m_y);
					
					DBL velX = -tmp.m_vel.m_x*sin(testangle) + tmp.m_vel.m_y*cos(testangle);
					
					tmp.m_pos.m_x = node_2_x - dist;
					tmp.m_pos.m_y = node_2_y;

					node.SetPoint(Math::C2DPoint(pMesh->GetBorderNode(nBoundaryNode).x - dist, pMesh->GetBorderNode(nBoundaryNode).y)*sin_phy);
				

					ALOGI("DD", AllToString(m_time) + CString(_T(" | "))
						+ AllToString(pMesh->GetBorderNode(nBoundaryNode).x) + CString(_T(" | "))
						+ AllToString(pMesh->GetBorderNode(nBoundaryNode).y) + CString(_T(" | "))
						+ AllToString(tmp.m_pos.m_x) + CString(_T(" | "))
						+ AllToString(tmp.m_pos.m_y) + CString(_T(" | "))
						+ AllToString(tmp.m_vel.m_x) + CString(_T(" | "))
						+ AllToString(tmp.m_vel.m_y) + CString(_T(" | "))
						+ AllToString(velX) + CString(_T(" | "))
						+ AllToString("SymX") + CString(_T(" | "))
						+ AllToString(angle_1*57.3) + CString(_T(" | "))
						+ AllToString((dRes)) + CString(_T(" | "))
						+ AllToString(number_of_element) + CString(_T(" | "))
						+ AllToString(dist) + CString(_T(" | "))
						+ AllToString(dSigmaNormal_x) + CString(_T(" | "))
						+ AllToString(dSigmaNormal_y) + CString(_T(" | "))
						+ AllToString(dSigmaInt) + CString(_T(" | "))
						+ AllToString(dEpsilonInt) + CString(_T(" | "))
						+ AllToString(closep)
						);
			

					
					//bc->at(nBoundaryNode).setSymX(velX, dRes, angle_1); //dRes*sin_phy		
					///////////////////////////////////

					//double ttt = atan(testangle);
					bc->at(nBoundaryNode).setSymX(velX, dRes, testangle);
				}
			//}
		}
		
	}
	return true;
}
	
		

//***********************************************************************************

//	DBL tmp2 = closestCurve->GelPerpendicularLength(pMesh->GetBorderNode(nBoundaryNode));
//	DBL tmp3 = closestCurve->GelPerpendicularLength(pMesh->GetBorderNode(nBoundaryNode2));

//number_of_element++;
//Угол между касательной и поверхностью инструмента
//DBL alpha_1 = Math::ToAngle(closestCurve->GetTangent());

//Угол между нормалью и поверхностью инструмента
//DBL beta = Math::ToAngle(closestCurve->GetNormal());

//	
//	DBL tmp4 = closestCurve->Lenth(pMesh->GetBorderNode(nBoundaryNode));
//	DBL tmp3 = tmp2 / tmp1;

//	DBL length = abs((node_2_x - node_1_x) - ((node_2_y - node_1_y))) /
//	(sqrt((node_2_x - node_1_x)*(node_2_x - node_1_x) + (node_2_y - node_1_y)*(node_2_y - node_1_y)));

//DBL l = closestCurve->GetClosestPoint(pMesh->GetBorderNode(nBoundaryNode), minim);
//	DBL dist;

//int nNode = m_shape.GetClosestNode(pMesh->GetBorderNode(nBoundaryNode), dist);
//DBL closep = m_shape.GetNode(nNode)->GetPoint();

//	DBL s = IsInside(pMesh->GetBorderNode(nBoundaryNode), m_shape.GetNode(nNode)->GetPoint(), EPS);
//CalcLength();

//DBL N = m_shape.GetClosestNode(pMesh->GetBorderNode(nBoundaryNode), length);

//*// cos_phy//
//cos(-beta + angle_2 +M_PI/2);
/*sin_phy;
sin(-beta + angle_2 + M_PI / 2)
*/
/*
if (abs(node_2_x - node_1_x) < abs(node_2_y - node_1_y))
{
phy = angle_2;

if (node_2_y - node_1_y < 0)
{
phy = 2 * M_PI - angle_2;
}
}
else
{
phy = angle_1;

if (node_2_x - node_1_x < 0)
{
phy = M_PI - angle_1;

if (phy < 0)
{
phy = 2 * M_PI + angle_1;
}
}


}*/
/*
if (0 < phy < M_PI / 2)
{
phy = M_PI / 2 + phy;
}
else if(M_PI/2 < phy < M_PI)
{
phy = 3*M_PI / 2 - phy;
}
else if (M_PI  < phy < 3*M_PI/2)
{
phy = abs(M_PI - phy);
}*/
//else if (3*M_PI / 2 < phy < 2*M_PI)
//		{
//		phy =  ;
//}



//-(-beta + angle_2) - угол нормали
//-(-beta + angle_2) - M_PI/2 -угол касательной

			//		if ((-beta + angle_2 >=M_PI/2) && angle_2>= M_PI / 2)
				//	{
					//	bc->at(nBoundaryNode).setSymX(velX, -dRes, 2 * M_PI - (M_PI + angle_2 - beta));
				//	}
				 //   else if ((-beta + angle_2 < M_PI/2) && angle_2 >= M_PI / 2)
				//	{
				//		bc->at(nBoundaryNode).setSymX(velX, -dRes, -beta + angle_2 );
				 //   }
				//	else if ((-beta + angle_2 >= M_PI / 2) && angle_2 < M_PI / 2)
				//	{
				//		bc->at(nBoundaryNode).setSymX(velX, -dRes, );
				//	}
				//	else if ((-beta + angle_2 < M_PI / 2) && angle_2 < M_PI / 2)
				//	{
				//		bc->at(nBoundaryNode).setSymX(velX, -dRes, );
				//	}



Math::C2DRect C2DRigid::GetBoundingBox()
{
	Math::C2DRect rect;

	rect.AddRect(m_shape.GetBoundingBox());
	rect.AddRect(m_motion.GetBoundingBox());
	
	return rect;
}

bool C2DRigid::GetBoundingBox(CRect2D &rect){
	rect.AddRect(m_shape.GetBoundingBox());
	rect.AddRect(m_motion.GetBoundingBox());
	return true;
}

void C2DRigid::DrawGL(GLParam &parameter){ 
	
	glColor3ub(138, 51, 36);
	m_shape.DrawGL(parameter);
	
	//m_motion.DrawGL(parameter);
}


void C2DRigid::DrawGL3D(GLParam &parameter) {

	//m_shape.DrawGL(parameter);
	DBL j, z, z1;
	DBL  step;

	DBL angle = parameter.GetAngle(); //угол поворота
	int n = parameter.GetStep3D();   // число разбиений
	step = angle / DBL(n);            // "шаг" отрисвки
									  //DBL angle_1;
									  //angle_1= angle*0;

									  //glColor3ub(138, 51, 36);
									  //m_shape.DrawGL(parameter);//Отрисовываем первую заглушку
									  //m_motion.DrawGL(parameter);
	if (parameter.IsDrawInstrument()) {

		if (m_shape.GetContourCount())
		{
			//проверка на наличие нулевого контура , если есть , то заполняем кэш
			C2DContour *pContour = m_shape.GetContour(0);
			if (!pContour->IsCache())
			{
				pContour->FillCache();
			}
			//дополнительно сохраняем кэш
			size_t f = pContour->GetCache().size();

			glColor3ub(138, 51, 36);
			//m_shape.DrawGL(parameter);
			glBegin(GL_POLYGON);
			for (size_t i = 0; i < f; i++)
				glVertex2dv(pContour->GetCache()[i]);
			glEnd();

			//отрисовка второй заглушки
			glMatrixMode(GL_MODELVIEW);
			//glLoadIdentity();
			glPushMatrix();
			GLfloat gl_angle = float(-angle);
			GLfloat	gl_x = 0.0;
			GLfloat gl_y = 1.0;
			GLfloat gl_z = 0.0;
			glRotatef(gl_angle, gl_x, gl_y, gl_z);
			glBegin(GL_POLYGON);
			for (size_t i = 0; i < f; i++)
				glVertex2dv(pContour->GetCache()[i]);
			glEnd();

			glPopMatrix();

			//glLoadIdentity();

			//циклы для отрисовки инструмента - аналогично C2DPlaneFEM

			glClearDepth(1.0f);
			glEnable(GL_DEPTH_TEST);
			glDepthFunc(GL_LEQUAL);
			glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

			for (size_t i = 0; i < f; i++)
			{
				Math::C2DPoint node = pContour->GetCache()[i];
				Math::C2DPoint node1 = pContour->GetCache()[(i + 1) % f];

				for (j = 0;j<n;j++) {

					z = j*step * M_PI_180;
					z1 = ((j + 1)*step) * M_PI_180;

					glBegin(GL_QUADS);              //отрисовываем полигоны  квадратами
					glColor3d(0.7, 0.7, 0.7);
					glVertex3d((node.x)*cos(z), (node.y), (node.x)*sin(z));
					glVertex3d((node1.x)*cos(z), (node1.y), (node1.x)*sin(z));
					glVertex3d((node1.x)*cos(z1), (node1.y), (node1.x)*sin(z1));
					glVertex3d((node.x)*cos(z1), (node.y), (node.x)*sin(z1));

					glEnd();

					//glLineWidth(0.7f);
					//glColor3d(0,0,0);
					//glBegin(GL_LINES);  //отрисовываем черным линии уровня //стыки между полигонами
					//  glVertex3d( (node.x)*cos(z),  (node.y),(node.x)*sin(z));
					//  glVertex3d( (node.x)*cos(z1),  (node.y),(node.x)*sin(z1));
					//   glEnd();
				}
			}

			for (size_t i = 0; i < m_shape.GetContour(0)->GetCurveCount(); i++)
			{
				
				Math::C2DPoint node = m_shape.GetContour(0)->GetCurve(i)->GetStartNode()->GetPoint();

				for (j = 0; j<n; j++) {

					z = j*step * M_PI_180;
					z1 = ((j + 1)*step) * M_PI_180;

					glLineWidth(0.7f);
					glColor3d(0, 0, 0);
					glBegin(GL_LINES);  //отрисовываем черным линии уровня //стыки между полигонами
					glVertex3d((node.x)*cos(z), (node.y), (node.x)*sin(z));
					glVertex3d((node.x)*cos(z1), (node.y), (node.x)*sin(z1));
					glEnd();

				}
			}

		}	
	}
}
//------------------------------------------------------	
//bool C2DRigid::IsDrawInstrument() const
//{
//	return m_DrawInstrument;
//}
//
//void C2DRigid::SetDrawInstrument(bool draw_inst)
//{
//	m_DrawInstrument = draw_inst;
//}
//-------------------------------------------------------


void C2DRigid::Preparations(const ITask *task){

}

C2DRigid::~C2DRigid()
{
}

/**************************
	Трение / Friction
*************************/

/*	fric_type - тип закона трения
	1 - Кулона-Амонтона
	2 - Зибеля
	3 - Леванова-Колмогорова
//*/
DBL C2DRigid::ChooseFrictionMethod(short fric_type, DBL dFrictCoeff, DBL dSigma_s, DBL dNormalPressure){
	
	switch(fric_type){
		case 1: return FrictLawKulon(dFrictCoeff, dNormalPressure); break;
		case 2: return FrictLawZibel(dFrictCoeff, dSigma_s); break;
		case 3: return FrictLawLevanov(dFrictCoeff, dSigma_s, dNormalPressure); break;
		default: return 0.0;
	}
}

DBL C2DRigid::FrictLawKulon(DBL dFrictCoeff, DBL dNormalPressure){
	
	// F = mu*P
	// P - нормальное давление
	DBL Ftr = dFrictCoeff*dNormalPressure;
	return Ftr;
}

DBL C2DRigid::FrictLawZibel(DBL dFrictCoeff, DBL dSigma_s){
	
	// F = mu*sigma_s 
	// sigma_s - интенсивность напряжений
	DBL Ftr = dFrictCoeff*dSigma_s;
	return Ftr;
}

//всё по y
DBL C2DRigid::FrictLawLevanov(DBL dFrictCoeff, DBL dSigma_s, DBL dNormalPressure){
	
	// F = K*(sigma_s/sqrt(3))*(1 - EXP(-1.25*P/sigma_s))
	// K - фактор трения

	if (fabs(dSigma_s) < EPS) return 0.0;	//чтобы не \ на 0
	DBL Ftr = dFrictCoeff*(dSigma_s/sqrt(3))*(1.0 - exp(-1.25*dNormalPressure/dSigma_s));
	return Ftr;
}


DBL C2DRigid::Friction(DBL dSn_x, DBL dSn_y, DBL dS_int, DBL dSquare){
			
		int nType = 3; //тип трения	
		
		//Трение
		DBL dSigmaNormal_x = dSn_x,	//нормальные напряжения
			dSigmaNormal_y = dSn_y,
			dSigmaInt = dS_int;		//интенсивность напряжения

		//Пропускаем первый шаг
		if (fabs(m_time) < EPS){
			dSigmaNormal_x = 0.0;
			dSigmaNormal_y = 0.0;
			dSigmaInt = 0.0;
		}
		
		// Пока считаем
		DBL dFrictForce_x = ChooseFrictionMethod(nType,m_FrictionCoeff,dSigmaInt,dSigmaNormal_x),
			dFrictForce_y = ChooseFrictionMethod(nType,m_FrictionCoeff,dSigmaInt,dSigmaNormal_y);
				
	    DBL dRes = dFrictForce_y*dSquare ;	//Если давим сверху, то сила трения по OY
		//DBL dRes = dFrictForce_y;
				
		//Указываем направление
		
		if (dSigmaNormal_y < 0.0){ 
			dRes = abs(dRes);
		}else{
			dRes = (-1)*abs(dRes);
		}
		
	return dRes;
}

void C2DRigid::WriteBCToLog() {

	
}