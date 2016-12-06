#pragma once
#include "I2DInteractiveObject.h"
#include "AllTypes.h"
#include "2DMotion.h"
#include "2DBCAtom.h"
#include "2DPoint.h"


//! Класс сборки Инструментов и Траекторий
class C2DRigid : public I2DInteractiveObject
{
	GETTYPE (T2DRIGID)
public:
	DBL m_time;					// Время прибытия в узел
	Math::C2DPoint m_oldpos;	// Предыдущая позиция

	C2DOutline m_shape; //< Форма инструмента
	IO::C2DPoint m_stickpoint; //< Точка привязки формы
	IO::CDouble m_FrictionCoeff; // Коэффициент трения формы
	CIO2DMotion m_motion; //< Движение формы

	bool m_IsStuck;			//Есть ли контакт с КЭ сеткой (трение) (init==false)
	bool m_IsFirstStuck;	//Первичен ли контакт с КЭ сеткой (сцепка/заделка) (init==true)

public:
	C2DRigid();

	void SetShape(C2DOutline* pshape)
	{
		m_shape.Copy(dynamic_cast<IO::CInterface*>(pshape));
	}

	void SetMotion(CIO2DMotion* pmotion)
	{
		m_motion.Copy(dynamic_cast<IO::CInterface*>(pmotion));
	}

	void SetStickPoint(const Math::C2DPoint &point){
		m_stickpoint() = point;
	}

	void SetFrictionCoeff(DBL coeff){
		m_FrictionCoeff() = coeff;
	}

	void ResetTime();

	//virtual bool Load(CStorage& file);
	//virtual bool Save(CStorage& file);

	virtual bool Init();
	virtual void Calc(double dt);
	virtual void Move(double dt);

	virtual void Preparations(const ITask *task);	// подготовка объектов для расчёта

	// 1 - точка point внутри или на границе, -1 - не внутри. closest - ближайшая, len - точность.
	// не учитывает углы и другие случаи
	virtual bool IsInside(const Math::C2DPoint& point, Math::C2DPoint& closest, DBL len);	
	
	virtual bool GetBC(const C2DMeshInterface *pMesh, size_t nBoundaryNode, C2DBCAtom& bc);
	virtual bool GetBC(const C2DMeshInterface *pMesh, std::vector<C2DBCAtom> *bc);

	virtual Math::C2DRect GetBoundingBox();
	virtual bool GetBoundingBox(CRect2D &rect);

	virtual void DrawGL(GLParam &parameter);
	virtual void DrawGL3D(GLParam &parameter);
 
	//ГУ трение
	DBL Friction(DBL dSn_x, DBL dSn_y, DBL dS_int, DBL dSquare);
	DBL ChooseFrictionMethod(short fric_type, DBL dFrictCoeff, DBL dSigma_s, DBL dNormalPressure);
	DBL FrictLawKulon(DBL dFrictCoeff, DBL dNormalPressure);	//1 Kulon
	DBL FrictLawZibel(DBL dFrictCoeff, DBL dSigma_s);			//2 Zibel
	DBL FrictLawLevanov(DBL dFrictCoeff, DBL dSigma_s, DBL dNormalPressure);	//3 Levanov


	virtual ~C2DRigid();
};

IODEF (C2DRigid)

