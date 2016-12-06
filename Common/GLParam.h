#pragma once
#include "2DPolygonBuilder.h"
#include "2DScaleVisualizer.h"

class GLParam{

protected:
	C2DScale *m_Scale;		//шкала
	C2DScale m_EmptyScale;	//пустая шкала ?
	int m_eField;			//номер поля для отображения

	DBL m_Angle;	//угол
	int m_Step3D;	// качество осесимметричной поверхности = количество разбиений угла m_Angle
	bool m_DrawInstrument;	//отображать ли инструмент
	bool m_FillColor;	//отображать ли цвета

public:	
	C2DPolygonBuilder m_PolygonBuilder;
	CScaleVisualizer m_ScaleVisualizer;

public:
	GLParam(void);

	GLParam(const GLParam &parameter);

	GLParam& operator = (const GLParam &parameter);

	C2DScale* GetScale();
	void SetScale(C2DScale &scale);

	int GetField(){
		return m_eField;
	}

	void SetField(int eField){
		m_eField = eField;
	}

	int GetStep3D()
	{
		return m_Step3D;
	}

	void SetStep3D(int step3d);

	DBL GetAngle()
	{
		return m_Angle;
	}

	void SetAngle(DBL angle);

	bool IsFillColor() const
	{
		return m_FillColor;
	}

	void SetFillColor(bool fill_color)
	{
		m_FillColor = fill_color;
	}


	bool IsDrawInstrument() const
	{
		return m_DrawInstrument;
	}

	void SetDrawInstrument(bool draw)
	{
		m_DrawInstrument = draw;
	}

};

