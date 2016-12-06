#pragma once
#include <vector>
#include "IOInterface.h"
#include "Color.h"	//class 2DColor
#include "ColorScale.h" // class 2DScaleColor

//---------- class ICallbackable ------------------
class ICallbackable{
public:
	virtual void Callback() = 0;
};

//---------- class C2DScale ------------------
class C2DScale : public IO::CClass
{
	GETTYPE(C2DSCALE)
public:
	C2DScale();
	C2DScale(size_t count, bool isDescrete, bool isLog);
	C2DScale(const C2DScale &scale);	
	virtual ~C2DScale() {}

	virtual bool Load(CStorage& file);
	bool SaveToPath(CString& path);

	bool IsInit() const;
	bool IsDescrete() const;
	bool IsLog() const;
	bool IsDrawLines() const;
	void SetDrawLines(bool draw);
	CString GetStringValue(size_t i, bool onlyFirstValue = false) const;	///Возвращает значение или диапазон значений в виде CString
	
	void Set(size_t count, bool isDescrete, bool isLog);					/// Устанавливает шкалу
	C2DColor GetColor(double value) const;								// Цвет по значению
	C2DColor GetColor(size_t index) const;								// Цвет по номеру значения в шкале
	void SetColor(size_t index, const C2DColor &color);

	bool IsScaleMark(double value) const;								// Определяет является ли значение отметкой на шкале
	void AddMark(double value);											// Добавляет отметку
	DBL GetValue(size_t index) const;										// Возвращает значение по порядковому номеру
	size_t GetCount() const;											// Кол-во значений в шкале. Всегда не меньше 2-ух
	void SetBounds(double minValue, double maxValue);					// Максимальное и минимальное значение элементов детали. Для GetStringValue
	C2DColor CalculateDefaultColor(double value) const;
	void SetCallback(ICallbackable *callbackable);

	void operator = (const C2DScale &scale);

protected:
	//void Split();													
	void ArrangeColors();									// Расставляет диапазон действия у цветов, если шкала дискретная
	C2DColor GetAppropriateColor(double value) const;		// Используется в случае непрерывной шкалы
	void RegisterMembers();

	double m_MinValue;
	double m_MaxValue;
	bool m_IsInit;
	bool m_DrawLines;
	ICallbackable *m_Callbackable;

// Далее идут поля, которые сохраняются и загружаются
public:
	IO::CString m_Name;	

protected:
	IO::CSimpleArray<C2DScaleColor, C2DSCALECLRRARRAY> m_Colors;
	IO::CBool m_IsDescrete;	
	IO::CBool m_IsLog;
};

IODEF(C2DScale)