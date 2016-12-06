#pragma once

#include <iostream>
#include "Element.h"
#include "Point.h"
#include "../DXF/CPoint.h"
#include "Operands.h"
#include "../../../../Common/2DOutline.h"

namespace MM {

	class CCurve : public CElement {
	public:
		bool Load(CReader & reader) { return xtd::error::report("Функция пока не работает"); }
	};

	class CLine : public CCurve {
	public:
		CLine() : m_startPoint(new CPoint), m_endPoint(new CPoint) { Init(); }
		CLine(const CPoint & startPoint, const CPoint & endPoint) :
			m_startPoint(new CPoint(startPoint)), m_endPoint(new CPoint(endPoint)) {
			Init();
		}
		CLine(const CUnitSptr & startPoint, const CUnitSptr & endPoint) : m_startPoint(startPoint), m_endPoint(endPoint) { Init(); }
	
		// ---- Доступ к переменным -------------------------------------------------------------------
	
		CPoint StartPoint() const { return dynamic_cast<CPoint &>(*m_startPoint->Result()); }
		CPoint   EndPoint() const { return dynamic_cast<CPoint &>(*m_endPoint->Result()); }
	
		CUnitSptr Result() const { return new CLine(m_startPoint->Result(), m_endPoint->Result()); }
	
		int ExportToOutline(C2DOutline * outline) const;
	
		// ---- Ввод-вывод ----------------------------------------------------------------------------
	
		bool Load(CReader & reader);
		void Save(std::ostream & stream) const;
	
	private:
		void Init() { m_definition = "line"; }

		CUnitSptr m_startPoint; //! Указатель на выражение, результат которого обязан указывать на экземлпляр класса CPoint.
		CUnitSptr m_endPoint;   //! Указатель на выражение, результат которого обязан указывать на экземлпляр класса CPoint.
	};

	class CCircle : public CCurve {
	public:
		CCircle() : m_center(new CPoint), m_radius(new CNumber) { Init(); }
		CCircle(const CPoint & center, const CNumber & radius) : m_center(new CPoint(center)), m_radius(new CNumber(radius)) { Init(); }
		CCircle(const CUnitSptr & center, const CUnitSptr & radius) : m_center(center), m_radius(radius) { Init(); }
	
		// ---- Доступ к переменным -------------------------------------------------------------------
	
		CPoint  Center() const { return dynamic_cast<CPoint  &>(*m_center->Result()); }
		CNumber Radius() const { return dynamic_cast<CNumber &>(*m_radius->Result()); }
	
		CUnitSptr Result() const { return new CCircle(m_center->Result(), m_radius->Result()); }

		//TODO: пока пусто (загрузка кругов)
		//int ExportToOutline(C2DOutline * outline) const;

		// ---- Ввод-вывод ----------------------------------------------------------------------------
	
		bool Load(CReader & reader);
		void Save(std::ostream & stream) const;

	private:
		void Init() { m_definition = "circle"; }

		CUnitSptr m_center; //! Указатель на выражение, результат которого обязан указывать на экземпляр класса CPoint.
		CUnitSptr m_radius;	//! Указатель на выражение, результат которого обязан указывать на экземпляр класса CNumber.
	};

	class CArc : public CCurve {
	public:
		CArc() : m_startPoint(new CPoint), m_endPoint(new CPoint), m_radius(new CNumber) { Init(); }
		CArc(const CPoint & startPoint, const CPoint & endPoint, const CNumber & radius) :
			m_startPoint(new CPoint(startPoint)), m_endPoint(new CPoint(endPoint)), m_radius(new CNumber(radius)) {
			Init();
		}
		CArc(const CUnitSptr & startPoint, const CUnitSptr & endPoint, const CUnitSptr & radius) :
			m_startPoint(startPoint), m_endPoint(endPoint), m_radius(radius) {
			Init();
		}
		
		// ---- Доступ --------------------------------------------------------------------------------
	
		CPoint  StartPoint() const { return dynamic_cast<CPoint  &>(*m_startPoint->Result()); }
		CPoint    EndPoint() const { return dynamic_cast<CPoint  &>(*m_endPoint->Result()); }
		CNumber     Radius() const { return dynamic_cast<CNumber &>(*m_radius->Result()); }
	
		CUnitSptr Result() const { return new CArc(m_startPoint->Result(), m_endPoint->Result(), m_radius->Result()); }
	
		int ExportToOutline(C2DOutline * outline) const;
	
		// ---- Ввод-вывод ----------------------------------------------------------------------------
	
		bool Load(CReader & reader);
		void Save(std::ostream & stream) const;

	private:
		void Init() { m_definition = "arc"; }
	
		CUnitSptr m_startPoint;
		CUnitSptr m_endPoint;
		CUnitSptr m_radius;
	};

} // namespace MM