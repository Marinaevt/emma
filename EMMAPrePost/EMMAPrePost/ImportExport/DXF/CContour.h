#pragma once

#include <vector>
#include "../../../../Common/2DOutline.h"
#include "CGraph.h"
#include "CCurve.h"
#include "CReader.h"

namespace DXF {

	class CContour;
	typedef std::vector<CContour> CContourArray;

	class CGraph;

	//! Контур.
	/*!
	Содержит массив указателей на разнообразные кривые, принадлежащие этому контуру.
	*/
	class CContour {
	public:
		CContour() {}
		~CContour() {}
	
		bool LoadFromDXF(CReader & reader);
		bool SaveToDXF(CReader & reader);
	
		void AddCurve(CCurve * curve) { m_curves.push_back(curve); }
	
		CCurvePtrArray::iterator CurvesBegin() { return m_curves.begin(); }
		CCurvePtrArray::iterator CurvesEnd() { return m_curves.end(); }
	
		CCurvePtrArray::const_iterator CurvesBegin() const { return m_curves.begin(); }
		CCurvePtrArray::const_iterator CurvesEnd() const { return m_curves.end(); }
	
		//! Регистрация кривых, содержащихся в контуре, в чертеже.
		/*!
		Регистрирует точки, содержащиеся в кривых из массива "m_curves", в чертеже "graph" и
		проверяет каждую из этих кривых на наличие в чертеже. Если кривая в чертеже уже есть, то
		указатель на неё освобождается и меняется на указатель на аналогичную кривую из чертежа.
		В противном случае в чертёж добавляется новая кривая.
		*/
		void RegisterCurvesInGraph(CGraph & graph);
	
		bool operator == (const CContour & contour) const;
	
		int ExportToOutline(C2DOutline * outline) const;

	private:
		CCurvePtrArray m_curves;
	};

} // namespace DXF
