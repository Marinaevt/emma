#pragma once

#include <string>
#include <vector>
#include <algorithm>

#include "Unit.h"
#include "../../../../Common/2DOutline.h"

namespace MM {

	//!		Геометрический элемент.
	/*!
	Базовый класс для всех геометрических элементов метамодели.
	*/
	class CElement : public CUnit {
	public:
		CElement() { Init(); }

		// ---- Модификация ---------------------------------------------------------------------------
	
		void SetName(const std::string & name) { m_name = name; }

		// ---- Доступ --------------------------------------------------------------------------------
	
		const std::string & Name() const { return m_name; }
	
		virtual int ExportToOutline(C2DOutline * outline) const;
	
		// ---- Индикаторы ----------------------------------------------------------------------------

		virtual bool NameIs(const std::string & string) const { return m_name == string; }
	
		// ---- Ввод-вывод ----------------------------------------------------------------------------
	
		typedef std::vector<CUnitSptr> parametersContainer;
	
		bool			  Load(CReader & reader) = 0;
		bool              Load(CReader & reader, parametersContainer & parameters);
		bool          LoadName(CReader & reader);
		bool PrepareForLoading(CReader & reader) const;
		bool     FinishLoading(CReader & reader);
		bool       UndoLoading(CReader & reader, parametersContainer & parameters);
	
		void SaveBeginning(std::ostream & stream) const;
		void SaveDelimeter(std::ostream & stream) const;
		void    SaveEnding(std::ostream & stream) const;
	
		void SaveAsExpression(std::ostream & stream) const;
	
	private:
		void Init() { m_name.clear(); }

	protected:
		std::string m_name; //! Имя, под которым элемент фигурирует в исходном файле.
	};

} // namespace MM