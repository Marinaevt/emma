#pragma once

#include <vector>

#include "xtl/xtl_stack.h"
#include "xtl/xtl_factory.h"
#include "xtl/xtl_token.h"

#include "Unit.h"

namespace MM {

	//!		Загрузчик.
	/*!
	Обеспечивает лексический анализ входного файла, разбивая его на лексемы, предоставляет доступ
	к этим лексемам, а также является хранилищем уже считанных, ожидающих и преостановленных юнитов в
	процессе загрузки метамодели.
	*/
	class CReader {
	public:
		CReader();
		CReader(const std::string & string);
	
		// ---- Чтение элементов ----------------------------------------------------------------------
	
		bool GetNext();		//чтение и разбор строки посимвольно, сравнение с вшитыми "лексемами"
		bool Read();
		bool Read(const std::string & string);
	
		bool HasToken() const { return m_currentPos != m_tokens.end(); }
		std::string PeekToken() const { return *m_currentPos; }
		std::string  GetToken() { return *m_currentPos++; }
	
		// ---- Модификация ---------------------------------------------------------------------------
	
		void Reset(const std::string & string);
	
		template <typename Type>
		void Put() {
			//m_tokenizer.add_recognizable_token(Type().Definition());
			m_factory.put<Type>();
		}
	
		void     AddReadyUnit(const CUnitSptr & unit) { m_readyUnits.push(unit); }
		void   AddWaitingUnit(const CUnitSptr & unit) { m_waitingUnits.push(unit); }
		void AddSuspendedUnit(const CUnitSptr & unit) { m_suspendedUnits.push(unit); }
	
		void AddReadyElement(const CUnitSptr & element) { m_readyElements.push(element); }
	
		// ---- Доступ --------------------------------------------------------------------------------
	
		const CUnitSptr &     GetReadyUnit() const { return     m_readyUnits.top(); }
		const CUnitSptr &   GetWaitingUnit() const { return   m_waitingUnits.top(); }
		const CUnitSptr & GetSuspendedUnit() const { return m_suspendedUnits.top(); }
	
		const CUnitSptr & GetReadyElement() const { return m_readyElements.top(); }
	
		CUnitSptr     ExtractReadyUnit() { return     m_readyUnits.pop(); }
		CUnitSptr   ExtractWaitingUnit() { return   m_waitingUnits.pop(); }
		CUnitSptr ExtractSuspendedUnit() { return m_suspendedUnits.pop(); }
	
		CUnitSptr ExtractReadyElement() { return m_readyElements.pop(); }
	
		bool     HasReadyUnit() const { return     !m_readyUnits.empty(); }
		bool   HasWaitingUnit() const { return   !m_waitingUnits.empty(); }
		bool HasSuspendedUnit() const { return !m_suspendedUnits.empty(); }
	
		bool HasReadyElement() const { return !m_readyElements.empty(); }
	
		//Возвращает текущее название токена
		const std::string & CurrentToken() const { return m_currentToken; }
	
		const xtd::stack<CUnitSptr>::container & GetReadyElements() const { return m_readyElements.genElNumlements(); }
	
	private:
		xtd::tokenizer m_tokenizer;
	
		typedef std::vector<std::string> tokensContainer;
		tokensContainer m_tokens;
		tokensContainer::const_iterator m_currentPos;
		std::string m_currentToken;		//текущее название токена
	
		xtd::factory<CUnit> m_factory;
	
		xtd::stack<CUnitSptr>     m_readyUnits;
		xtd::stack<CUnitSptr>   m_waitingUnits;
		xtd::stack<CUnitSptr> m_suspendedUnits;
	
		xtd::stack<CUnitSptr> m_readyElements; // Сюда записываются только элементы, не выражения.
	};

} // namespace MM