#pragma once

#include "xtl/xtl_memory.h"
#include "xtl/xtl_functional.h"
#include "xtl/xtl_single_dispatcher.h"

#include "Unit.h"
#include "Arithmetic.h"

namespace MM {

	//!		Диспетчер унарных функций.
	/*!
	Выполняет диспетчеризацию унарных функций, оперирующих с указателем на "CUnit".
	Для этого в нём регистрируется пара, состоящая из идентификатора оперируемого элемента и конкретной функции,
	которая должна быть вызвана, если придёт запрос на соответствующий идентификатор.
	*/
	template <
		template <
		typename ArgumentType,
		typename ResultType
		> class UnaryPredicate
	>
	class CUnitSingleDispatcher {
		typedef xtd::smart_ptr<CUnit> CUnitSptr;

	public:
		static CUnitSptr Dispatch(const CUnitSptr & operand) { return Instance().DoDispatch(operand); }

		template <typename Type, typename ResultType>
		static void Put() { Instance().DoPut<Type, ResultType>(); }

		template <typename Type>
		static void Put() { Instance().DoPut<Type>(); }

		CUnitSingleDispatcher& operator = (const CUnitSingleDispatcher& p) {
			m_dispatcher = p.m_dispatcher;
			instance = p.instance;
			return *this;
		};

	private:
		CUnitSingleDispatcher() { Init(); }
		CUnitSingleDispatcher(const CUnitSingleDispatcher &);

		CUnitSptr DoDispatch(const CUnitSptr & operand) const { return m_dispatcher.dispatch(operand); }

		template <typename Type, typename ResultType>
		void DoPut() { m_dispatcher.template put<Type, ResultType>(); }

		template <typename Type>
		void DoPut() { m_dispatcher.template put<Type>(); }

		//!		Инициализация.
		/*!
		В этой функции регистрируются операции с объектами, разрешённые по-умолчанию.
		Дополнительные операции можно зарегистрировать извне с помощью вызова
		"CUnitSingleDispatcher<операция>::Put<список_типов>();"
		*/
		void Init() { /* инициализации должны быть специализированы */ }

		static CUnitSingleDispatcher & Instance() { return instance; }

	private:
		typedef xtd::single_dispatcher<CUnit, UnaryPredicate> DispatcherType;
		DispatcherType m_dispatcher;

		static CUnitSingleDispatcher instance;
	};

} // namespace MM