#pragma once

#include "xtl/xtl_memory.h"
#include "xtl/xtl_functional.h"
#include "xtl/xtl_double_dispatcher.h"

#include "Unit.h"
#include "Arithmetic.h"

namespace MM {

	//!		Диспетчер бинарных функций.
	/*!
	Выполняет двойную диспетчеризацию бинарных функций, оперирующих с указателями на "CUnit".
	Для этого в нём регистрируются пары, состоящие из пары строк — идентификаторов оперируемых элементов — и конкретной функции,
	которая должна быть вызвана, если придёт запрос на соответствующую пару идентификаторов.
	*/
	template <
		template <
		typename FirstArgumentType,
		typename SecondArgumentType,
		typename ResultType
		> class BinaryPredicate
	>
	class CUnitDoubleDispatcher {
		typedef xtd::smart_ptr<CUnit> CUnitSptr;

	public:
		static CUnitSptr Dispatch(const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) {
			return Instance().DoDispatch(leftOperand, rightOperand);
		}

		template <typename FirstType, typename SecondType, typename ResultType>
		static void Put() { Instance().DoPut<FirstType, SecondType, ResultType>(); }

		template <typename Type>
		static void Put() { Instance().DoPut<Type>(); }

		CUnitDoubleDispatcher& operator = (const CUnitDoubleDispatcher& p) {
			m_dispatcher = p.m_dispatcher;
			instance = p.instance;
			return *this;
		};

	private:
		CUnitDoubleDispatcher() { Init(); }
		CUnitDoubleDispatcher(const CUnitDoubleDispatcher &);

		CUnitSptr DoDispatch(const CUnitSptr & leftOperand, const CUnitSptr & rightOperand) const {
			return m_dispatcher.dispatch(leftOperand, rightOperand);
		}

		template <typename FirstType, typename SecondType, typename ResultType>
		void DoPut() { m_dispatcher.template put<FirstType, SecondType, ResultType>(); }

		template <typename Type>
		void DoPut() { m_dispatcher.template put<Type>(); }

		//!		Инициализация.
		/*!
		В этой функции регистрируются разрешённые по-умолчанию операции с объектами.
		Дополнительные операции можно зарегистрировать извне с помощью вызова
		"CUnitDoubleDispatcher<операция>::Put<список_типов>();"
		*/
		void Init() { /* инициализации должны быть специализированы */ }

		static CUnitDoubleDispatcher & Instance() { return instance; }

	private:
		typedef xtd::double_dispatcher<CUnit, BinaryPredicate> DispatcherType;
		DispatcherType m_dispatcher;

		static CUnitDoubleDispatcher instance;
	};

} // namespace MM