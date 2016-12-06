#ifndef _XTL_MEMORY_H
#define _XTL_MEMORY_H

#include <memory>
#include <map>

namespace xtd {

//!		Умный указатель.
/*!
		Внешне выглядит и ведёт себя как обычный указатель, но на самом деле является контейнером для
	настоящего указателя и следит за внешними ссылками на этот указатель, автоматически освобождая
	память, когда этих ссылок становится ровно ноль.	
		Умеет работать только с одной иерархией, т.е. только с классом "type" и его наследниками. К другим
	типам не преобразуется.
*/
template <typename type>
class smart_ptr {
private:
	// ---- Определения ---------------------------------------------------------------------------
	
	class traits {
	public:
		traits (type * some_pointer) : pointer(some_pointer), reference_count(1) {}
		~traits () { delete pointer; }
		
		void retain () { ++reference_count; }
		bool release () { return --reference_count == 0; }
		
		type * get () const { return pointer; }
		
	private:	
		type * const pointer;
		int reference_count;
	};
	
public:
	// ---- Создание и уничтожение ----------------------------------------------------------------
	
	smart_ptr (type * some_pointer = nullptr) : m_traits(some_pointer ? get_traits_for_pointer(some_pointer) : nullptr) {}
	smart_ptr (const smart_ptr & some_smart_pointer) : m_traits(nullptr) { copy(some_smart_pointer); }
	~smart_ptr () { release(); }
	
	// ---- Модификация ---------------------------------------------------------------------------
	
	smart_ptr & operator = (const smart_ptr & some_smart_pointer) {
		copy(some_smart_pointer);
		
		return *this;
	}
	
	// ---- Операции сравнения --------------------------------------------------------------------
	
	bool operator == (const smart_ptr & some_smart_pointer) const { return m_traits == some_smart_pointer.m_traits; }
	bool operator != (const smart_ptr & some_smart_pointer) const { return m_traits != some_smart_pointer.m_traits; }
	
	// ---- Доступ --------------------------------------------------------------------------------
	
	type & operator * () const { return *m_traits->get(); }
	type * operator -> () const { return m_traits->get(); }

protected:
	// ---- Служебные функции ---------------------------------------------------------------------
	
	void retain () const { if (m_traits) m_traits->retain(); }
	
	void release () const {
		if (m_traits && m_traits->release()) {
			observer.erase(m_traits->get());
			
			delete m_traits;
		}
	}
	
	void copy (const smart_ptr & some_smart_pointer) {
		release();
		m_traits = some_smart_pointer.m_traits;
		retain();
	}
	
private:
	traits * get_traits_for_pointer (type * some_pointer) {
		typename observer_type::iterator i = observer.find(some_pointer);
		traits * new_traits;
		
		if (i == observer.end()) {
			new_traits = new traits(some_pointer);
			observer.insert(std::make_pair(some_pointer, new_traits));
		} else {
			new_traits = (*i).second;
			new_traits->retain();
		}
		
		return new_traits;
	}
	
private:
	traits * m_traits;
	
	// Ассоциативный массив, устанавливающий биекцию между указателем и числом ссылок на него, предотвращая
	// таким образом создание нескольких различных указателей на один и тот же объект.
	typedef std::map<type *, traits *> observer_type;
	static observer_type observer;
};

// Инициализация статической переменной.
template <typename type> typename smart_ptr<type>::observer_type smart_ptr<type>::observer;

}

#endif