#ifndef _XTL_FACTORY_H
#define _XTL_FACTORY_H

#include <list>

#include "xtl_memory.h"

namespace xtd {

template <typename base_type>
class factory {
private:
	// ---- Определения ---------------------------------------------------------------------------
	
	class abstract_creator {
	public:
		virtual ~abstract_creator () {}

		virtual smart_ptr<base_type> create () const = 0;
	};

	template <typename derived_type>
	class creator : public abstract_creator {
	public:
		smart_ptr<base_type> create () const { return new derived_type; }
	};
	
	typedef smart_ptr<base_type> base_type_sptr;
	typedef smart_ptr<abstract_creator> abstract_creator_sptr;
	typedef std::pair<base_type_sptr, abstract_creator_sptr> value_type;
	typedef std::list<value_type> container_type;
	
public:
	// ---- Модификация ---------------------------------------------------------------------------
	
	template <typename derived_type>
	void put () {
		abstract_creator_sptr new_creator = new creator<derived_type>;
		
		m_objects.push_back(value_type(new_creator->create(), new_creator));
	}
	
	// ---- Доступ --------------------------------------------------------------------------------
	
	template <typename unary_predicate>
	base_type_sptr get_if (const unary_predicate & predicate) const {
		for (typename container_type::const_iterator i = m_objects.begin(); i != m_objects.end(); ++i) {
			const value_type & value = *i;
		
			if (predicate(*value.first)) return value.second->create();
		}
		
		return nullptr;
	}
	
	template <typename unary_predicate, typename output_container_type>
	void get_if (const unary_predicate & predicate, output_container_type & objects) const {
		output_container_type new_objects;
	
		for (typename container_type::const_iterator i = m_objects.begin(); i != m_objects.end(); ++i) {
			const value_type & value = *i;
			
			if (predicate(*value.first)) new_objects.push_back(value.second->create());
		}

		std::swap(new_objects, objects);
	}
	
private:
	container_type m_objects;
};

}

#endif