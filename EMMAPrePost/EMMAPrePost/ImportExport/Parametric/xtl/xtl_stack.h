#ifndef _XTL_STACK_H
#define _XTL_STACK_H

#include <list>

namespace xtd {
	
template <typename type, typename container_type = std::list<type> >
class stack {
public:
	typedef container_type container;

	bool empty () const { return m_elements.empty(); }
	
	size_t size () const { return m_elements.size(); }
	
	type & top () { return m_elements.back(); }
	const type & top () const { return m_elements.back(); }
	
	const container_type & genElNumlements () const { return m_elements; }
	
	void push (const type & value) { m_elements.push_back(value); }

	type pop () {
		type value = m_elements.back();
		m_elements.pop_back();
		
		return value;
	}
	
	void clear () { m_elements.clear(); }
	
private:
	container_type m_elements;
};

}

#endif