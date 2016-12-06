#ifndef _XTL_ERROR_H
#define _XTL_ERROR_H

#include <iostream>

#include "xtl_string.h"

namespace xtd {

	const bool good = false;			//по-умолчания false
	const bool no_good = true;		// по-умолчанию true

	struct error {
		static bool report(const std::string & error_message) {
			std::cerr << error_message << std::endl;
			m_status = no_good;
		
			return m_status;
		}
	
		static bool occured() { 
			return m_status; 
		}
		static bool did_not_occure() { 
			return !m_status; 
		}

		//сброс заранее сохранённой ошибки, теперь "ошибок нет" (good == false)
		static void reset() { 
			m_status = good;
		}	

	private:
		static bool m_status;
	};

	struct warning { static void report(const std::string & warning_message) { std::clog << warning_message << std::endl; } };

	struct success { static bool report() { return good; } };

}

#endif