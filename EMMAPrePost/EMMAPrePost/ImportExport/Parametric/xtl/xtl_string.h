#ifndef _XTL_STRING_H
#define _XTL_STRING_H

#include <string>

namespace xtd {

inline bool contains (const std::string & string, char c) {
	return string.find_first_of(c) != std::string::npos;
}

inline bool contains_all_of (const std::string & string, const std::string & characters) {
	return characters.find_first_not_of(string) == std::string::npos;
}

inline bool is_shorter_or_less (const std::string & s1, const std::string & s2) {
	if (s1.size() < s2.size()) return true;
	else if (s1 < s2) return true;
	else return false;
}

}

#endif