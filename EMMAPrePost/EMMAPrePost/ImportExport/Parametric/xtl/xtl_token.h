#pragma once

#include <string>
#include <list>

#include "xtl_memory.h"

namespace xtd {

inline bool is_digit (char c) { return std::string("1234567890").find(c) != std::string::npos; }
inline bool is_letter (char c) { return std::string("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ").find(c) != std::string::npos; }
inline bool is_punctuation_mark (char c) { return std::string("!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~").find(c) != std::string::npos; }
inline bool is_space (char c) { return std::string(" \t\n\v\f\r").find(c) != std::string::npos; }
inline bool is_letter_or_digit (char c) { return is_letter(c) || is_digit(c); }

class tokenizer {
private:
	class abstract_rule {
	public:
		virtual ~abstract_rule () {}
		
		virtual bool allows_to_start_with (char c) const = 0;
		
		virtual size_t find_next_token (const std::string & string, size_t pos = 0) const = 0;
	};
	
	typedef smart_ptr<abstract_rule> abstract_rule_sptr;
	typedef std::list<abstract_rule_sptr> rules_container;
	
	class number_rule : public abstract_rule {
	public:
		bool allows_to_start_with (char c) const { return is_digit(c); }
		
		size_t find_next_token (const std::string & string, size_t pos = 0) const {
			size_t new_pos = std::string::npos;
			
			while (is_digit(string[pos])) ++pos;
			
			if (std::string("eE").find_first_of(string[pos]) != std::string::npos)
				new_pos = find_next_token_in_exponential_form(string, pos + 1);
			else if (string[pos] == '.')
				new_pos =  find_next_token_as_decimal(string, pos + 1);
			
			if (new_pos != std::string::npos) return new_pos;
			else return pos;
		}
		
	public:
		size_t find_next_token_in_exponential_form (const std::string & string, size_t pos) const {
			if (std::string("-+").find_first_of(string[pos]) != std::string::npos) ++pos;
			else return std::string::npos;
			
			if (is_digit(string[pos])) while (is_digit(string[++pos]));
			else return std::string::npos;
			
			return pos;
		}
		
		size_t find_next_token_as_decimal (const std::string & string, size_t pos) const {
			if (is_digit(string[pos])) while (is_digit(string[++pos]));
			else return std::string::npos;
			
			if (std::string("eE").find_first_of(string[pos]) != std::string::npos)
				return find_next_token_in_exponential_form(string, pos + 1);
			else return pos;
		}
	};
	
	class literal_rule : public abstract_rule {
	public:
		bool allows_to_start_with (char c) const { return is_letter(c); }
		
		size_t find_next_token (const std::string & string, size_t pos = 0) const {
			while (is_letter_or_digit(string[pos])) ++pos;
			
			return pos;
		}
	};
    
    class comment_rule : public abstract_rule {
    public:
        bool allows_to_start_with (char c) const { return c == '/' || c == '*'; }
        
        size_t find_next_token (const std::string & string, size_t pos = 0) const {
            if ((string[pos] == '/' && string[pos + 1] == '*')
                ||
                (string[pos] == '*' && string[pos + 1] == '/')
            ) return pos + 2;
            
            return pos;
        }
    };
	
	class punctuation_rule : public abstract_rule {
	public:
		bool allows_to_start_with (char c) const { return is_punctuation_mark(c) != 0; }
		
		size_t find_next_token (const std::string & string, size_t pos = 0) const {
			if (is_punctuation_mark(string[pos])) ++pos;
			
			return pos;
		}
	};
	
	class spacing_rule : public abstract_rule {
	public:
		bool allows_to_start_with (char c) const { return is_space(c); }
		
		size_t find_next_token (const std::string & string, size_t pos = 0) const {
			while (is_space(string[pos])) ++pos;
			
			return pos;
		}
	};
    
    class unrecognized_rule : public abstract_rule {
    public:
        bool allows_to_start_with (char c) const { return true; }
        
        size_t find_next_token (const std::string & string, size_t pos = 0) const {
            while (!is_space(string[pos])) ++pos;
            
            return pos;
        }
    };
	
public:
	tokenizer () { init(); }

	template <typename output_container>
	bool tokenize (const std::string & string, output_container & container) const {
		output_container tokens(container);
		
		for (size_t begin = 0, end = find_next_token(string); end != std::string::npos; begin = end, end = find_next_token(string, begin)) {
			tokens.push_back(string.substr(begin, end - begin));
		}
			
		if (tokens.size() == 0) {
			return false;
		}
		
		std::swap(tokens, container);
		return true;
	}
	
	size_t find_next_token (const std::string & string, size_t pos = 0) const {
		if (pos < string.size()) {
            rules_container::const_iterator i = next_rule_for_character(string[pos], m_rules.begin());
            
			while (i != m_rules.end()) {
                size_t new_pos = (*i)->find_next_token(string, pos);
                
                if (new_pos != pos) return new_pos;
                
                i = next_rule_for_character(string[pos], ++i);
            }
            
            xtd::warning::report("Символ \"" + std::string(1, string[pos++]) + "\" неизвестен токенайзеру. Лексема будет проигнорирована.");
            return unrecognized_rule().find_next_token(string, pos);
        }
		
		return std::string::npos;
	}
	
private:
	rules_container::const_iterator next_rule_for_character (char c, rules_container::const_iterator start) const {
		for (rules_container::const_iterator i = start; i != m_rules.end(); ++i)
			if ((*i)->allows_to_start_with(c)) return i;
		
		return m_rules.end();
	}
	
	void init () {
		m_rules.push_back(new number_rule());
		m_rules.push_back(new literal_rule());
        m_rules.push_back(new comment_rule());
		m_rules.push_back(new punctuation_rule());
		m_rules.push_back(new spacing_rule());
	}
	
private:
	rules_container m_rules;
};

}