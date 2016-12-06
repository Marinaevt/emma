#pragma once

#include <iostream>
#include <vector>

// Служебные функции.

//! Освобождает массив указателей.
template <class T>
void deleteVector (std::vector<T *> & vector) {
	for (typename std::vector<T *>::iterator i = vector.begin(); i < vector.end(); ++i) delete *i;
}

//! Записывает в консоль ошибки, предупреждения и прочие сообщение о работе программы.
inline void logMessage (const char * message) { std::cout << message << std::endl; }