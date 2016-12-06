#pragma once

#include "..\..\Common\SimpleManager.h"
//#include <string>
//#include <vector>

//!  Шаблон Фабрика (Factory). 
/*!
  Упрощенный шаблон Фабрика (Factory) построен на основе шаблона Менеджер (SimpleManager),
  предназначен для производства объектов из определенной иерархии классов с интерфейсом
  (базовым абстрактным классом) Interface. Хранит ассоциативный массив зарегистрированных 
  указателей на функции-создатели объектов иерархии, ассоциированных с некими идентификаторами,
  например строковыми названиями классов иерархии.
*/


template <typename Type, typename Interface>
class Factory : public SimpleManager<typename Type, typename Interface* (*)() >
{
	//std::vector<Interface*> created_;

	public:

	typedef Interface* (*Creator)();

	Interface* createObject(const Type& type)
	{
		if (checkData(type)) {
			Interface* ptr = getData(type)();
			//created_.push_back(ptr);
			return ptr;
		}
		return nullptr;
	}

	Factory() {}
	~Factory()
	{
		/* удалять будем самостоятельно
		for (size_t i = 0; i < created_.size(); i++) {
			if (created_[i]) delete created_[i];
		}
		*/
	}
};

