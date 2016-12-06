#pragma once

#include <map>
#include <vector>

//!  Шаблон Менеджер (SimpleManager). 
/*!
  Упрощенный шаблон Менеджер (SimpleManager) - по сути некая оболочка для ассоциативного
  массива std::map, предназначенная для удобства.
*/


template <typename Id, typename Data>
class SimpleManager
{
public:

	typedef std::map<Id, Data> DataMap;
	typedef std::vector<Id> IdVector;

protected:

	DataMap map_;

public:

	bool registerData(const Id& id, Data data)
	{
		return map_.insert(DataMap::value_type(id, data)).second;
	}

	bool unregisterData(const Id& id)
	{
		return map_.erase(id) == 1;
	}

	bool checkData(const Id& id) const
	{
		return (map_.count(id) == 1);
	}

	Data getData(const Id& id) const
	{
		DataMap::const_iterator i = map_.find(id);

		if (i == map_.end())
		{
			// по хорошему нужно вызывать исключение
			//
			return Data();
		}

		return i->second;
	}

	IdVector getRegisteredIds() const
	{
		DataMap::const_iterator i = map_.begin();
		IdVector temp;

		while (i != map_.end())
		{
			temp.push_back(i->first);
			++i;
		}

		return temp;
	}
};
