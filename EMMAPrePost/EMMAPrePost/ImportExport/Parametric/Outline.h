#pragma once

#include <iostream>
#include <list>

#include "Unit.h"
#include "Parameter.h"
#include "Compare.h"
#include "Contour.h"
#include "../../../../Common/2DOutline.h"
#include "../../../../Common/Startup.h"

/*!
	Пространство имён параметрической модели.
	Сокращение от слова метамодель (MetaModel).
*/
namespace MM {

class COutline {
public:
	// ---- Определения ---------------------------------------------------------------------------
	
	typedef std::pair<std::string, double> ValueType;

	// ---- Модификация ---------------------------------------------------------------------------
	
	/*!
			Установка новых графических элементов чертежа.
			Копирует указатели на новые элементы в массив графических элементов "m_elements", уничтожая старые элементы.
			Указатели на элементы-параметры дублируются в массив "m_parameters" для удобства дальнейшней работы с ними.
	*/
	template <typename InputContainer>
	void SetElements (const InputContainer & elements) {
		ElementsContainer newElements, newParameters;

		for (typename InputContainer::const_iterator i = elements.begin(); i != elements.end(); ++i) {
			CUnitSptr pUintSPtr = *i;
			newElements.push_back(pUintSPtr);
			// Копирование указателей на параметры в отдельный массив.
			if (pUintSPtr->DefinitionIs(CParameter().Definition())) {
				newParameters.push_back(pUintSPtr);
			}
		}
			
		std::swap(newElements, m_elements);
		std::swap(newParameters, m_parameters);
	}
	
	//!		Установка новых значений параметров.
	/*!
			Принимает массив пар «имя-значение», где имя и значение соответствуют одному из параметров
		модели, ищет параметры по имени и присваивает им новые значения.
	*/
	template <typename InputContainer>
	void SetParametersValues (InputContainer & container) {
		for (ElementsContainer::const_iterator i = m_parameters.begin(); i != m_parameters.end(); ++i) {
			CParameter & parameter = dynamic_cast<CParameter &>(**i);
			
			for (typename InputContainer::const_iterator j = container.begin(); j != container.end(); ++j) {
				if (parameter.NameIs(j->first)) {
					parameter.SetValue(j->second);
					break;
				}
			}
		}
	}
	
	//!		Установка нового значения параметра.
	/*!
			Присваивает параметру по имени "name" значение "value".
	*/
	void SetParameterValue (const std::string & name, double value);

	// ---- Доступ --------------------------------------------------------------------------------
	
	//!		Взятие имён и значений параметров модели.
	/*!
			Записывает в массив "container" пары «имя-значение», соответствующие параметрам модели.
	*/
	template <typename OutputContainer>
	void GetParameters (OutputContainer & container) const {
		OutputContainer parameters;

		for (ElementsContainer::const_iterator i = m_parameters.begin(); i != m_parameters.end(); ++i) {
			const CParameter & parameter = dynamic_cast<const CParameter &>(**i);
			
			parameters.push_back(ValueType(parameter.Name(), parameter.DoubleValue()));
		}
		
		std::swap(container, parameters);
	}
	
	//!		Конкретизация метамодели.
	/*!
			Создание геометрии метамодели на основе заданных на данный момент значений параметров
		и экспорт получившейся геометрии в класс "C2DOutline".
	*/
	C2DOutline * Instantiate () const;
	
	// ---- Ввод-вывод ----------------------------------------------------------------------------
	
	bool Load (std::istream & stream);
	bool LoadFile (const CString &pathToFile);
	
	void Save (std::ostream & stream) const;
	void SaveFile(const char * pathToFile) const;
    
private:
    void InitReader (CReader & reader) const;
    void InitReaderExpressions (CReader & reader) const;
    void InitReaderElements (CReader & reader) const;

private:
	typedef std::list<CUnitSptr> ElementsContainer;
	
	ElementsContainer m_elements;	//массив ВСЕХ элементов (т.е. всего)
	ElementsContainer m_parameters; //! В этом массиве хранятся только параметры. Для удобства доступа к ним.
};

} // namespace MM