#include "StdAfx.h"
#include <fstream>

#include "xtl/xtl_log.h"

#include "Outline.h"
#include "Reader.h"

#include "Operands.h"
#include "Operators.h"
#include "Functions.h"
#include "SpecialExpressions.h"

#include "Parameter.h"
#include "Point.h"
#include "Curves.h"
#include "Contour.h"
#include "SpecialElements.h"
#include "Compare.h"

namespace MM {

	bool COutline::Load(std::istream & stream) {
		std::string string;
		CReader reader;

		InitReader(reader);

		//получаем весь файл в строку, которую будем разбивать на лексемы
		getline(stream, string, '\0');

		if (!string.empty()) {
			bool bRead = reader.Read(string);
			if (reader.HasReadyElement()) {
				const xtd::stack<CUnitSptr>::container ReaderContainer = reader.GetReadyElements();
				SetElements(ReaderContainer);
			} else {
				bool bErrRep = xtd::error::report("Не удалось загрузить файл.");	//+ устанавливает status ошибки
				xtd::error::reset();	//сбрасываем сохранённую ошибку для другой загрузки
				return bErrRep;	//true - не загрузили
			}
		}
		else {
			bool bErrRep = xtd::error::report("Файл пуст.");	//+ устанавливает status ошибки
			xtd::error::reset();	//сбрасываем сохранённую ошибку для другой загрузки
			return bErrRep;	//true - файл пуст
		}

		bool bErrRep = xtd::success::report();
		return bErrRep;	//false - успешно загрузили
	}

	//Загружаем чертёж из параметрического файла .meta
	bool COutline::LoadFile(const CString &pathToFile) {
		// конвертируем строку в ANSI
		CT2CA pszConvertedAnsiString(pathToFile);

		std::ifstream stream(pszConvertedAnsiString);
		
		bool bResLoad = false;		//result of Load (success or error)
		if (stream.is_open()) {
			
			bResLoad = Load(stream);
			return bResLoad;
		}
		
		// иначе ошибка - файла нет
		bResLoad = xtd::error::report("Файла \"" + std::string(pszConvertedAnsiString) + "\" не существует.");
		return bResLoad;

	}

	void COutline::Save(std::ostream & stream) const {
		for (ElementsContainer::const_iterator i = m_elements.begin(); i != m_elements.end(); ++i) {
			(*i)->Save(stream);
			stream << ";\n";
		}
	}

	void COutline::SaveFile(const char * pathToFile) const {
		std::ofstream stream(pathToFile);

		if (stream.is_open()) {
			Save(stream);
		}
		else {
			xtd::error::report("Не удалось создать файл \"" + std::string(pathToFile) + "\".");
		}
	}

	void COutline::SetParameterValue(const std::string & name, double value) {
		for (ElementsContainer::iterator i = m_parameters.begin(); i != m_parameters.end(); ++i) {
			CParameter & parameter = dynamic_cast<CParameter &>(**i);

			if (parameter.NameIs(name)) {
				parameter.SetValue(value);
				break;
			}
		}
	}

	void COutline::InitReader(CReader & reader) const {
		InitReaderExpressions(reader);
		InitReaderElements(reader);
	}

	void COutline::InitReaderElements(CReader & reader) const {
		reader.Put<CParameter>();
		reader.Put<CPoint>();

		reader.Put<CLine>();
		reader.Put<CCircle>();
		reader.Put<CArc>();

		reader.Put<CContour>();

		reader.Put<CVariable>();
		reader.Put<CVariableLoader>();

		reader.Put<CElementEnd>();
	}

	void COutline::InitReaderExpressions(CReader & reader) const {
		reader.Put<CNumber>();
		reader.Put<CPi>();
		reader.Put<CExponent>();

		reader.Put<CNegation>();
		reader.Put<CUnaryPlus>();

		reader.Put<CSum>();
		reader.Put<CDifference>();
		reader.Put<CProduct>();
		reader.Put<CQuotient>();
		reader.Put<CPower>();

		//?	reader.Put<CFactorial>();

		reader.Put<CSine>();
		reader.Put<CCosine>();
		reader.Put<CLogarithm>();
		// reader.Put<CSquareRoot>();
		// reader.Put<CSymmetryX>();

		reader.Put<COpeningBracket>();
		reader.Put<CClosingBracket>();
		reader.Put<CComma>();
		//	Put<CComment>();
		reader.Put<CCommentStart>();
		reader.Put<CCommentEnd>();
		reader.Put<CSpace>();
	}

	C2DOutline * COutline::Instantiate() const {

		C2DOutline * outline = dynamic_cast<C2DOutline *>(IO::CreateObject(C2DOUTLINE));

		int nExportResult = -1;	//результат экспорта
		for (ElementsContainer::const_iterator i = m_elements.begin(); i != m_elements.end(); ++i) {
			if ((*i)->DefinitionIs(CContour().Definition())) {
				nExportResult = static_cast<const CElement &>(**i).ExportToOutline(outline);
			}
		}

		//если где-то возникла ошибка, удаляем, возвращаем нуль
		if (nExportResult < 0) {
			delete outline;
			return nullptr;
		}

		return outline;
	}

} // namespace MM