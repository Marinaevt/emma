#include "StdAfx.h"
#include "xtl/xtl_log.h"

#include "Reader.h"
#include "SpecialExpressions.h"
#include "Compare.h"

namespace MM {

	CReader::CReader() {
		m_readyUnits.clear();
		m_waitingUnits.clear();
		m_suspendedUnits.clear();
		m_readyElements.clear();
	}

	CReader::CReader(const std::string & string) {
		Reset(string);
	}

	void CReader::Reset(const std::string & string) {
		bool bTokenize = m_tokenizer.tokenize(string, m_tokens);
		// Уничтожение пробелов.
		m_tokens.erase(std::remove_if(m_tokens.begin(), m_tokens.end(), FStringIsDefinitionOf(CSpace())), m_tokens.end());

		m_currentPos = m_tokens.begin();

		m_readyUnits.clear();
		m_waitingUnits.clear();
		m_suspendedUnits.clear();

		m_readyElements.clear();
	}

	bool CReader::GetNext() {
		std::list<CUnitSptr> elements;

		m_factory.get_if(FUnitDefinitionIs(m_currentToken = *m_currentPos++), elements);

		if (!elements.empty()) {
			for (std::list<CUnitSptr>::const_iterator i = elements.begin(); i != elements.end(); ++i) {
				CUnitSptr pUnitSptr = *i;
				bool status = pUnitSptr->Load(*this);

				if (xtd::error::did_not_occure()) {
					return status;
				}
				else {
					xtd::error::reset();
				}
			}

			return xtd::error::report("Лексема \"" + m_currentToken + "\" использована неправильно.");
		}
		
		//иначе 
		return xtd::error::report("Лексемы \"" + m_currentToken + "\" нет в фабрике юнитов.");
		
	}

	bool CReader::Read() {
		while (m_currentPos != m_tokens.end() || HasSuspendedUnit()) {
			bool status;

			if (HasSuspendedUnit()) {
				status = ExtractSuspendedUnit()->Load(*this);
			}
			else {
				status = GetNext();
			}

			if (xtd::error::occured()) {
				return xtd::no_good;
			}
			// Если статус равен "xtd::no_good", но ошибки не произошло, то это означает, что загрузка выражения была
			// приостановлена, поэтому необходимо прервать чтение и вернуть вершину стека ожидающему выражению.
			else if (status == xtd::no_good) {
				break;
			}
		}

		if (HasReadyUnit()) {
			return xtd::good;
		}
		else {
			return xtd::no_good;
		}
	}

	bool CReader::Read(const std::string & string) {
		Reset(string);

		bool bRead = Read();
		return bRead;
	}

} // namespace MM