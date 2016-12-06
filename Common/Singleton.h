#pragma once

// Упрощенный шаблон [синглтон] обеспечивает
// единственность и глобальный доступ к объекту
// класса T. 

template <class T>
class CSingleton
{
public:

	static T& Instance()
	{
		if (!m_pInstance) {
			m_pInstance = new T();

			// ОБЯЗАТЕЛЬНО для вызова деструктора ~T()
			// при завершении программы
			std::atexit(Destroy);
		}

		return *m_pInstance;
	}

private:

	static void Destroy()
	{ 
		//if (m_pInstance) 
		delete m_pInstance; 
	}

	CSingleton();
	CSingleton(const CSingleton&);
	CSingleton& operator=(const CSingleton&);
	~CSingleton();

	static T* m_pInstance;
};

template<class T> T* CSingleton<T>::m_pInstance = 0;
