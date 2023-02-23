/*
 *  Created on: Feb 19, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKDICTIONARY_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKDICTIONARY_HPP_


#include <list>
#include <memory>
#include <typeinfo>
#include <sweet/ErrorBase.hpp>
#include <sweet/ProgramArguments.hpp>
#include <sweet/shacks/ShackInterface.hpp>

namespace sweet
{

/*
 * A dictionary using class types as key.
 *
 * It also integrates parsing of program arguments.
 */
class ShackDictionary
{
public:
	sweet::ErrorBase error;

private:
	std::list<ShackInterface*> _list;

	bool _registerationOfClassInstanceFinished;
	bool _getFinished;


public:
	void registrationFinished()
	{
		_registerationOfClassInstanceFinished = true;
	}


public:
	void getFinished()
	{
		_getFinished = true;
	}


public:
	ShackDictionary()
	{
		clear();
	}


public:
	void clear()
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			delete *i;
		}
		_list.clear();

		error.reset();

		_registerationOfClassInstanceFinished = false;
		_getFinished = false;
	}


public:
	~ShackDictionary()
	{
		clear();
	}


public:
	template<typename T>
	bool registerFirstTime()
	{
		if (_registerationOfClassInstanceFinished)
		{
			const std::string& tname = typeid(T).name();
			error.set("Registration already finished (type '"+tname+"')");
			return false;
		}

		if (exists<T>())
		{
			const std::string& tname = typeid(T).name();
			error.set("Class of type '"+tname+"' already exists");
			return false;
		}

		T* newClass = new T();
		_list.push_back(newClass);
		return true;
	}

public:
	template<typename T>
	bool exists()
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			// check whether generic interface can be casted to type T
			T* derived = dynamic_cast<T*>(*i);

			if (derived != nullptr)
				return true;
		}

		return false;
	}

public:
	template<typename T>
	T* get(bool i_auto_registration = false)
	{
		if (!i_auto_registration)
		{
			if (!_registerationOfClassInstanceFinished)
			{
				const std::string& tname = typeid(T).name();
				error.set("Registration of class instances needs to be finished first (type '"+tname+"')");
				return nullptr;
			}
		}

		if (_getFinished)
		{
			const std::string& tname = typeid(T).name();
			error.set("Getting an element class already finished (type '"+tname+"')");
			return nullptr;
		}

		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			// check whether generic interface can be casted to type T
			T* derived = dynamic_cast<T*>(*i);

			if (derived != nullptr)
				return derived;
		}

		const std::string& tname = typeid(T).name();
		error.set("Type '"+tname+"' not found in dictionary");
		return nullptr;
	}

	/*
	 * Auto registrate this particular class if it doesn't exist and return an instance
	 */
public:
	template<typename T>
	T* getAutoRegistration()
	{
		if (!exists<T>())
			registerFirstTime<T>();

		return get<T>(true);
	}


public:
	void printProgramArguments(const std::string& i_prefix = "")
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			(*i)->printProgramArguments(i_prefix);
		}
	}

public:
	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			if (!((*i)->processProgramArguments(i_pa)))
			{
				error.forward((*i)->error);
				return false;
			}
		}
		return true;
	}

public:
	void printShackData(const std::string& i_prefix = "")
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			(*i)->printShack(i_prefix);
		}
	}

};

}

#endif
