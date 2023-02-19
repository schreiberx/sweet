/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include <list>
#include <memory>
#include <typeinfo>
#include <sweet/SWEETError.hpp>
#include <sweet/ProgramArguments.hpp>
#include <sweet/variables/VariablesClassInterface.hpp>
class VariablesClassDictionary
{
public:
	sweet::ErrorBase error;

private:
	std::list<VariablesClassDictionaryInterface*> _list;

	bool registerationClosed;
	bool getVariableClassClosed;


public:
	void closeRegistration()
	{
		registerationClosed = true;
	}

public:
	void closeGetVariableClass()
	{
		getVariableClassClosed = true;
	}


public:
	VariablesClassDictionary()	:
		registerationClosed(false),
		getVariableClassClosed(false)
	{

	}
public:
	~VariablesClassDictionary()
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			delete *i;
		}
	}

public:
	template<typename T>
	bool registerParameterClass()
	{
		if (registerationClosed)
		{
			const std::string& tname = typeid(T).name();
			error.errorSet("Registration already closed (type '"+tname+"')");
			return false;
		}

		if (parameterClassExists<T>())
		{
			const std::string& tname = typeid(T).name();
			error.errorSet("Class of type '"+tname+"' already exists");
			return false;
		}

		T* newClass = new T();
		_list.push_back(newClass);
		return true;
	}

public:
	template<typename T>
	bool parameterClassExists()
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
	T* getVariableClassClass()
	{
		if (getVariableClassClosed)
		{
			const std::string& tname = typeid(T).name();
			error.errorSet("Getting a variable class already closed (type '"+tname+"')");
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
		error.errorSet("Type '"+tname+"' not found in dictionary");
		return nullptr;
	}

public:
	void outputProgramArguments(std::string i_prefix = "")
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			(*i)->outputProgramArguments(i_prefix);
		}
	}

public:
	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			if (!(*i)->processProgramArguments(i_pa))
			{
				error.errorForward((*i)->getError());
				return false;
			}
		}
		return true;
	}

public:
	void outputVariables(std::string i_prefix = "")
	{
		for (auto i = _list.begin(); i != _list.end(); i++)
		{
			(*i)->outputVariables(i_prefix);
		}
	}

};



#include <iostream>
#include "swe_sphere_variables/PDESWESphereParameters.hpp"
#include "swe_sphere_variables/IODataParameters.hpp"


int
main(int i_argc, char *i_argv[])
{
	/*
	 * We start by setting up the class to parse the program arguments
	 */
	std::cout << " + ProgramArguments()" << std::endl;
	sweet::ProgramArguments pa;
	if (!pa.setup(i_argc, i_argv))
	{
		std::cout << "Error: " << pa.error.errorGet() << std::endl;
		return 1;
	}


	{
		/*
		 * Now we warmup the dictionary which can be stored to store arbitrary
		 * classes.
		 *
		 * This is specialized for processing our program arguments and using it
		 * later on in the SWEET programs.
		 */
		std::cout << " + VariablesClassDictionary()" << std::endl;
		VariablesClassDictionary varClassDict;


		/*
		 * Register new classes
		 */
		std::cout << "   + registerParameterClass<PDESWEParametersSphere>()" << std::endl;
		varClassDict.registerParameterClass<PDESWEParametersSphere>();
		varClassDict.registerParameterClass<IODataParameters>();

		/*
		 * Now we close the registration
		 *
		 * This will avoid performance bugs!
		 */
		varClassDict.closeRegistration();


		/*
		 * After registering all classes, we can check whether we should output the help information
		 */
		if (pa.argumentWithKeyExists("h") || pa.argumentWithKeyExists("help"))
		{
			varClassDict.outputProgramArguments();
			return EXIT_FAILURE;
		}

		/*
		 * Now its time to process all program arguments with all registered classes
		 */
		varClassDict.processProgramArguments(pa);

		/*
		 * Get handler to new class PDESWEParametersSphere
		 */
		PDESWEParametersSphere *sweParametersSphere = varClassDict.getVariableClassClass<PDESWEParametersSphere>();
		if (sweParametersSphere == nullptr)
		{
			std::cerr << "Not a SWEET error: " << varClassDict.error.errorGet() << std::endl;
			return EXIT_FAILURE;
		}

		/*
		 * Get handler to new class PDESWEParametersSphere
		 */
		IODataParameters *ioDataParameters = varClassDict.getVariableClassClass<IODataParameters>();
		if (ioDataParameters == nullptr)
		{
			std::cerr << "Not a SWEET error: " << varClassDict.error.errorGet() << std::endl;
			return EXIT_FAILURE;
		}

		/*
		 * Now we close getting the parameter class
		 *
		 * This will avoid performance bugs!
		 */
		varClassDict.closeGetVariableClass();


		/*
		 * Now its time to process all program arguments with all registered classes
		 */
		std::cout << " + varClassDict.outputVariables()" << std::endl;
		varClassDict.outputVariables("    ");


		/*
		 * And we can also print them individually
		 */
		std::cout << " + sweParametersSphere->outputVariables()" << std::endl;
		sweParametersSphere->outputVariables("    ");

		/*
		 * And we can also print them individually
		 */
		std::cout << " + ioDataParameters->outputVariables()" << std::endl;
		ioDataParameters->outputVariables("    ");
	}

	return 0;
}

