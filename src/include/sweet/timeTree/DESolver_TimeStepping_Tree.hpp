#ifndef SRC_PROGRAMS_SIMDATA_PDESOLVER_TIMESTEPPING_TREE_HPP_
#define SRC_PROGRAMS_SIMDATA_PDESOLVER_TIMESTEPPING_TREE_HPP_


namespace sweet
{

class DESolver_TimeStepping_Tree
{
public:
	DESolver_TimeStepping_Tree()
	{
	}

	~DESolver_TimeStepping_Tree()
	{
		clear();
	}


public:
	class Argument;
	class Function;


public:
	std::shared_ptr<DESolver_TimeStepping_Tree::Function> mainFunction;


	void clear()
	{
		mainFunction.reset();
	}

	/**
	 * This represents an entire function of the form
	 *
	 * "FunctionNameAnd123SomeNumbers(these, are=parameters, and, can, also, be, a, recursive, function(yadda))"
	 */
public:
	class Function
	{
	public:
		/// Function name itself
		std::string function_name;

	public:
		/// Arguments of function
		std::vector<std::shared_ptr<Argument>> arguments;

	private:
		/// String which contains information about the function to debug
		/// (in case of an error in parsing the arguments in the time stepper)
		std::string _debugMessage;

	public:
		Function(
				const std::string &i_function_name
		)	:
			function_name(i_function_name)
		{
		}

		~Function()
		{
			clear();
		}

	public:
		void clear()
		{
			function_name = "";
			arguments.clear();
		}

	public:
		/*
		 * Output information about this function and its arguments including
		 * recursive calls to other functions
		 */
		void print(const std::string &i_prefix_str = "")
		{
			std::cout << i_prefix_str << function_name << "(" << std::endl;
			std::string new_prefix_str = i_prefix_str + "  ";

			for (std::size_t i = 0; i < arguments.size(); i++)
				arguments[i]->print(new_prefix_str);


			std::cout << i_prefix_str << ")" << std::endl;
		}

public:
		void setDebugMessage(const std::string &i_debugMessage)
		{
			_debugMessage = i_debugMessage;
		}
		std::string getDebugMessage()
		{
			return _debugMessage;
		}
		std::string getNewLineDebugMessage()
		{
			return "\n"+_debugMessage;
		}
	};

public:
	/**
	 * Class which represents one particular argument of a function
	 */
	class Argument
	{
	public:
		ErrorBase error;

		enum ArgumentType
		{
			ARG_INVALID,
			ARG_TYPE_FUNCTION,	// just a function
			ARG_TYPE_KEY_VALUE,	// a key and a value (2 strings)
			ARG_TYPE_KEY_FUNCTION,	// a key and a function as a value
			ARG_TYPE_VALUE,	// only a value
		};
		ArgumentType argType;

		std::shared_ptr<Function> function;
		std::string key;
		std::string value;

		/*
		 * Argument debugging message which can be printed in case there's something
		 * wrong with this argument
		 */
private:
		std::string _debugMessage;

		bool _argumentParsedAndAccessed;


public:
		Argument()	:
			argType(ARG_INVALID),
			_argumentParsedAndAccessed(false)
		{
		}

		void reset()
		{
			key = "";
			value = "";
			_argumentParsedAndAccessed = false;
		}

		bool getValue(bool &o_value)
		{
			if (value == "true")
			{
				o_value = true;
				return true;
			}

			if (value == "false")
			{
				o_value = false;
				return true;
			}

			/*
			 * Interpret 0 as false and everything else as true
			 */
			int value_int;
			try
			{
				value_int = std::stoi(value);
			}
			catch (const std::exception &e)
			{
				return error.set("Exception caught during conversion of value '"+value+"' to integer: "+e.what());
			}

			o_value = (value_int != 0);
			return true;
		}

		bool getValue(int &o_value)
		{
			try
			{
				o_value = std::stoi(value);
			}
			catch (const std::exception &e)
			{
				error.set("Exception caught during conversion of value '"+value+"' to integer: "+e.what());
				return false;
			}

			return true;
		}

		bool getValue(double &o_value)
		{
			try
			{
				o_value = std::stod(value);
			}
			catch (const std::exception &e)
			{
				error.set("Exception caught during conversion of value '"+value+"' to integer: "+e.what());
				return false;
			}

			return true;
		}

	public:
		void print(const std::string &i_prefix_str = "")
		{
			std::string new_prefix_str = i_prefix_str + "  ";

			switch(argType)
			{
			case ARG_TYPE_FUNCTION:
				function->print(new_prefix_str);
				break;

			case ARG_TYPE_KEY_VALUE:
				std::cout << i_prefix_str << "'" << key << "' => '" << value << "'" << std::endl;
				break;

			case ARG_TYPE_KEY_FUNCTION:
				std::cout << i_prefix_str << "'" << key << "' => FUNCTION" << std::endl;
				function->print(new_prefix_str);
				break;

			case ARG_TYPE_VALUE:
				std::cout << i_prefix_str << "'" << value << "'" << std::endl;
				break;

			case ARG_INVALID:
				break;
			}
		}

public:
		void setDebugMessage(const std::string &i_debugMessage)
		{
			_debugMessage = i_debugMessage;
		}
		std::string getDebugMessage()
		{
			return _debugMessage;
		}
		std::string getNewLineDebugMessage()
		{
			return "\n"+_debugMessage;
		}
	};


	/**
	 * Proxy printing
	 */
	void print(const std::string i_prefix_str = "")
	{
		mainFunction->print(i_prefix_str);
	}

};

}

#endif
