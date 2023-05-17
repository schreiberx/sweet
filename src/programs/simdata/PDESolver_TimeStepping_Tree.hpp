#ifndef SRC_PROGRAMS_SIMDATA_PDESOLVER_TIMESTEPPING_TREE_HPP_
#define SRC_PROGRAMS_SIMDATA_PDESOLVER_TIMESTEPPING_TREE_HPP_


namespace sweet
{

class PDESolver_TimeStepping_Tree
{
public:
	PDESolver_TimeStepping_Tree()
	{
	}

	~PDESolver_TimeStepping_Tree()
	{
		clear();
	}


public:
	class Argument;
	class Function;


public:
	std::shared_ptr<PDESolver_TimeStepping_Tree::Function> mainFunction;


	void clear()
	{
		mainFunction.reset();
	}

	/*
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
		std::string debug_message;

		bool argumentParsedAndAccessed;


		Argument()	:
			argType(ARG_INVALID),
			argumentParsedAndAccessed(false)
		{
		}

		void reset()
		{
			key = "";
			value = "";
			argumentParsedAndAccessed = false;
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
