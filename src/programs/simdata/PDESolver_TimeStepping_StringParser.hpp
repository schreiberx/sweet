#ifndef SRC_INCLUDE_SWEET_PDESOLVER_AMARGUMENT_PARSER_HPP_
#define SRC_INCLUDE_SWEET_PDESOLVER_AMARGUMENT_PARSER_HPP_

#include <string>
#include <vector>
#include <ostream>
#include <memory>

#include "PDESolver_TimeStepping_Tree.hpp"
#include <sweet/core/ErrorBase.hpp>

namespace sweet
{
/*
 * Parses strings of the form
 *
 * function_name_or_NAME([arguments])
 *
 * where [arguments] should be parsed by ParserArguments
 */
class PDESolver_TimeSteppingStringParser
{
public:
	ErrorBase error;


	/**
	 * A string parser with simple operations for parsing
	 */
	class StringParser
	{
	public:
		// The time stepping string itself handed over during the setup
		std::string _timestepping_string;

		// reading position of the time stepping string
		int _reading_position;

		// maximum reading position
		int _max_reading_position;

	public:
		StringParser()
		{
		}

	public:
		void setup(const std::string &i_timestepping_string)
		{
			_timestepping_string = i_timestepping_string;
			_reading_position = 0;
			_max_reading_position = _timestepping_string.length();
		}

		bool isAllParsed()
		{
			return _reading_position == (int)_timestepping_string.length();
		}

		/**
		 * Skip whitespaces
		 */
	public:
		void parse_skipWhitespace()
		{
			std::string &str = _timestepping_string;
			int &pos = _reading_position;

			while (true)
			{
				if (pos >= _max_reading_position)
					break;

				char c = str[pos];

				if (c != ' ')
					break;

				pos++;
			}
		}

		/**
		 * Parse next identifier
		 *
		 * [a-Z0-9_]*
		 */
	public:
		bool parse_nextIdentifier(
				std::string &o_token
		)
		{
			parse_skipWhitespace();

			std::string &str = _timestepping_string;
			int &pos = _reading_position;

			// accummulator for token
			std::string acc;

			while (true)
			{
				char c = str[pos];

				// Check for string
				if (	(c >= 'a' && c <= 'z') ||
						(c >= 'A' && c <= 'Z') ||
						(c >= '0' && c <= '9') ||
						(c == '_')
					)
				{
					// add to accummulator
					acc += c;
					pos++;
					continue;
				}

				// Character not determined => Stop
				break;
			}

			o_token = acc;
			return o_token.length() != 0;
		}

		/**
		 * Parse for open parenthesis '
		 */
	public:
		bool parse_particularCharacter(char i_character)
		{
			parse_skipWhitespace();

			std::string &str = _timestepping_string;
			int &pos = _reading_position;

			if (pos >= _max_reading_position)
				return false;

			if (str[pos] != i_character)
				return false;

			pos++;

			return true;
		}

		/**
		 * Return some error information string which can be printed to the console
		 * for further debugging information where the error occurred.
		 */
		std::string getErrorInfo(
				const std::string &i_message = "Error location"
		)
		{
			std::ostringstream oss;

			oss << _timestepping_string << std::endl;

			for (int i = 0; i < _reading_position; i++)
				oss << " ";
			oss << "^~ " << i_message << std::endl;

			return oss.str();
		}
	};

	/*
	 * Handler to string parser
	 */
	StringParser stringParser;

	/*
	 * The main function which was parsed.
	 * All time stepping strings need to have one main function
	 */
public:
	bool _setupFinished;
	bool _errorForDoubleParsing;
	bool _errorForDuplicateKeysInParsing;


public:
	PDESolver_TimeSteppingStringParser(
			bool i_errorForDoubleParsing = true,	///< Trigger an error if an argument is parsed twice
			bool i_errorForDuplicateKeysInParsing = true		///< Trigger an error if there are duplicate keys
	)	:
		_setupFinished(false),
		_errorForDoubleParsing(i_errorForDoubleParsing),
		_errorForDuplicateKeysInParsing(i_errorForDuplicateKeysInParsing)
	{
	}


	void clear()
	{
		_setupFinished = false;
	}

	bool genTimeSteppingTree(
			const std::string &i_timeSteppingString,
			PDESolver_TimeStepping_Tree &o_timeSteppingTree
	)
	{
		std::string tmp;

		stringParser.setup(i_timeSteppingString);

		if (!stringParser.parse_nextIdentifier(tmp))
			return error.set("Failed to find first token in string - is it maybe empty?\n"+stringParser.getErrorInfo());

		std::string debug_message = stringParser.getErrorInfo("Location or error");

		o_timeSteppingTree.clear();
		o_timeSteppingTree.mainFunction = std::make_shared<PDESolver_TimeStepping_Tree::Function>(tmp);

		// store potential debugging information in time stepping string to use it later
		o_timeSteppingTree.mainFunction->setDebugMessage(debug_message);

		if (!stringParser.parse_particularCharacter('('))
			return error.set("Open parenthesis '(' expected!\n"+stringParser.getErrorInfo());

		if (!_setup_functionArguments(o_timeSteppingTree.mainFunction))
			return false;

		if (!stringParser.parse_particularCharacter(')'))
			return error.set("Closing parenthesis ')' expected!\n"+stringParser.getErrorInfo());

		stringParser.parse_skipWhitespace();

		if (!stringParser.isAllParsed())
			return error.set("Parsing failed at this position!\n"+stringParser.getErrorInfo());
		return true;
	}

	/*
	 * This setup routine is called once a function has been detected.
	 *
	 * This includes that also the opening parenthesis is parsed.
	 * The closing one will be parsed from the caller
	 */
	bool _setup_functionArguments(
			std::shared_ptr<PDESolver_TimeStepping_Tree::Function> io_function
	)
	{
		while (true)
		{
			std::string debug_message = stringParser.getErrorInfo("Location or error");

			// Search for first identifier
			// if it doesn't exist, assume that there are no further arguments
			std::string identifier;
			if (!stringParser.parse_nextIdentifier(identifier))
				return true;

			/*
			 * Determine kind of argument:
			 *
			 * key-value, key-function
			 *   or
			 * value-only
			 *   or
			 * function
			 */

			// Prepare argument
			std::shared_ptr<PDESolver_TimeStepping_Tree::Argument> arg = std::make_shared<PDESolver_TimeStepping_Tree::Argument>();
			io_function->arguments.push_back(arg);

			// store potential debugging information in time stepping string to use it later
			arg->setDebugMessage(debug_message);

			/*
			 * Check for key-value or key-function
			 */
			if (stringParser.parse_particularCharacter('='))
			{
				std::string &key = identifier;

				std::string value;
				if (!stringParser.parse_nextIdentifier(value))
					return error.set("Identifier expected!\n"+stringParser.getErrorInfo());

				// Determine function
				if (stringParser.parse_particularCharacter('('))
				{
					// key-function
					arg->argType = PDESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION;
					arg->key = key;
					arg->function = std::make_shared<PDESolver_TimeStepping_Tree::Function>(value);

					// store potential debugging information in time stepping string to use it later
					arg->function->setDebugMessage(debug_message);

					if (!_setup_functionArguments(arg->function))
						return false;

					if (!stringParser.parse_particularCharacter(')'))
						return error.set("Closing parenthesis ')' expected!\n"+stringParser.getErrorInfo());
				}
				else
				{
					// key-value
					arg->argType = PDESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE;
					arg->key = key;
					arg->value = value;
				}
			}
			else if (stringParser.parse_particularCharacter('('))
			{
				// key-function
				arg->argType = PDESolver_TimeStepping_Tree::Argument::ARG_TYPE_FUNCTION;
				arg->function = std::make_shared<PDESolver_TimeStepping_Tree::Function>(identifier);

				// store potential debugging information in time stepping string to use it later
				arg->function->setDebugMessage(debug_message);

				if (!_setup_functionArguments(arg->function))
					return false;

				if (!stringParser.parse_particularCharacter(')'))
					return error.set("Closing parenthesis ')' expected!\n"+stringParser.getErrorInfo());
			}
			else
			{
				// No parenthesis => value

				// value
				arg->argType = PDESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE;
				arg->value = identifier;
			}

			// If ',' is not found, this was the last argument
			if (!stringParser.parse_particularCharacter(','))
				break;
		}

		return true;
	}

#if 0
	friend
	std::ostream&
	operator<<(std::ostream &io_os, const PDESolver_TimeSteppingStringParser &i_pa)
	{
		for (std::size_t i = 0; i < i_pa._arguments.size(); i++)
		{
			const Argument &a = i_pa._arguments[i];
			std::cout << "'" << a.key << "' => '" << a.value << "' (processed=" << a.argumentParsedAndAccessed << ")" << std::endl;
		}
		return io_os;
	}
#endif
};

}

#endif
