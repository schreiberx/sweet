/*
 * Parareal_CoutPrefix.hpp
 *
 *  Created on: 18 Apr 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_CONSOLEPREFIX_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_CONSOLEPREFIX_HPP_


#include <iostream>

/**
 * Prefix the std::cout output with a given string
 */
class Parareal_ConsolePrefix
{
#if 1
private:
	/*
	 * Prefix std::cout with string
	 *
	 * Source: http://stackoverflow.com/questions/27336335/c-cout-with-prefix
	 */
	class PrefixBuffer
	    : public std::streambuf
	{
	    std::string     prefix;
	    std::streambuf* sbuf;
	    bool            need_prefix;

	    int sync()
	    {
	        return this->sbuf->pubsync();
	    }

	    int overflow(int c)
	    {
	        if (c != std::char_traits<char>::eof()) {
	            if (this->need_prefix
	                && !this->prefix.empty()
	                && (int)this->prefix.size() != this->sbuf->sputn(&this->prefix[0], this->prefix.size())) {
	                return std::char_traits<char>::eof();
	            }
	            this->need_prefix = c == '\n';
	        }
	        return this->sbuf->sputc(c);
	    }

	public:
	    PrefixBuffer()	:
	        sbuf(0),
	        need_prefix(true)
		{
	    }

	public:
	    void setPrefix(const std::string& i_prefix)
	    {
	    	prefix = i_prefix;
	    }

	public:
	    void setPrefix(const char *i_prefix)
	    {
	    	prefix = i_prefix;
	    }


	public:
	    void setStreamBuf(std::streambuf* i_sbuf)
	    {
	    	sbuf = i_sbuf;
	    }
	};

	PrefixBuffer coutPrefixBuffer, cerrPrefixBuffer;
#endif

	std::streambuf *coutbuf, *cerrbuf;

public:
	void start(
			int i_number		///< prefix number will be extended to "[NUMBER] "
	)
	{
		std::ostringstream ss;
		ss << "[" << i_number << "] ";

		coutPrefixBuffer.setPrefix(ss.str());
		cerrPrefixBuffer.setPrefix(ss.str());

		std::cout.rdbuf(&coutPrefixBuffer);
		std::cerr.rdbuf(&cerrPrefixBuffer);
	}


	void start(
			const char *i_prefix	///< prefix string
	)
	{
		coutPrefixBuffer.setPrefix(i_prefix);
		std::cout.rdbuf(&coutPrefixBuffer);

		cerrPrefixBuffer.setPrefix(i_prefix);
		std::cerr.rdbuf(&cerrPrefixBuffer);
	}



	void end()
	{
		std::cout.rdbuf(coutbuf);
		std::cerr.rdbuf(cerrbuf);
	}


	Parareal_ConsolePrefix()
	{
		coutbuf = std::cout.rdbuf();
		cerrbuf = std::cerr.rdbuf();

		coutPrefixBuffer.setStreamBuf(coutbuf);
		cerrPrefixBuffer.setStreamBuf(cerrbuf);
	}


	~Parareal_ConsolePrefix()
	{
		end();
	}

};




#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_CONSOLEPREFIX_HPP_ */
