/*
 * Copyright 2010 Martin Schreiber
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


/*
 * CGlRenderOStream.hpp
 *
 *  Created on: Mar 22, 2010
 *      Author: martin
 */

#ifndef CGLRENDEROSTREAM_HPP_
#define CGLRENDEROSTREAM_HPP_

#include <sstream>
#include <ostream>
#include "libgl/hud/CGlFreeType.hpp"

/**
 * this class offers the usual operator<< overloaded string output handling
 */
class CGlRenderOStream	:
	public std::ostream
{
	/**
	 * abstraction for stringbuf to catch sync() function for drawing
	 */
	class CGlRenderStreamBuf	:
				public std::stringbuf
	{
		CGlFreeType &free_type;

public:
		CGlRenderStreamBuf(	CGlFreeType &p_free_type)	:
				free_type(p_free_type)
		{
		}


        virtual int sync()
        {
			free_type.renderString(str().c_str());
			str("");
			return 0;
        }
	};


	CGlFreeType &free_type;
	CGlRenderStreamBuf streambuf;

public:
	CGlRenderOStream(CGlFreeType &p_free_type)	:
			std::ostream(&streambuf),
			free_type(p_free_type),
			streambuf(free_type)
	{

	}
};

#endif /* CGLRENDEROSTREAM_HPP_ */
