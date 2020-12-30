#ifndef SWEET_FILE_OPERATIONS_HPP__
#define SWEET_FILE_OPERATIONS_HPP__

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


class FileOperations
{
public:
	static
	bool file_exists(const std::string& i_filename)
	{
		struct stat buffer;
		return (stat(i_filename.c_str(), &buffer) == 0);
	}

	static
	bool file_is_readable(const std::string& i_filename)
	{
		if (FILE *file = fopen(i_filename.c_str(), "r"))
		{
			fclose(file);
			return true;
		}
		else
		{
			return false;
		}
	}
};


#endif
