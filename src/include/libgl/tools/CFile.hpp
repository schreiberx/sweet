#ifndef FILE_H
#define FILE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <errno.h>


/**
 * \brief read the contents of a file or store contents to a file
 */
class FileContents
{
public:

	/**
	 * read the contents to string
	 */
	static bool getTextFileContent(
			const std::string &filePath,	///< filePath to source file
			std::string &p_data,			///< data where to store the filecontents
			std::string &errorLog			///< error output
	)	{
		int length;
		char *data;

		FILE *file = fopen(filePath.c_str(), "rb");
		if (file == nullptr)
		{
			std::cerr << "ERROR: getTextFileContent(" << filePath << ") fopen: " << strerror(errno) << std::endl;
			return false;
		}

		fseek(file, 0, SEEK_END);
		length = ftell(file);
		fseek(file, 0, SEEK_SET);

		data = new char[length];

		if (!data)
		{
			std::cerr << "ERROR: getTextFileContent(" << filePath << ") out of memory!" << std::endl;
			return false;
		}

		int readLength;
		readLength = fread(data, 1, length, file);
		fclose(file);

		if (readLength != length)
		{
			std::cerr << "ERROR: getTextFileContent(): readLength != fileLength" << std::endl;
			return false;
		}

		p_data.assign(data, readLength);

		delete[] data;

		return true;
	}


	/**
	 * store contents to a file
	 */
	static
	bool storeContentToFile(	const char *filePath,	///< path to file to store content to
						void *data,				///< pointer to dataset
						size_t length			///< length of dataset in bytes
	)
	{
		FILE *file = fopen(filePath, "w");
		if (file == nullptr)
		{
			std::cerr << "fileContent(" << filePath << ") fopen: " << strerror(errno) << std::endl;
			return 0;
		}

		size_t x = fwrite(data, length, 1, file);
		fclose(file);

		// error check
		if (x < length)
			return false;

		return true;
	}


	/**
	 * load data from file
	 * \param filePath	file to load
	 * \param data		storage area
	 * \param max_length	maximum length
	 * \return		-1: error	else: loaded bytes
	 */
	static
	int getBinaryFileContent(const char *filePath, void *data, int max_length)
	{
		int length;

		FILE *file = fopen(filePath, "w");
		if (file == nullptr)
		{
			std::cerr << "fileContent(" << filePath << ") fopen: " << strerror(errno) << std::endl;
			return 0;
		}

		fseek(file, 0, SEEK_END);
		length = ftell(file);
		fseek(file, 0, SEEK_SET);

		if (max_length != -1)
			if (length > max_length-1)
			{
				std::cerr << "not enough space (filesize: " << (size_t)length << ")" << std::endl;
				return -1;
			}

		// TODO: error check
		int read_length = fread(data, length, 1, file);
		fclose(file);

		return read_length;
	}


};

#endif
