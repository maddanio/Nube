// File: Exception.CPP
// Purpose: Implementation file for exception handling using standard exceptions

#include "Exception.h"

using namespace nimbus;

/**
Constructor
@param errMsg Error message
@param file File where the error was located
@param function Class method where the error was located
@param line Sourc code line where the error was encountered
*/

NimbusException::NimbusException(const char* errMsg, const char* file, const char* function, int line)
{
	this->errMsg = errMsg;
	this->file = file;
	this->function = function;
	this->line = std::to_string(line);
}

/**
Main method for error information
*/
void NimbusException::printError() const throw()
{
	std::cout << "Error in file: " << file << " in function: " << function << " in line: " << line << std::endl;
	std::cout << "Cause: " << errMsg << std::endl;
}