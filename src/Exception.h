// File: Exception.h
// Purpose: Header file for exception handling using standard exceptions

#pragma once

#include <iostream>
#include <exception>
#include <fstream>
#include <string>

namespace nimbus
{

	/**
	Exception handling class using STL standard exceptions
	*/
	class NimbusException : public std::exception
	{
	private:
		std::string errMsg; // Error message
		std::string file; // File where the error was produced
		std::string function; // Function where the error was produced
		std::string line;	  // Line number where the error was produced
	public:
		NimbusException(const char* errMsg, const char* file, const char* function, int line);
		void printError() const throw();
		~NimbusException() _NOEXCEPT {}
	};
}