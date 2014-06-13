#ifndef _UTILITY_
#define _UTILITY_

#include <string>
#include <fstream>
#include <boost/filesystem.hpp>

using std::string;
bool FileExist( const std::string& Name);

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> split(const std::string &s, char delim); 
/* file operations */

std::istream& Ignoreline(std::ifstream& in, std::ifstream::pos_type& pos);
std::string GetLastLine(std::ifstream& in);
int print_addresses(const int domain, string& ipAddr);
#endif
