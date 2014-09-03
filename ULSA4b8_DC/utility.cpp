#include "utility.h"
#include <limits>
#include <sstream>
#include <boost/filesystem.hpp>

bool FileExist( const std::string& Name)
{
  return boost::filesystem::exists(Name); 
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}
std::istream& Ignoreline(std::ifstream& in, std::ifstream::pos_type& pos)
{
  pos = in.tellg();
  return in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

std::string GetLastLine(std::ifstream& in)
{
  std::ifstream::pos_type pos = in.tellg();

  std::ifstream::pos_type lastPos;
  while (in >> std::ws && Ignoreline(in, lastPos))
    pos = lastPos;

  in.clear();
  in.seekg(pos);

  std::string line;
  std::getline(in, line);
  in.clear();
  in.seekg(std::ios::beg);
  return line;
}
