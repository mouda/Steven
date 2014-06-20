#include "fileHandler.h"

FileHandler::FileHandler(const std::string fileName)
{
  m_file.open(fileName.c_str(), std::fstream::out | std::fstream::app);

}
FileHandler::FileHandler(const std::string fileName, const std::string option)
{
  if ( option == "app") {
    m_file.open(fileName.c_str(), std::fstream::app );
  }
  else if ( option == "out") {
    m_file.open(fileName.c_str(), std::fstream::out );
  }
  else {
    std::cerr << "Error Flags setted" << std::endl;
  }
}

FileHandler::~FileHandler()
{
  m_file.close();

}

void
FileHandler::WriteString( const std::string& str )
{
  m_file << str << std::endl;
}

void
FileHandler::WriteInt( const int rhs) 
{
  m_file << rhs << std::endl;
}

void
FileHandler::WriteDouble( const double rhs)
{
  m_file << rhs << std::endl;
}
