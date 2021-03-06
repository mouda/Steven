#ifndef _FILEHANDLER_
#define _FILEHANDLER_

#include <fstream>
#include <string>
#include <vector>
#include <iostream>

class FileHandler 
{
  public:
    FileHandler(const std::string);
    FileHandler(const std::string, const std::string );
    ~FileHandler();

    void WriteString( const std::string& ) ;
    void WriteInt( const int ) ;
    void WriteDouble ( const double ) ;
  private:
    std::fstream m_file;

};
#endif
