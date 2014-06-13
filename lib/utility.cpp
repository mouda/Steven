#include "utility.h"
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <linux/netdevice.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <unistd.h>
#include <unistd.h>
#include <limits>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <boost/filesystem.hpp>

using std::string;
using std::cout;
using std::endl;
using std::cerr;

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

int print_addresses(const int domain, string& ipAddr)
{
  int s;
  struct ifconf ifconf;
  struct ifreq ifr[50];
  int ifs;
  int i;
  string strInterfaceName; 
  string strIpAddr;

  s = socket(domain, SOCK_STREAM, 0);
  if (s < 0) {
    cerr << "socket" << endl;
    return 0;
  }

  ifconf.ifc_buf = (char *) ifr;
  ifconf.ifc_len = sizeof ifr;

  if (ioctl(s, SIOCGIFCONF, &ifconf) == -1) {
    cerr << "ioctl" << endl;
    return 0;
  }

  ifs = ifconf.ifc_len / sizeof(ifr[0]);
  for (i = 0; i < ifs; i++) {
    char ip[INET_ADDRSTRLEN];
    struct sockaddr_in *s_in = (struct sockaddr_in *) &ifr[i].ifr_addr;

    if (!inet_ntop(domain, &s_in->sin_addr, ip, sizeof(ip))) {
      cerr << "net_ntop" << endl;
      return 0;
    }
    strInterfaceName.assign(ifr[i].ifr_name);
    strIpAddr.assign(ip);
    if ( strInterfaceName == "eth0" ) {
      break;
    }

  }
  string dlim("192.168.21.");
  ipAddr = strIpAddr.substr(dlim.length(),strIpAddr.length());

  close(s);

  return 1;
}
