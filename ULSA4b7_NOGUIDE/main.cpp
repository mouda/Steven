#include<iostream>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<cassert>

#include "ULSA4b7_NOGUIDE.h"

#define SA_INI_TEMP 3.0
#define SA_FIN_TEMP 0.5
using namespace std;

int print_addresses(const int domain, string& ipAddr);

int main(int argc, char* argv[])
{
  //********************************************//
  //User-assigned parameters                    //
  //********************************************//

  if(argc!=14)
  {
    cout <<"Usage: ULSA_DC_Estimate "         << endl;
    cout <<"\t 1 [ totalNodes         ]" << endl;
    cout <<"\t 2 [ maxChNum           ]" << endl;
    cout <<"\t 3 [ quantizationBits   ]" << endl;
    cout <<"\t 4 [ PowerMax(dbm)      ]" << endl;
    cout <<"\t 5 [ Bandwidth(kHz)     ]" << endl;
    cout <<"\t 6 [ FilePath           ]" << endl;
    cout <<"\t 7 [ Correlation Factor ]" << endl;
    cout <<"\t 8 [ MapIndex           ]" << endl;
    cout <<"\t 9 [ outputIndex        ]" << endl;
    cout <<"\t 10[ Simulation Time    ]" << endl;
    cout <<"\t 11[ Fidelity Ratio     ]" << endl;
    cout <<"\t 12[ Structure and Detail: 0 to turn off; 1 to turn on ]" << endl;
    cout <<"\t 13[ Number of Iteration]"  << endl;


    cout<<"Example:Usage: winULSA4b7_NOGUIDE 195 10 10 10 180 mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_2.txt 10000 2 2 3"<<endl;
    cout<<"Output 1: Struc, Detail, Matrics"<<endl;
    cout<<"Output 2: Struc, Detail, Metrics, Metrics Every Round"<<endl;
    cout<<"IniFlag: kmeans , HeadLimites"<<endl;
    return 0;
  }
  const int maxChNum= atoi(argv[2]);
  const int totalNodes=atoi(argv[1]);
  const double quantizationBits=atof(argv[3]);
  const double powerMaxDbm =atof(argv[4]);
  const double bandwidthKhz =atof(argv[5]);
  //const double correlationFactor=atof(argv[7]);
  const double compressionRatio=atof(argv[7]);
  const int outputControl=atoi(argv[9]);
  const int simulationTime=atoi(argv[10]);
  const double fidelityRatio=atof(argv[11]);

  const int isDetailOn=atoi(argv[12]);

  //********************************************//
  // Set internal parameters                    //
  //********************************************//
  srand (time(NULL));
  //srand(10);
  clock_t begin, end;
  //---------------------//
  // Control Constant    //
  //---------------------//
  char iniFlag[]="kmedoids_distance";
  const float powerMaxWatt = pow(10,((float)powerMaxDbm)/10) /1000;//we don't need to divide the BW(bandwidthKhz*1000);//Unit := Watt
  const int SAIter=atoi(argv[13]);
  //****************************//
  //Exception Handling Variable //
  //****************************//
  int totalNodesCheck = -1;

  //---------------------------//
  //Input Variable             //
  //---------------------------//
  double radius=0;

  //---------------------------//
  //@Data Parsing              //
  //---------------------------//
  FILE *fid;
  fid = fopen(argv[6], "r");
  if(fid==NULL)
  {
    cout << "Read map file failed!!\n";
    return -1;
  }
  fscanf(fid, "%d %lf", &totalNodesCheck,&radius);
 // double density= static_cast<double>(totalNodes)*1e6/(3.1415*radius*radius);
  if (totalNodes != totalNodesCheck)
  {
    cout<<"Mapfile Error. \
    Please check the identication between content and file name"<<endl;
    return -1;
  }

  //---------------------------//
  //@Get IP address            //
  //---------------------------//
  
  int domains[] = { AF_INET, AF_INET6 };
  string ipAddr;
  if (!print_addresses(domains[0], ipAddr))
    return 1;

  double alpha =pow (10, -log10(SA_INI_TEMP/SA_FIN_TEMP)/SAIter);

  //ULSA4b7_NOGUIDE *toolSA = new ULSA4b7_NOGUIDE(fid, totalNodes, maxChNum, SAIter,outputControl, SA_INI_TEMP, alpha, correlationFactor);
  //ULSA4b7_NOGUIDE toolSA(fid, totalNodes, maxChNum, SAIter,outputControl, SA_INI_TEMP, alpha, correlationFactor);
  ULSA4b7_NOGUIDE toolSA(fid, totalNodes, maxChNum, SAIter,outputControl, isDetailOn,SA_INI_TEMP, alpha, compressionRatio,ipAddr );


  toolSA.radius=radius;
  begin = clock();
  for(int i=0;i<simulationTime;i++)
  {
    if(!toolSA.setSystem(powerMaxWatt, quantizationBits,bandwidthKhz, fidelityRatio))
    {cout<<"Set parameter failed"<<endl;return -1;}
    if (!toolSA.setInitialStucture (iniFlag))
    {cout<<"The "<<argv[5]<<" can't be initialized by"<<iniFlag<<endl;continue;}

    toolSA.startCool();
  }

  end = clock();
  fclose(fid);
  printf("Computing Time %f \n",((float)(end-begin))/CLOCKS_PER_SEC);
  toolSA.releaseMemory();
  return 1;
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
