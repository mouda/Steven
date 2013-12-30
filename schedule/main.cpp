#include <iostream>
#include <string>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include "map.h"
#include "mapFactory.h"
#include "clusterStructure.h"
#include "csFactory.h"
#include "scheduler.h"
#include "schedFactory.h"
#include "simulator.h"

#define SA_INI_TEMP 3.0
#define SA_FIN_TEMP 0.5

using namespace std;
int main(int argc, char *argv[])
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


    cout<<"Output 1: Struc, Detail, Matrics"<<endl;
    cout<<"Output 2: Struc, Detail, Metrics, Metrics Every Round"<<endl;
    cout<<"IniFlag: kmeans , HeadLimites"<<endl;
    return 0;
  }

  const int maxChNum= atoi(argv[2]);
  const int totalNodes=atoi(argv[1]);
  const double quantizationBits=atof(argv[3]);
  const double powerMaxDbm =atof(argv[4]);
  const float powerMaxWatt = pow(10,((float)powerMaxDbm)/10) /1000;
  //we don't need to divide the BW(bandwidthKhz*1000);//Unit := Watt
  const double bandwidthKhz =atof(argv[5]);
  const double compressionRatio=atof(argv[7]);
  const int outputControl=atoi(argv[9]);
  const int simulationTime=atoi(argv[10]);
  const double fidelityRatio=atof(argv[11]);
  const int isDetailOn=atoi(argv[12]);
  string mapFileName(argv[6]);


  MapFactory myMapFactory(mapFileName, powerMaxWatt, compressionRatio, maxChNum, totalNodes);
  Map* myMap = myMapFactory.CreateMap();
  CORRE_MA_OPE* myMatrxComputer = myMapFactory.CreateMatrixComputer();
  CsFactory myCsFactory;
  SchedulerFactory mySchedFactory;

  ClusterStructure* myCS = myCsFactory.CreateClusterStructure();
  Scheduler* myScheduler = mySchedFactory.CreateScheduler();

  Simulator mySimulator(*myMap, *myCS, *myScheduler);
  mySimulator.Run();
  
  return 0;
}
