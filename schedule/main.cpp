#include <iostream>
#include <string>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <boost/program_options.hpp>

#include "map.h"
#include "mapFactory.h"
#include "clusterStructure.h"
#include "csFactory.h"
#include "scheduler.h"
#include "maxSNRScheduler.h"
#include "schedFactory.h"
#include "simulator.h"

#define SA_INI_TEMP 3.0
#define SA_FIN_TEMP 0.5

namespace po = boost::program_options;

using namespace std;
int main(int argc, char *argv[])
{
  //********************************************//
  //User-assigned parameters                    //
  //********************************************//
  int maxChNum;
  int totalNodes;
  double quantizationBits;
  double powerMaxDbm;
  //we don't need to divide the BW(bandwidthKhz*1000);//Unit := Watt
  double bandwidthKhz;
  double compressionRatio;
  double txTimePerSlot;
  int SAIter;
  double fidelityRatio;
  string mapFileName;
  string strAlgFlag;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "Produce help message")
      ("bandwidth,b",     po::value<double>(),  "Bandwidth (kHz) ")
      ("power,p",         po::value<double>(),  "Maximum Power (dbm) ")
      ("quantization,q",  po::value<double>(),  "Bits of quantization")
      ("map,m",           po::value<string>(),  "Map file name")
      ("algorithm,G",     po::value<string>(),  "Scheduling algorithm Baseline|Algorithm ")
      ("nodes,n",         po::value<int>(),     "Initial number of nodes")
      ("heads,H",         po::value<int>(),     "Initial number of heads")
      ("time,t",          po::value<double>(),  "Transmission time per slot (ms) ")
      ("fidelity,f",      po::value<double>(),  "Fidelity ratio")
      ("iteration,i",     po::value<int>(),     "Number of iteration Simulate Annealing")
      ("correlation,c",   po::value<double>(),  "Correlation level");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.size() == 0 || vm.count("help")) {
      cout << desc << "\n";
      return 0;
    } else if(vm.size() == 10 ) {

      totalNodes =        vm["nodes"].as<int>();
      maxChNum =          vm["heads"].as<int>();
      quantizationBits =  vm["quantization"].as<double>();
      powerMaxDbm =       vm["power"].as<double>();
      bandwidthKhz =      vm["bandwidth"].as<double>();
      txTimePerSlot =     vm["time"].as<double>();
      compressionRatio =  vm["correlation"].as<double>();
      fidelityRatio =     vm["fidelity"].as<double>();
      mapFileName =       vm["map"].as<string>();
      strAlgFlag =        vm["algorithm"].as<string>();

      double powerMaxWatt = pow(10,(powerMaxDbm)/10) /1000;

      MapFactory myMapFactory(mapFileName, powerMaxWatt, 
          compressionRatio, quantizationBits, bandwidthKhz, maxChNum, totalNodes);
      Map* myMap = 0;
      CORRE_MA_OPE* myMatComputer  = 0;
      myMap = myMapFactory.CreateMap();
      myMatComputer = myMapFactory.CreateMatrixComputer();
      if (!myMap ) {
        cerr << "Error: Failed to initialize map" << endl;
        return 1;
      }
      if (!myMatComputer) {
        cerr << "Error: Failed to initialize correlation comulter" << endl;
        return 1;
      }
      ClusterStructure* myCS = 0;
      CsFactory myCsFactory(myMap, myMatComputer);
      myCS = myCsFactory.CreateClusterStructure();

      if (!myCS) {
        cerr << "Error: Failed to initalize cluster structure" << endl;
        return 1;
      }

      SchedulerFactory mySchedFactory(0.01, bandwidthKhz, myMap, myMatComputer, myCS);
      Scheduler* myScheduler = 0;
      myScheduler = mySchedFactory.CreateScheduler(strAlgFlag);
      if (!myScheduler) {
        cerr << "Error: Failed to initialize scheduler" << endl;
        return 1;
      }
      Simulator mySimulator(myMap, myCS, myScheduler, myMatComputer);
      mySimulator.SelfCheck();
      mySimulator.Run(1);
    }
    else {
      cout << desc << "\n";
      return 0;
    }

  }
  catch(std::exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!\n";
  }

  
  return 0;
}
