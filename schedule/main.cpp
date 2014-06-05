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
#include "kmeansCsFactory.h"
#include "minResCsFactory.h"
#include "minPowerSACluster.h"
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
  int     maxChNum;
  int     totalNodes;
  double  quantizationBits;
  double  powerMaxDbm;
  //we don't need to divide the BW(bandwidthKhz*1000);//Unit := Watt
  double  bandwidthKhz;
  double  spatialCompressionRatio;
  double  temporalCorrFactor;
  double  txTimePerSlot;
  //double  endTime;
  int     tier2NumSlot;
  int     SAIter;
  double  fidelityRatio;
  string  mapFileName;
  /* output file name string */
  string  CSFName;
  string  entropyFName;
  string  totalEntropyFName;
  string  MSEFName;
  string  solutionFName;
  string  supportFName;
  string  powerFName; 
  string  strAlgFlag;
  string  CSFormation;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "Produce help message")
      ("bandwidth,b",             po::value<double>(),  "Bandwidth (kHz) ")
      ("power,p",                 po::value<double>(),  "Maximum Power (dbm) ")
      ("quantization,q",          po::value<double>(),  "Bits of quantization")
      ("map,m",                   po::value<string>(),  "Map file name")
      ("algorithm,A",             po::value<string>(),  "Scheduling algorithm Baseline|Algorithm ")
      ("nodes,n",                 po::value<int>(),     "Initial number of nodes")
      ("heads,H",                 po::value<int>(),     "Initial number of heads")
      ("txTime,t",                  po::value<double>(),  "Transmission time per slot (ms) ")
      ("fidelity,f",              po::value<double>(),  "Fidelity ratio")
      ("iteration,i",             po::value<int>(),     "Number of iteration Simulate Annealing")
      ("spatialCorrelation,c",    po::value<double>(),  "Spatial Correlation level")
      ("temporalCorrelation,T",   po::value<double>(), "Temproal Correlation factor")
      ("tier2NumSlot,N",          po::value<int>(), "Number of tier-2 slots")
      ("ClusterStructureOutput,C",po::value<string>(), "Cluster structure output file name")
      ("TotalEntropy,O",          po::value<string>(), "Total entropy per slot")
      ("EntropyOutput,E",         po::value<string>(), "Entropy output file name")
      ("SupportOutput,U",         po::value<string>(), "Support number output file name")
      ("MSEOutput,M",             po::value<string>(), "MSE output file name")
      ("SolutionOutput,S",        po::value<string>(), "Solution output file name")
      ("PowerOutput,P",           po::value<string>(), "Power output file name")
      ("ClusterFormation,F",      po::value<string>(), "Cluster Formation Algorithm");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.size() == 0 || vm.count("help")) {
      cout << desc << "\n";
      return 0;
    } else if(vm.size() > 13 && vm.size() <= 22 ) {

      totalNodes =              vm["nodes"].as<int>();
      maxChNum =                vm["heads"].as<int>();
      quantizationBits =        vm["quantization"].as<double>();
      powerMaxDbm =             vm["power"].as<double>();
      bandwidthKhz =            vm["bandwidth"].as<double>();
      txTimePerSlot =           vm["txTime"].as<double>();
      spatialCompressionRatio = vm["spatialCorrelation"].as<double>();
      temporalCorrFactor =      vm["temporalCorrelation"].as<double>();
      fidelityRatio =           vm["fidelity"].as<double>();
      mapFileName =             vm["map"].as<string>();
      strAlgFlag =              vm["algorithm"].as<string>();
      tier2NumSlot   =          vm["tier2NumSlot"].as<int>();
      CSFormation =             vm["ClusterFormation"].as<string>();

      /* output control */
      if (vm.count("ClusterStructureOutput")) {
        CSFName = vm["ClusterStructureOutput"].as<string>();
      }
      if (vm.count("EntropyOutput")) {
        entropyFName = vm["EntropyOutput"].as<string>();
      }
      if (vm.count("MSEOutput")) {
        MSEFName = vm["MSEOutput"].as<string>();
      }
      if (vm.count("SolutionOutput")) {
        solutionFName = vm["SolutionOutput"].as<string>();
      }
      if (vm.count("SupportOutput")) {
        supportFName = vm["SupportOutput"].as<string>();
      }
      if (vm.count("TotalEntropy")) {
        totalEntropyFName = vm["TotalEntropy"].as<string>();
      }
      if (vm.count("PowerOutput")) {
        powerFName = vm["PowerOutput"].as<string>();
      }

      double powerMaxWatt = pow(10,(powerMaxDbm)/10) /1000;

      MapFactory myMapFactory(
          mapFileName, 
          powerMaxWatt, 
          spatialCompressionRatio, 
          temporalCorrFactor,
          quantizationBits, 
          bandwidthKhz, 
          maxChNum, 
          totalNodes
          );

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
      CsFactory* myCsFactory = 0;
      
      if (CSFormation == "Kmeans") {
        myCsFactory = new KmeansCsFactory(myMap, myMatComputer);
      }
      else if (CSFormation == "MinRes"){
        myCsFactory = new MinResCsFactory(myMap, myMatComputer); 
        (dynamic_cast<MinResCsFactory*>(myCsFactory))->SetCompressionRatio(spatialCompressionRatio);
        (dynamic_cast<MinResCsFactory*>(myCsFactory))->SetMapFileName(mapFileName);
        (dynamic_cast<MinResCsFactory*>(myCsFactory))->SetFidelityRatio(fidelityRatio);
      }
      else if (CSFormation == "MinPower") {
      }

      myCS = myCsFactory->CreateClusterStructure();

      if (!myCS) {
        cerr << "Error: Failed to initalize cluster structure" << endl;
        return 1;
      }

      SchedulerFactory mySchedFactory(txTimePerSlot, tier2NumSlot, bandwidthKhz, myMap, myMatComputer, myCS);
      Scheduler* myScheduler = 0;
      myScheduler = mySchedFactory.CreateScheduler(strAlgFlag);
      if (!myScheduler) {
        cerr << "Error: Failed to initialize scheduler" << endl;
        return 1;
      }
      Simulator mySimulator(
          myMap, 
          myCS, 
          myScheduler, 
          myMatComputer, 
          entropyFName, 
          MSEFName, 
          solutionFName, 
          supportFName,
          totalEntropyFName,
          powerFName 
          );
      mySimulator.SelfCheck();
      mySimulator.SequentialRun(tier2NumSlot);

      /* output control */
      if (vm.count("ClusterStructureOutput")) {
        mySimulator.WriteCS( CSFName );
      }
      if (vm.count("EntropyOutput")) {
        mySimulator.WriteEntropy();
      }
      if (vm.count("MSEOutput")) {
        mySimulator.WriteMSE();
      }
      if (vm.count("SolutionOutput")) {
        mySimulator.WriteSolution();
      }
      if (vm.count("SupportOutput")) {
        mySimulator.WriteSupport();
      }
      if (vm.count("TotalEntropy")) {
        mySimulator.WriteTotalEntropy();
      }
      if (vm.count("PowerOutput")) {
        mySimulator.WriteVecPower();
      }

      delete myCsFactory;
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
