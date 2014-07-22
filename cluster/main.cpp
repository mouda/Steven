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
#include "minPowerCsFactory.h"
#include "minPowerImageCsFactory.h"
#include "baselineImageCsFactory.h"
#include "kmeansCsFactory.h"
#include "simulator.h"
#include "fileHandler.h"
#include "minResCsFactory.h"
#include "nonGuidedCsFactory.h"

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
  double  tier1TxTime;
  //double  endTime;
  int     tier2NumSlot;
  int     SAIter;
  double  fidelityRatio;
  string  mapFileName;
  /* output file name string */
  string  CSFName;
  string  MetisFName;
  string  CSFormation;
  string  labelFName;
  string  WeightMatrixFName;
  bool    imageFlag = false;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "Produce help message")
      ("bandwidth,b",             po::value<double>(),  "Bandwidth (kHz) ")
      ("power,p",                 po::value<double>(),  "Maximum Power (dbm) ")
      ("quantization,q",          po::value<double>(),  "Bits of quantization")
      ("map,m",                   po::value<string>(),  "Map file name")
      ("nodes,n",                 po::value<int>(),     "Initial number of nodes")
      ("heads,H",                 po::value<int>(),     "Initial number of heads")
      ("txTime,t",                po::value<double>(),  "Transmission time per slot (ms) ")
      ("tier1TxTime,I",           po::value<double>(),  "Tier 1 tx time (ms)")
      ("fidelity,f",              po::value<double>(),  "Fidelity ratio")
      ("iteration,i",             po::value<int>(),     "Number of iteration Simulate Annealing")
      ("spatialCorrelation,c",    po::value<double>(),  "Spatial Correlation level")
      ("temporalCorrelation,T",   po::value<double>(),  "Temproal Correlation factor")
      ("tier2NumSlot,N",          po::value<int>(),     "Number of tier-2 slots")
      ("ClusterStructureOutput,C",po::value<string>(),  "Cluster structure output file name")
      ("WeightedMatrix,W",        po::value<string>(),  "Matrix form of Weght matrix")
      ("ClusterFormation,F",      po::value<string>(),  "Cluster Formation Algorithm")
      ("iterationlog,l",          "Iteration log");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.size() == 0 || vm.count("help")) {
      cout << desc << "\n";
      return 0;
    } else if(vm.size() == 14  || vm.size() == 15 ) {

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
      tier2NumSlot   =          vm["tier2NumSlot"].as<int>();
      tier1TxTime    =          vm["tier1TxTime"].as<double>();
      CSFormation =             vm["ClusterFormation"].as<string>();

      /* output control */
      if (vm.count("ClusterStructureOutput")) {
        CSFName = vm["ClusterStructureOutput"].as<string>();
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
      ImageSource*  myImageSource = 0;
      myMap = myMapFactory.CreateMap(imageFlag);
      if (CSFormation != "ImageSource" && CSFormation != "ImageBaseline" ) {
        myMatComputer = myMapFactory.CreateMatrixComputer();
      }
      if (!myMap ) {
        cerr << "Error: Failed to initialize map" << endl;
        return 1;
      }
      if (CSFormation != "ImageSource" && CSFormation != "ImageBaseline" && !myMatComputer) {
        cerr << "Error: Failed to initialize correlation comulter" << endl;
        return 1;
      }
      ClusterStructure* myCS = 0;
      CsFactory* myCsFactory = 0;
      
      if (vm.count("ClusterStructureOutput")) {
        if (CSFormation == "MinRes"){
          myCsFactory = new MinResCsFactory(myMap, myMatComputer); 
          (dynamic_cast<MinResCsFactory*>(myCsFactory))->SetCompressionRatio(spatialCompressionRatio);
          (dynamic_cast<MinResCsFactory*>(myCsFactory))->SetMapFileName(mapFileName);
          (dynamic_cast<MinResCsFactory*>(myCsFactory))->SetMapFileName(mapFileName);
          (dynamic_cast<MinResCsFactory*>(myCsFactory))->SetTier2NumSlot(tier2NumSlot);
          (dynamic_cast<MinResCsFactory*>(myCsFactory))->SetTier1TxTime(tier1TxTime);
        } else if (CSFormation == "MinPower") {
          myCsFactory = new MinPowerCsFactory(myMap, myMatComputer);
          (dynamic_cast<MinPowerCsFactory*>(myCsFactory))->SetCompressionRatio(spatialCompressionRatio);
          (dynamic_cast<MinPowerCsFactory*>(myCsFactory))->SetMapFileName(mapFileName);
          (dynamic_cast<MinPowerCsFactory*>(myCsFactory))->SetFidelityRatio(fidelityRatio);
          (dynamic_cast<MinPowerCsFactory*>(myCsFactory))->SetTier2NumSlot(tier2NumSlot);
          (dynamic_cast<MinPowerCsFactory*>(myCsFactory))->SetTier1TxTime(tier1TxTime);
          if (vm.count("iterationlog")) {
            (dynamic_cast<MinPowerCsFactory*>(myCsFactory)->SetIterationLog(true));
          }
        } else if (CSFormation == "NonGuided") {
          myCsFactory = new NonGuidedCsFactory(myMap, myMatComputer);
          (dynamic_cast<NonGuidedCsFactory*>(myCsFactory))->SetCompressionRatio(spatialCompressionRatio);
          (dynamic_cast<NonGuidedCsFactory*>(myCsFactory))->SetMapFileName(mapFileName);
          (dynamic_cast<NonGuidedCsFactory*>(myCsFactory))->SetFidelityRatio(fidelityRatio);
          (dynamic_cast<NonGuidedCsFactory*>(myCsFactory))->SetTier2NumSlot(tier2NumSlot);
          (dynamic_cast<NonGuidedCsFactory*>(myCsFactory))->SetTier1TxTime(tier1TxTime);
          if (vm.count("iterationlog")) {
            (dynamic_cast<NonGuidedCsFactory*>(myCsFactory)->SetIterationLog(true));
          }
        } else if (CSFormation == "ImageSource") {
          myImageSource = myMapFactory.CreateImageSource();
          myCsFactory = new MinPowerImageCsFactory(myMap, myImageSource);
          (dynamic_cast<MinPowerImageCsFactory*>(myCsFactory))->SetCompressionRatio(spatialCompressionRatio);
          (dynamic_cast<MinPowerImageCsFactory*>(myCsFactory))->SetMapFileName(mapFileName);
          (dynamic_cast<MinPowerImageCsFactory*>(myCsFactory))->SetFidelityRatio(fidelityRatio);
          (dynamic_cast<MinPowerImageCsFactory*>(myCsFactory))->SetTier2NumSlot(tier2NumSlot);
          (dynamic_cast<MinPowerImageCsFactory*>(myCsFactory))->SetTier1TxTime(tier1TxTime);
          if (vm.count("iterationlog")) {
            (dynamic_cast<MinPowerImageCsFactory*>(myCsFactory)->SetIterationLog(true));
          }
        } else if (CSFormation == "ImageBaseline") {
          myImageSource = myMapFactory.CreateImageSource();
          myCsFactory = new BaselineImageCsFactory(myMap, myImageSource);
          (dynamic_cast<BaselineImageCsFactory*>(myCsFactory))->SetCompressionRatio(spatialCompressionRatio);
          (dynamic_cast<BaselineImageCsFactory*>(myCsFactory))->SetMapFileName(mapFileName);
          (dynamic_cast<BaselineImageCsFactory*>(myCsFactory))->SetFidelityRatio(fidelityRatio);
          (dynamic_cast<BaselineImageCsFactory*>(myCsFactory))->SetTier2NumSlot(tier2NumSlot);
          (dynamic_cast<BaselineImageCsFactory*>(myCsFactory))->SetTier1TxTime(tier1TxTime);
          if (vm.count("iterationlog")) {
            (dynamic_cast<BaselineImageCsFactory*>(myCsFactory)->SetIterationLog(true));
          }
        }
        if (!myCsFactory) {
          cerr << "Error: Failed to initalize cluster structure factory" << endl;
          return 1;
        }

        myCS = myCsFactory->CreateClusterStructure();

        if (!myCS) {
          cerr << "Error: Failed to initalize cluster structure" << endl;
          return 1;
        }


        Simulator mySimulator(
            myMap, 
            myCS, 
            myMatComputer
            );

        /* output control */
        mySimulator.WriteCS( CSFName );

        if (CSFormation == "ImageBaseline") {
          mySimulator.WriteWorseCaseTier2Power("ImageBaselineTier2Power.out", txTimePerSlot , powerMaxWatt); 
        }
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
