#ifndef SABASE_H
#define SABASE_H
#include "../../commonLibrary/SimSystem.h"
#include <cstdio>
#include <cstdlib>
class SABASE{
  public:
  SABASE(){};
  SABASE(FILE *fileReadCursor, int inputTotalNodes, int inputMaxChNum,int inputSAFac,  \
               int inOutputControl, double inputTemprature, double InputSaAlpha, \
               double inCorrelationFactor);
  bool setSystem(float inPowerMaxWatt, int inQuantizationBits,double inBandwidthKhz, double fidelity);
  bool setInitialStucture(char* inputFlag);
  bool setIniStruKmeans();//not public but related to setIniStrucKmeans
  bool setIniHeadLimited();
  //@setby Constructor
  int SAIter; //Number of iterations
  int outCtrl;
  int maxChNum;
  int totalNodes;
  double temparature;
  double constantIniTemprature;
  double alpha;
  double C2;
};
#endif
