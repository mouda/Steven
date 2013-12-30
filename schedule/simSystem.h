#ifndef SimSystem_H
#define SimSystem_H
#include <cmath>
#include "ULAGENT.h"
class SimSystem
{
  public:
    //SimSystem();
    double returnChannelGain_BS(ULAGENT &node1);
    double returnChannelGain_2Nodes(ULAGENT &node1, ULAGENT &node2);
    double returnChannelGainByPos( double lhsX, double lhsY, double rhsX, double rhsY);
    double returnRate_BS(ULAGENT &node1, double bandwidthKhz,double PowerWatt);
    double returnInBandThermalNoise(double bandwidthKhz);




  private:
  const static float BSx=0,BSy=0;
  //const static float n0 = 3.981071705534985e-21;// = pow(10,-17.4)/1000;Watt, ref:Wikipwdia Thermal Noise
  const static float n0 = 1e-19;// = pow(10,-17.4)/1000;Watt, ref:Wikipwdia Thermal Noise

  const static float pathLoss0_MacroUE = 131.1;// no unit ref:LTE release 9
  const static float pathLossAlpha_MacroUE =42.68; //no unit ref:LTE release 9
  const static float pathLoss0_MacroUE_Gib = 131.1;// For simulation Adjustment
  const static float pathLossAlpha_MacroUE_Gib =42.68; //
};




#endif
