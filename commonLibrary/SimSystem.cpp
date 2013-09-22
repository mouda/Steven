//@2013/1/11 SimSystem is simply providing some simulation parameter constuction utility

#include<iostream>
#include "SimSystem.h"
using namespace std;

/*SimSystem::SimSystem(){
 cout<<"  sd"<<endl;
 cout<<"Simsystem: Thermal Noise "<<n0<<" Watt/Hz"<<endl;
}*/

double SimSystem::returnChannelGain_BS(ULAGENT &nodes1){
    float tempDib = pow(BSx - nodes1.locX,2) +pow(BSy - nodes1.locY,2);
    float pathLossDBMacroUE = pathLoss0_MacroUE + pathLossAlpha_MacroUE * (0.5*log10(tempDib) - 3);
    double channelGain = pow(10,-(pathLossDBMacroUE/10));
    return channelGain;
}
double SimSystem::returnChannelGain_2Nodes(ULAGENT &nodes1, ULAGENT &nodes2){
    float tempDib = pow(nodes1.locX - nodes2.locX,2) +pow(nodes1.locY - nodes2.locY,2);
    float pathLossDBMacroUE = pathLoss0_MacroUE + pathLossAlpha_MacroUE * (0.5*log10(tempDib) - 3);
    double channelGain = pow(10,-(pathLossDBMacroUE/10));
    return channelGain;
}

double SimSystem::returnRate_BS(ULAGENT &nodes1, double bandwidthKhz, double powerWatt){
    float tempDib = pow(BSx - nodes1.locX,2) +pow(BSy - nodes1.locY,2);
    float pathLossDBMacroUE = pathLoss0_MacroUE_Gib + pathLossAlpha_MacroUE_Gib * (0.5*log10(tempDib) - 3);
    double channelGain = pow(10,-(pathLossDBMacroUE/10));
    double rate = bandwidthKhz*1e3*log2(1+powerWatt*channelGain/(bandwidthKhz*1e3*n0));//(bps)
    return rate;
}
double SimSystem::returnInBandThermalNoise(double bandwidthKhz){
    double inBandNoise = bandwidthKhz*1e3*n0;
    return inBandNoise;
}
