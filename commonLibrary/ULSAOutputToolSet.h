#ifndef ULSAOutputToolSet_H
#define ULSAOutputToolSet_H
#include <cstdio>
#include <string>
#include "../commonLibrary/ULCS1b.h"
template <typename T>
class ULSAOutputToolSet{

    public:
    //ULSAOutputToolSet(T *inMyObj);
    //T myObj;
    void writeIniHead(char * filename1, T &obj);
    void writeBestStru(char* filename3, T &obj);
    void writeBestStru_V2(char* filename3, T &obj);
    void writeAll(char* filename4, double density,T &obj);
    void writeClusterInfo (char *filename5, T &obj,  char * time);
    void writeEnergy(char* filename3, T &obj);
    void summaryNwrite2tiers(char*filename, T&obj);
    void writePeformance_MinResors_with2ndPowerControl(char *filename,T &myobj, double best2nd_ms, double fidelityRatio);
    void writePeformance_MinResors_with2ndPowerControl_4b(char *filename,T &myobj, double best2nd_ms, double fidelityRatio,int bestChNum);
    void summaryNwrite2tiers_MinResors_with2ndPowerControl(char *filename,T &myobj, double best2nd_ms);
    void showIndEntropy(T &myobj);

    std::string getIpAddr();

    double best2ndJoule;
    double bestTotalJoule;
    double bestTotal_ms;
};



#endif
