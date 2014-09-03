#include <cmath>
#include <iostream>
#include <cstdio>
#include <cassert>
#include<cfloat>
using namespace std;
#include "ULConstraintSolver.h"

ULConstraintSolver::ULConstraintSolver(int inChNum,int inTotalNodes,float inPowerMax, float ininBandNoise,float inC2,vector<int>& inVecHeadName,float** inGij, \
                                       double* inputNextNodePower,list<list<int> >* inListCluMem)
{

    chNum = inChNum;
    totalNodes = inTotalNodes;
    powerMax=inPowerMax*scale;
    inBandNoise = ininBandNoise*scale;
    C2 = inC2;
    vecHeadName = &inVecHeadName;
    Gij =inGij;
    bandwidthKhz = -1;//Not using this parameter if using this constructor
    maIndexInterference = new int* [chNum];
    for (int i=0; i<chNum; i++)maIndexInterference[i] = new int [chNum];
    maStrengthInterference = new double* [chNum];
    for (int i=0; i<chNum; i++)maStrengthInterference[i] = new double [chNum];

    nextNodePower = inputNextNodePower;
    listCluMember = inListCluMem;


    powerDifference = new double [totalNodes];
    for(int i=0; i<totalNodes; i++)powerDifference[i] = 0.0;
    powerDifferenceRatio = new double [totalNodes];
    for(int i=0; i<totalNodes; i++)powerDifferenceRatio[i] = 0.0;

}

//Constructor for solving minimize T2 find
ULConstraintSolver::ULConstraintSolver(int inChNum,int inTotalNodes,float inPowerMax, double ininBandNoise,double inBandwidthKhz, double inData \
                                       , vector<int>&inVecHeadName ,float** inGij,double* inNextNodePower,list<list<int> >*inListCluMem)
{

    chNum = inChNum;
    totalNodes = inTotalNodes;
    powerMax=inPowerMax*scale;
    inBandNoise = ininBandNoise*scale;
    bandwidthKhz = inBandwidthKhz;

    data_homo=inData;
    vecHeadName = &inVecHeadName;
    Gij =inGij;

    maIndexInterference = new int* [chNum];
    for (int i=0; i<chNum; i++)maIndexInterference[i] = new int [chNum];
    maStrengthInterference = new double* [chNum];
    for (int i=0; i<chNum; i++)maStrengthInterference[i] = new double [chNum];

    nextNodePower = inNextNodePower;
    listCluMember = inListCluMem;


    powerDifference = new double [totalNodes];
    for(int i=0; i<totalNodes; i++)powerDifference[i] = 0.0;
    powerDifferenceRatio = new double [totalNodes];
    for(int i=0; i<totalNodes; i++)powerDifferenceRatio[i] = 0.0;

}



ULConstraintSolver::~ULConstraintSolver()
{
    delete [] powerDifference;
    delete [] powerDifferenceRatio;

    for(int i=0; i<chNum; i++)delete[] maIndexInterference[i];
    delete [] maIndexInterference;
    for(int i=0; i<chNum; i++)delete[] maStrengthInterference[i];
    delete [] maStrengthInterference;
}

/*
    The main content of solve the existed Structure
    -Do power update
*/
bool ULConstraintSolver::solve()
{

    statusFlag=-1;
    chNum = listCluMember->size();
    for(int i=0; i<totalNodes; i++)powerDifference[i] = 0.0;
    for(int i=0; i<totalNodes; i++)powerDifferenceRatio[i] = 0.0;
    for(int i=0; i<totalNodes; i++)nextNodePower[i] =0;
    updateInterference();
    //------------------------------------//
    //add scale to avoid computation error//
    //------------------------------------//
    int loopCounter = 0;//counter to count how many round we updated Power
    exceedPc = false;
    while(loopCounter <2||(!exceedPc))
    {
        //cout<<"loop "<<loopCounter<<endl;
        //cout<<"avgR"<< avgRatio<<endl;
        if (loopCounter>2&&checkConverged())
        {
            //cout<<"Ratio Converged"<<endl;
            //break;
        }
        if (loopCounter>2&&checkDifference())
        {
            //cout<<"Difference is small"<<endl;
            break;
        }
        //Change the uplink node power cluster by cluster, which means that we will change all the members undet cluster i and go i+1
        updateInterference();// Need to update before change all member power to avoid the result from last time
        changeAllMemberPower();
        loopCounter++;
         //cout<<loopCounter<<"-th round"<<endl;

    }

    if (avgRatio >=1)//average ratio =1 means converged already
    {
        //cout<<"diverge avgRatio= " <<avgRatio<<endl;    //PowerDiverged
        for (int i=0; i<totalNodes; i++)
        {
            nextNodePower[i]/=scale;
            if (nextNodePower[i]>(powerMax+1e-6)){// cout<<nextNodePower[i]<<endl;
                exceedPc=true;
            }
            //cout<<"node "<<i <<": "<<nextNodePower[i]<<endl;
        //cout<<"  =>"<<nextNodePower[i]<<endl;
        }
        if(exceedPc){
            statusFlag=2;
        }
        else {statusFlag=3;}
        //cout<<"=========================="<<endl;
        return false;
    }

    for (int i=0; i<totalNodes; i++)
    {
        //cout<<"node i :"<<i <<" Ori= "<<nextNodePower[i];
        if(avgRatio<1){
          nextNodePower[i]=nextNodePower[i] + avgRatio * powerDifference[i] /(1-avgRatio);
          nextNodePower[i]/=scale;
        }
        //cout<<"node "<<i <<": "<<nextNodePower[i]<<endl;
        if (nextNodePower[i]>(powerMax+1e-6)){// cout<<nextNodePower[i]<<endl;
        exceedPc=true;}

        //cout<<"  =>"<<nextNodePower[i]<<endl;
    }
    //showVerificationResult();
    if (exceedPc == true)
    {
        statusFlag=2;
        //cout<<"exceedPc"<<endl;    // Infeasible
        //cout<<"================="<<endl;
        return false;
    }
    //cout<<loopCounter<<endl;
    //cout<<"vvvvvvvvvvvvvvvvvvvvvFeasiblevvvvvvvvvvvvvvv"<<endl;
    statusFlag=1;
    return true;
}

double ULConstraintSolver::solve_withT2Adj_BinerySearch(double inIniT2){
  double powerConverge=1e-6;
  double tmpPowerMax;
  double minFeasibleT2=-1;
  double maxInFeasiT2=0;
  int tmpNodeIndex=0;
  best2nd_ms=inIniT2;
  while(1){
    C2=data_homo/best2nd_ms/bandwidthKhz;
//    cout<<"C2 = "<<C2<<endl;
    if (solve()){

      returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax);
      if (tmpNodeIndex==-1){best2nd_ms=0;return best2nd_ms;}
//    cout<<"Ori Node: "<<tmpNodeIndex<<" ori power : "<<tmpPowerMax<<endl;
      if (abs(powerMax-tmpPowerMax)<powerConverge){
        //cout<<"power converge"<<endl;
            //cout<<"Time "<<(minFeasibleT2+maxInFeasiT2)/2<<endl;
        break;
      }
      else if((abs(minFeasibleT2-maxInFeasiT2)<1e-5)){
        //cout<<"Time converge"<<endl;
        int pi=-1;
        double pw=0;returnMaxNextNodePower(pi,pw);
           // cout<<"Time "<<(minFeasibleT2+maxInFeasiT2)/2<<endl;
        //cout<<"power max="<<pw<<endl;
        break;
      }
      minFeasibleT2=best2nd_ms;
      best2nd_ms=(best2nd_ms+maxInFeasiT2)/2;
    }
    else {

      returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax);
      maxInFeasiT2=best2nd_ms;
      if(minFeasibleT2==-1)best2nd_ms*=2;
      else{best2nd_ms=(best2nd_ms+minFeasibleT2)/2;}


//      cout<<x"Ori Node: "<<tmpNodeIndex<<" ori power : "<<tmpPowerMax<<endl;
      //cout<<"infeasible"<<endl;


//      cout<<endl;
    }
    //cout<<maxInFeasiT2<<" ~ "<<minFeasibleT2<<endl;
    //cout<<"next T2 is = "<<best2nd_ms<<endl;
   // cout<<"nodeIndex = "<<tmpNodeIndex<<"  nodePower"<<tmpPowerMax<<endl;
  }
  //showVerification
  criticalNode=tmpNodeIndex;
  criticalHeadIndex=returnHeadIndex_ByNodeName(criticalNode);


  return best2nd_ms;
}

double ULConstraintSolver::solve_withT2Adj_BinerySearch_2 (double inIniT2){
  double powerConverge=1e-6;
  double tmpPowerMax;
  bool firstFeasibleFound=false;
  double minT2=0;
  double maxT2=inIniT2;
  int tmpNodeIndex=0;
  double tmp=(minT2+maxT2)/2;
  double minFeasibleT=DBL_MAX;
  while(1){
    C2=data_homo/tmp/bandwidthKhz;
    //cout<<"C2 = "<<C2<<endl;
    bool feasFlag=solve();
    if (feasFlag){
      firstFeasibleFound=true;
      returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax);
      if(tmp<minFeasibleT)
        minFeasibleT=tmp;

      if (tmpNodeIndex==-1){best2nd_ms=0;return best2nd_ms;}
//    cout<<"Ori Node: "<<tmpNodeIndex<<" ori power : "<<tmpPowerMax<<endl;
      if (abs(powerMax-tmpPowerMax)<powerConverge){
        //cout<<"power converge"<<endl;
        //cout<<"Time "<<(minT2+maxT2)/2<<endl;
        break;
      }
      else if((abs(minT2-maxT2)<1e-5)){
        //cout<<"Time converge"<<endl;
        int pi=-1;
        double pw=0;returnMaxNextNodePower(pi,pw);
        //cout<<"Time "<<(minT2+maxT2)/2<<endl;
        //cout<<"power max="<<pw<<endl;
        break;
      }
      maxT2=tmp;
    }
    else {
      if(firstFeasibleFound){
         if(statusFlag==3){
            maxT2=tmp;
         }
         else{
            minT2=tmp;
         }
      }
      else{
          maxT2*=2;
      }

      if((abs(minT2-maxT2)<1e-6)){
        //Feasible solution lost because of precision
        tmp=minFeasibleT;
        C2=data_homo/tmp/bandwidthKhz;
        solve();
        break;
      }

    }
   returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax);
    //cout<<"T2="<<tmp<<" Power = "<<tmpPowerMax<<" Feasibility="<<feasFlag<<endl;
    //cout<<"             status "<<statusFlag<<endl;
    //cout<<"T_u="<<maxT2<<" T_l="<<minT2<<endl;

    tmp=(minT2+maxT2)/2;
  }
  //showVerification
  criticalNode=tmpNodeIndex;
  criticalHeadIndex=returnHeadIndex_ByNodeName(criticalNode);
  best2nd_ms=tmp;
  returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax);
  //cout<<"abcnxfdsfadfdf"<<endl;

  /*FILE *fid=fopen("ConvergePoweraGathered.txt","a+");
  fprintf(fid,"%f\n",tmpPowerMax);
  fclose(fid);
*/
  return best2nd_ms;
}














double ULConstraintSolver::solve_withT2Adj(double inIniT2){
  double powerConverge=1e-6;
  double tmpPowerMax;
  int tmpNodeIndex=0;
  best2nd_ms=inIniT2;

  while(1){
    cout<<"data "<<data_homo;
    C2=data_homo/best2nd_ms/bandwidthKhz;
    cout<<"C2 = "<<C2<<endl;
    if (solve()){
      updateInterference();
      if (abs(powerMax-tmpPowerMax)<powerConverge) break;
      returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax);
      cout<<"Ori Node: "<<tmpNodeIndex<<" ori power : "<<tmpPowerMax<<endl;
      updateT2_DataHomo(tmpNodeIndex,tmpPowerMax,true);

    }
    else {
        updateInterference();
        cout<<"infeasible"<<endl;
        returnMaxNextNodePower(tmpNodeIndex,tmpPowerMax);
        cout<<"Ori Node: "<<tmpNodeIndex<<" ori power : "<<tmpPowerMax<<endl;

        updateT2_DataHomo(tmpNodeIndex,tmpPowerMax,false);
    }

    cout<<"T2 is = "<<best2nd_ms<<endl;
    //cout<<"nodeIndex = "<<tmpNodeIndex<<"  nodePower"<<tmpPowerMax<<endl;

  }

  return best2nd_ms;
}
void ULConstraintSolver::returnMaxNextNodePower(int &nodeIndex, double &nodePower){
  nodePower=0;
  nodeIndex=-1;
  for(int i=0;i<totalNodes;i++){
    if (nextNodePower[i]>nodePower){
      nodePower=nextNodePower[i];
      nodeIndex=i;
    }
  }
  //assert(nodeIndex!=-1);
}
int ULConstraintSolver::returnHeadIndex_ByNodeName(int nodeName){
    list<list<int> > ::iterator it_LiInt=listCluMember->begin();
    for(int i=0; it_LiInt!=listCluMember->end();it_LiInt++,i++){
        list <int> ::iterator it_Int=it_LiInt->begin();
        for(;it_Int!=it_LiInt->end();it_Int++){
            if(*it_Int==nodeName){
                return i;
            }
        }
    }
    cout<<"Head Name Not Found in ULConstraintSolver!"<<endl;
    assert(1);
    return -1;
}

void ULConstraintSolver::updateT2_DataHomo(int &nodeIndex, double &nodePower, bool feasibility){
  list< list<int> >::iterator it_LiInt=listCluMember->begin();
  double ampConstant;
  if (feasibility)
    ampConstant=powerMax/nodePower;
  else{
    ampConstant=1;
  }
  double interferenceOfMaxPowerNode=0;
  for(int i=0;it_LiInt!=listCluMember->end();it_LiInt++,i++){
    list<int>::iterator it_Int=it_LiInt->begin();
    int size=it_LiInt->size();
    for (;it_Int!=it_LiInt->end();it_Int++){
      if(*it_Int==nodeIndex)break;
    }
    if(*it_Int==nodeIndex){
      for (int j=0;j<chNum;j++)interferenceOfMaxPowerNode+=maStrengthInterference[i][j];
      interferenceOfMaxPowerNode*=ampConstant;
       cout<<"max Interference: "<<interferenceOfMaxPowerNode<<endl;
      best2nd_ms=data_homo*(size-1)/(bandwidthKhz*log2(1+(powerMax*Gij[nodeIndex][(*vecHeadName)[i]]/(interferenceOfMaxPowerNode+inBandNoise))));

    }
  }
}
//solve with homogeneous power
void ULConstraintSolver::solve_withPowerHomo(double &totalT2){

    for(int i=0; i<totalNodes; i++)nextNodePower[i] =powerMax;
    for(int i=0; i<totalNodes;i++)(*nodes)[i].T2ndtier=0;
    updateInterference();
    totalT2=0;
    list <list <int> >::iterator it1 = listCluMember->begin();
    for(int headIndex=0; it1!=listCluMember->end(); headIndex++,it1++)
    {
        //Compute the Interference headIndex received.
        double accuInterference = inBandNoise;
        for(unsigned int i=0; i<(*vecHeadName).size(); i++)
        {
            if(maIndexInterference[headIndex][i]!=-1)
            {
                accuInterference+= maStrengthInterference[headIndex][i];
            }
        }
        //cout<<"Cluster: "<<headIndex<<" ->"<<accuInterference<<" C2 "<<C2<<endl;
        list <int>::iterator it2 =it1->begin();
        double T2ndtierOfSingleClu=0;
        for(;it2!=it1->end();it2++){
            (*nodes)[*it2].T2ndtier=data_homo/(bandwidthKhz*log2(1+powerMax*Gij[*it2][(*vecHeadName)[headIndex]]/(inBandNoise+accuInterference)));
            T2ndtierOfSingleClu+=(*nodes)[*it2].T2ndtier;
        }
        if(T2ndtierOfSingleClu>totalT2)totalT2=T2ndtierOfSingleClu;
    }
}


void ULConstraintSolver::Do1stTierPowerMax_ForBest(double &returnJoule, double &return1sttierMS, \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS,vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float  *Gib,vector<int>&vecBestHeadName ) {

    assert(bandwidthKhz!=-1);//if assrt fail, you use the wrong constructor or you forget to set bandwidth
    returnJoule=0;
    return1sttierMS=0;
    for (unsigned int i=0; i<vecBestHeadName.size(); i++) {
        vecClusterHeadMS[i]=(vecClusterHeadBits[i]/rateib_MaxPower[vecBestHeadName[i]])*1000;
        vecClusterHeadWatt[i]=powerMax;
        return1sttierMS+=vecClusterHeadMS[i];
        returnJoule+=(vecClusterHeadWatt[i]*vecClusterHeadMS[i]/1000);
        /*cout<<"SolverHERE "<<vecClusterHeadBits[i]<<"bits "<<rateib_MaxPower[vecBestHeadName[i]]<<"bps ,"\
        <<(vecClusterHeadBits[i]/rateib_MaxPower[vecBestHeadName[i]])*1000<<"ms"<<endl;*/
    }
}



void ULConstraintSolver::Do1stTierPowerControl(double &returnJoule, double &return1sttierMS, \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS,vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float  *Gib) {

    assert(bandwidthKhz!=-1);//if assrt fail, you use the wrong constructor
    need1stResorsIndividualMs =-1;
    //Find the longest resource needed
    for (unsigned int i=0; i<vecHeadName->size(); i++) {
        if ((1000*vecClusterHeadBits[i]/rateib_MaxPower[(*vecHeadName)[i]])>need1stResorsIndividualMs) {
            need1stResorsIndividualMs = (vecClusterHeadBits[i]/rateib_MaxPower[(*vecHeadName)[i]])*1000;
        }


        /*  cout<<"SolverHERE "<<vecClusterHeadBits[i]<<"bits "<<rateib_MaxPower[(*vecHeadName)[i]]<<"bps ,"\
          <<(vecClusterHeadBits[i]/rateib_MaxPower[(*vecHeadName)[i]]*1000)<<"ms"<<endl;*/

    }
    returnJoule = 0;
    return1sttierMS = need1stResorsIndividualMs * vecHeadName->size();
    for (unsigned int i=0; i<vecHeadName->size(); i++) {
        vecClusterHeadMS[i] = need1stResorsIndividualMs;
        vecClusterHeadWatt[i] = inBandNoise/scale *(pow(2,(vecClusterHeadBits[i]/(bandwidthKhz*need1stResorsIndividualMs)))-1) /Gib[(*vecHeadName)[i]];
        returnJoule+=(vecClusterHeadWatt[i]*need1stResorsIndividualMs)/1000;
    }
}


void ULConstraintSolver::Do1stTierPowerControl_ForBest(double &returnJoule, double &return1sttierMS, \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS,vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float  *Gib,vector<int>&vecBestHeadName ) {

    assert(bandwidthKhz!=-1);//if assrt fail, you use the wrong constructor
    double need1stResorsIndividualMs =-1;
    //Find the longest resource needed
    for (unsigned int i=0; i<vecBestHeadName.size(); i++) {
        if ((1000*vecClusterHeadBits[i]/rateib_MaxPower[vecBestHeadName[i]])>need1stResorsIndividualMs) {
            need1stResorsIndividualMs = (vecClusterHeadBits[i]/rateib_MaxPower[vecBestHeadName[i]])*1000;
        }
        /*cout<<"SolverHERE "<<vecClusterHeadBits[i]<<"bits "<<rateib_MaxPower[vecBestHeadName[i]]<<"bps ,"\
        <<(vecClusterHeadBits[i]/rateib_MaxPower[vecBestHeadName[i]])*1000<<"ms"<<endl;*/

    }
    returnJoule = 0;
    return1sttierMS = need1stResorsIndividualMs * vecBestHeadName.size();
    for (unsigned int i=0; i<vecBestHeadName.size(); i++) {
        vecClusterHeadMS[i] = need1stResorsIndividualMs;
        vecClusterHeadWatt[i] = inBandNoise/scale *(pow(2,(vecClusterHeadBits[i]/(bandwidthKhz*need1stResorsIndividualMs)))-1) /Gib[vecBestHeadName[i]];
        returnJoule+=(vecClusterHeadWatt[i]*need1stResorsIndividualMs/1000);
    }
}





void ULConstraintSolver::Do1stTierPowerControl_ForBest_PibAdjust(double &returnJoule, double &return1sttierMS, \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS,vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float  *Gib,vector<int>&vecBestHeadName,double pibAdjust_Scale ) {
    assert(bandwidthKhz!=-1);//if assrt fail, you use the wrong constructor
    double need1stResorsIndividualMs =-1;
    //Find the longest resource needed
    for (unsigned int i=0; i<vecBestHeadName.size(); i++) {
        if ((1000*vecClusterHeadBits[i]/rateib_MaxPower[vecBestHeadName[i]])>need1stResorsIndividualMs) {
            need1stResorsIndividualMs = (vecClusterHeadBits[i]/rateib_MaxPower[vecBestHeadName[i]])*1000;
        }
        /*cout<<"SolverHERE "<<vecClusterHeadBits[i]<<"bits "<<rateib_MaxPower[vecBestHeadName[i]]<<"bps ,"\
        <<(vecClusterHeadBits[i]/rateib_MaxPower[vecBestHeadName[i]])*1000<<"ms"<<endl;*/

    }
    returnJoule = 0;
    return1sttierMS = need1stResorsIndividualMs * vecBestHeadName.size();
    for (unsigned int i=0; i<vecBestHeadName.size(); i++) {
        vecClusterHeadMS[i] = need1stResorsIndividualMs;
        vecClusterHeadWatt[i] =pibAdjust_Scale* inBandNoise/scale *(pow(2,(vecClusterHeadBits[i]/(bandwidthKhz*need1stResorsIndividualMs)))-1) /Gib[vecBestHeadName[i]];
        returnJoule+=(vecClusterHeadWatt[i]*need1stResorsIndividualMs/1000);
    }
}




bool ULConstraintSolver::Do1stTierPowerControl_1stLimited_EqDiv(double &returnJoule,  \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS,vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float  *Gib,double timeMs1st) {

    assert(bandwidthKhz!=-1);//if assrt fail, you use the wrong constructor
    bool powerFeasible=true;
    int headNum = vecHeadName->size();
    need1stResorsIndividualMs =timeMs1st/static_cast<double>(headNum);

    returnJoule = 0;
    for (unsigned int i=0; i<vecHeadName->size(); i++) {
        vecClusterHeadMS[i] = need1stResorsIndividualMs;
        vecClusterHeadWatt[i] = inBandNoise/scale *(pow(2,(vecClusterHeadBits[i]/(bandwidthKhz*need1stResorsIndividualMs)))-1) /Gib[(*vecHeadName)[i]];
        if( vecClusterHeadWatt[i]>powerMax)powerFeasible=false;
        returnJoule+=(vecClusterHeadWatt[i]*need1stResorsIndividualMs)/1000;
    }
    return powerFeasible;
}


bool ULConstraintSolver::Do1stTierPowerControl_ForBest_1stLimited_EqDiv(double &returnJoule,  \
        vector<double>&vecClusterHeadWatt, vector<double> &vecClusterHeadMS,vector<double> &vecClusterHeadBits, \
        double *rateib_MaxPower,float  *Gib,vector<int>&vecBestHeadName,double timeMs1st ) {

    assert(bandwidthKhz!=-1);//if assrt fail, you use the wrong constructor
    bool powerFeasible=true;
    int headNum = vecHeadName->size();
    double need1stResorsIndividualMs =timeMs1st/static_cast<double>(headNum);
    returnJoule = 0;
    for (unsigned int i=0; i<vecBestHeadName.size(); i++) {
        vecClusterHeadMS[i] = need1stResorsIndividualMs;
        vecClusterHeadWatt[i] = inBandNoise/scale *(pow(2,(vecClusterHeadBits[i]/(bandwidthKhz*need1stResorsIndividualMs)))-1) /Gib[vecBestHeadName[i]];
        if( vecClusterHeadWatt[i]>powerMax)powerFeasible=false;
        returnJoule+=(vecClusterHeadWatt[i]*need1stResorsIndividualMs/1000);
    }
    return powerFeasible;
}








void ULConstraintSolver::updateInterference()
{
    list<list <int> >::iterator itl = listCluMember->begin();
    for(int i=0; i<chNum; i++,itl++) //For each CH we will see there member 1-by-1
    {
        if((*vecHeadName)[i]==-1){continue;}//Cluster Not exist
        list<list <int> >:: iterator itl2 = listCluMember->begin(); //search other cluster from start for each CH
        for (int j=0; j<chNum; j++,itl2++) //Search other cluster for each member in cluster (*vecHeadName)[i]
        {
            int interferentSource = -1;
            float maxInterference = -99;
            if (i==j)//Same headIndex which means this is self head.
            {
                maIndexInterference[i][j]=-1;
                maStrengthInterference[i][j]=0;
                continue;

            }
            else if (itl2->size()<=1)//empty list ->0
            {
                maIndexInterference[i][j] = -1;
                maStrengthInterference[i][j] = 0;
                continue;
            }
            //Find the strongest interference from other group
            list<int> ::iterator itl3 = itl2->begin();
            for(; itl3!=itl2->end(); itl3++) //search every node(itl3) under each cluster(itl2) != self cluster
            {

                float tempInterference = nextNodePower[(*itl3)] * Gij[(*vecHeadName)[i]][(*itl3)];
                if((*vecHeadName)[j]==(*itl3))continue;//In Downlink Scenario othe head won't transmit power
                else if (tempInterference>maxInterference )
                {
                    interferentSource = (*itl3);
                    maxInterference = tempInterference;
                }
            }
            maIndexInterference[i][j] = interferentSource;
            maStrengthInterference[i][j] = maxInterference;//Interference form j(HeadIndex) to i(HeadIndex)
        }
    }
}

/*
   Update the members power according to the power last round
*/
void ULConstraintSolver::changeAllMemberPower()
{
    list <list <int> >::iterator it1 = listCluMember->begin();
    for(int headIndex=0; it1!=listCluMember->end(); headIndex++,it1++)
    {
        //Compute the Interference headIndex received.
        if ((*vecHeadName)[headIndex]==-1)continue;

        double accuInterference = inBandNoise;
        for(unsigned int i=0; i<(*vecHeadName).size(); i++)
        {
            if ((*vecHeadName)[i]==-1)continue;
            if(maIndexInterference[headIndex][i]!=-1)
            {
                accuInterference+= maStrengthInterference[headIndex][i];
            }
        }
        //cout<<"Cluster: "<<headIndex<<" ->"<<accuInterference<<" C2 "<<C2<<endl;
        list <int>::iterator it2 =it1->begin();
        // update all the power  in the cluster we are interested.
        int sizeOfCluter = it1->size();
        for(int i=0; it2!=it1->end(); i++,it2++)
        {
            if(sizeOfCluter==1||(*vecHeadName)[headIndex]==(*it2))//if the cluster has only one head or the *it2 is the head(same Name)
            {
                powerDifference[(*it2)]=0;
                powerDifferenceRatio[(*it2)] =0;
                nextNodePower[(*it2)] = 0;
                continue;
            }

            else
            {
                float powerCursor = 0;
                double tempDifference = 0;
                cout<<"c2 "<<C2*(sizeOfCluter-1)<<endl;
                cout<<"Channel Gain from"<<(*vecHeadName)[headIndex]<<" "<<(*it2)<<" "<<Gij[(*vecHeadName)[headIndex]][(*it2)];
                powerCursor = (accuInterference *(pow((float)2.0, (sizeOfCluter-1)*C2)-1))/Gij[(*vecHeadName)[headIndex]][(*it2)];
                cout<<"Power "<<powerCursor/scale<< " accu "<< accuInterference/scale<<endl;
                tempDifference = powerCursor -nextNodePower[(*it2)];

                if(powerDifference[(*it2)]!=0)powerDifferenceRatio[(*it2)] = (double)tempDifference / powerDifference[(*it2)];
                powerDifference[(*it2)] = tempDifference;
                nextNodePower[(*it2)] = powerCursor;
                /*
                cout<< "update power "<<*it2<<" "<<nextNodePower[(*it2)]<<" and P*G"<< Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)]  <<endl;
                cout.precision(5);
                cout<<scientific<<"Node "<<*it2<<" SINR "<< Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)] /accuInterference<< \
                 "  Rate:"<<log2(1+Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)] /accuInterference)<<endl;
                */

                if(nextNodePower[(*it2)]>powerMax)
                {
                    exceedPc=true;
                }
            }
            //cout<<endl;
        }
    }
}



/*
     Check whether if the ratio is converged don't check with avg check with random 1;
*/
 bool ULConstraintSolver::checkConverged()
{
    bool convergence = true;
    double thre = 1e-5;
    avgRatio=0; // find the avgRatio
    double avgMember =0;

    //compute average power ratio
    for (int i =0; i<totalNodes; i++)
    {
        if (powerDifferenceRatio[i]!=0)
        {
            avgRatio+=powerDifferenceRatio[i];
            avgMember++;
        }
    }
    if(avgMember!=0) avgRatio/=avgMember;
    //cout<<"AvGr"<<avgRatio<<endl;
    //check if all the power ratio close to the average

    for (int i =0; i<chNum; i++)
    {
        assert(avgRatio>-1);
        if (abs(powerDifferenceRatio[i]-avgRatio)>thre&&powerDifferenceRatio[i]!=0)
        {
            convergence = false;
            break;
        }
    }
    return convergence;
}

/*
   Check if the power difference is converged.

*/
bool ULConstraintSolver::checkDifference()  //check absolute difference
{
    double thre = 1e-12;
    bool allTooSmall = true;
    for(int i=0; i<totalNodes; i++)
    {
        //cout<<"power Difference "<<powerDifference[i]<<endl;
        if(abs(powerDifference[i]) > thre &&powerDifference[i]!=0)
        {
            allTooSmall = false;
            break;
        }
    }
    return allTooSmall;
} //check if difference is too small


/*
  Use if the structure change lead to infeasible.

  -returen the infesibility of the "whole system" to compare



*/
double ULConstraintSolver::returnInfeasibility()
{
    avgRatio=0; // find the avgRatio
    double avgMember =0;

    //compute average power ratio
    for (int i =0; i<totalNodes; i++)
    {
        if (powerDifferenceRatio[i]!=0)
        {
            avgRatio+=powerDifferenceRatio[i];
            avgMember++;
        }
    }
    if(avgMember==0) avgRatio = 0;
    else avgRatio/=avgMember;
    //cout<<"avgRatio:"<<avgRatio<<endl;
    if (avgRatio<1)return 0.0;
    else return 3*log(avgRatio);
}





void ULConstraintSolver::printPowerRAll()
{
    for(int i=0; i<totalNodes; i++)
    {
        cout<<i<<" "<<powerDifferenceRatio[i]<<endl;
    }
}

void ULConstraintSolver::printPowerAll()
{
    for(int i=0; i<totalNodes; i++)
    {
        cout<<i<<" "<<nextNodePower[i]<<endl;
    }
}

void ULConstraintSolver::printPowerDifAll()
{
    for(int i=0; i<totalNodes; i++)
    {
        cout<<i<<" "<<powerDifference[i]<<endl;
    }
}
void ULConstraintSolver::showVerificationResult(vector <double> &vecBestReceivedInterference,vector<double> &vecBestSINR_forVerification, \
                                    vector<double> &vecBestBpsHz_forVerification)
{
    vecBestReceivedInterference.clear();
    updateInterference();
    list <list <int> >::iterator it1 = listCluMember->begin();
    for(int headIndex=0; it1!=listCluMember->end(); headIndex++,it1++)
    {
        //Compute the Interference headIndex received.
        if ((*vecHeadName)[headIndex]==-1)continue;
        double accuInterference = inBandNoise;
        for(int i=0; i<chNum; i++)
        {
            if(maIndexInterference[headIndex][i]!=-1)
            {
                accuInterference+= maStrengthInterference[headIndex][i];
            }
        }
        vecBestReceivedInterference.push_back(accuInterference-inBandNoise);
        //cout<<"Cluster: "<<headIndex<<" BpsHz "<<C2<<endl;
        list <int>::iterator it2 =it1->begin();
        // update all the power  in the cluster we are interested.
        int sizeOfCluter = it1->size();
        for(int i=0; it2!=it1->end(); i++,it2++)
        {
            vecBestSINR_forVerification[*it2]=(Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)] /accuInterference);
            vecBestBpsHz_forVerification[*it2]=(log2(1+Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)] /accuInterference)/(sizeOfCluter-1));
            /*
            cout.precision(5);
            cout<<scientific<<"Normal: Node "<<*it2<<" SINR "<< Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)] /accuInterference<< \
            "  BpsHz:"<<log2(1+Gij[(*vecHeadName)[headIndex]][(*it2)]*nextNodePower[(*it2)] /accuInterference)/(sizeOfCluter-1)<<endl;
            */
        }
//        cout<<"---------------"<<endl;
    }
}

int ULConstraintSolver::returnCriticalNode(){
    return criticalNode;
}
int ULConstraintSolver::returnCriticalHeadIndex(){
    return criticalHeadIndex;
}

void ULConstraintSolver::setinBand(double in){
    bandwidthKhz=in;
}

//-------------------------------------------------------------//
/*

    After change the power of each cluster member for one cluster
    update the interference list for that cluster


void ULConstraintSolver::updateSingleClusterInterference(int headIndex,list <list<int> >::iterator it1)
{
    //update the interferece
    for (int i=0; i<chNum; i++)
    {
        if(headIndex==i)continue;
        list <int> ::iterator it2 = it1->begin();
        for(; it2!=it1->end(); it2++)
        {
            if (maStrengthInterference[i][headIndex]<(Gij[(*it2)][(*vecHeadName)[i]]* nextNodePower[*it2]))
            {
                maStrengthInterference[i][headIndex] = (Gij[(*it2)][(*vecHeadName)[i]]* nextNodePower[*it2]);
                maIndexInterference[i][headIndex] = *it2;
            }
        }
    }
}


    Change the power of each cluster member for cluster( headIndex)

void ULConstraintSolver::changeSingleClusterMemberPower(int headIndex,list <list<int> >::iterator it1)
{
    double accuInterference=inBandNoise;
    //compute all the interference for one head
    for(int i=0; i<chNum; i++)
    {
        if(maStrengthInterference[headIndex][i]!=-1)
            accuInterference+= maStrengthInterference[headIndex][i];
    }
    list <int>::iterator it2 =it1->begin();
    int sizeOfCluter = it1->size();
    // update all the power  of the cluster we are interested.)
    for(int i=0; i<sizeOfCluter; i++,it2++)
    {
        if(sizeOfCluter==1||(*vecHeadName)[headIndex]==(*it2))//special case of power update
        {
            powerDifference[(*it2)]=0;
            powerDifferenceRatio[(*it2)] =0;
            nextNodePower[(*it2)] = 0;
            continue;
        }
        else
        {
            float powerCursor = 0;
            double tempDifference = 0;
            powerCursor = (accuInterference *(pow((float)2.0, (sizeOfCluter-1)*C2)-1))/Gij[(*vecHeadName)[headIndex]][(*it2)];
            tempDifference = powerCursor -nextNodePower[(*it2)];

            if(powerDifference[(*it2)]!=0)powerDifferenceRatio[(*it2)] = (double)tempDifference / powerDifference[(*it2)];
            powerDifference[(*it2)] = tempDifference;
            nextNodePower[i] = powerCursor;


            if(nextNodePower[(*it2)]>powerMax)
            {
                exceedPc=true;
            }
        }
    }
}
*/





