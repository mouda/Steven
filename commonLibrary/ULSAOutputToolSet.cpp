#include<iostream>
#include<cstdio>
#include<list>
#include<vector>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <linux/netdevice.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <unistd.h>

using namespace std;
#include"ULSAOutputToolSet.h"


/*
template<class  T> ULSAOutputToolSet<T>::ULSAOutputToolSet(T *inMyobj){
  myObj = *inMyobj;// Programming like these need to construct extra operator overload and it will construct on more object
  //Need to be further discussed
}*/
template<class  T> void ULSAOutputToolSet<T>::writeIniHead(char* filename, T &myObj)
{

    FILE* fid;
    fid = fopen(filename,"w");
    fprintf(fid,"%zu\n",myObj.cSystem->vecHeadName.size());

    for(unsigned int i=0; i<myObj.cSystem->vecHeadName.size(); i++)
    {
        fprintf(fid,"%zu\n",myObj.cSystem->vecHeadName[i]);
    }
    fclose(fid);
}


template<class  T> void ULSAOutputToolSet<T>::writeBestStru(char *filename,T &myObj)
{
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=myObj.best2nd_ms+myObj.best1st_ms;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myObj.listCluMemBest->begin();
    for(; it1!=myObj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        for(; it2!=it1->end(); it2++)
        {
            if(myObj.powerBest[*it2]>0) {
                bestTotalJoule+=(myObj.powerBest[*it2]*myObj.best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myObj.powerBest[*it2]*myObj.best2nd_ms/(tempSize-1))/1000;
            }
        }
    }
    bestTotalJoule +=myObj.best1st_Joule;



    FILE *fid;
    fid = fopen(filename,"a");
    cout<<filename<<endl;
    fprintf(fid,"%d\n",myObj.totalNodes);
    fprintf(fid,"%d\n",myObj.maxChNum);
    fprintf(fid,"%e\n",myObj.powerMax);
    fprintf(fid,"%e\n",myObj.C2);
    fprintf(fid,"%d\n", myObj.bestFeasibleSupNum-myObj.maxChNum);
    fprintf(fid,"%d\n", myObj.SAIter);
    //fprintf(fid,"%f\n", myObj.computingTimes);
    fprintf(fid,"%5e\n", myObj.best1st_ms);
    fprintf(fid,"%5e\n", myObj.best2nd_ms);
    fprintf(fid,"%5e\n", myObj.best1st_Joule);
    fprintf(fid,"%5e\n", best2ndJoule);


    for (int i= 0; i<myObj.maxChNum; i++)
    {
        for( int j=0; j<myObj.totalNodes; j++)
        {
            fprintf(fid,"%d ",myObj.bestMaClusterStru[i][j]);
        }
        fprintf(fid,"\n");
    }

    for (int i= 0; i<myObj.maxChNum; i++)
    {
        fprintf(fid,"%d ",myObj.vecHeadNameBest[i]+1);
    }
    fprintf(fid,"\n");
    for (int i= 0; i<myObj.totalNodes; i++)
    {
        fprintf(fid,"%E ", myObj.powerBest[i]);
    }
    fprintf(fid,"\n");
    it1=myObj.listCluMemBest->begin();
    for (; it1!=myObj.listCluMemBest->end(); it1++ )
    {
        list <int>::iterator it2 = it1->begin();
        for (; it2!=it1->end(); it2++)
            fprintf(fid,"%d ",(*it2)+1);
        fprintf(fid,"\n");
    }
    fprintf(fid,"\n");
    fprintf(fid,"============Find Best @ Round %d================\n",myObj.roundBest);
    list <list <int> >::iterator it11 = myObj.listCluMemBest->begin();
    for(int headIndex=0; it11!=myObj.listCluMemBest->end(); headIndex++,it11++)
    {
        //Compute the Interference headIndex received.
        fprintf(fid,"==========Cluster %d, RcvInterf=%.5e============\n",headIndex,myObj.vecBestReceivedInterference[headIndex]);
        list <int>::iterator it22 =it11->begin();
        for(int i=0; it22!=it11->end(); i++,it22++)
        {
            fprintf(fid,"Node %d, SINR:%.5e, BpsHz:%.5e\n", *it22,myObj.vecBestSINR_forVerification[*it22],myObj.vecBestBpshz_forVerification[*it22]);
        }
        fprintf(fid,"---------------\n");
    }






    fclose(fid);
}


template<class  T> void ULSAOutputToolSet<T>::writeBestStru_V2(char *filename,T &myObj)
{
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=myObj.best2nd_ms+myObj.best1st_ms;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myObj.listCluMemBest->begin();
    for(; it1!=myObj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        if (tempSize==1)continue;
        for(; it2!=it1->end(); it2++)
        {
            if(myObj.powerBest[*it2]>0) {
                bestTotalJoule+=(myObj.powerBest[*it2]*myObj.best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myObj.powerBest[*it2]*myObj.best2nd_ms/(tempSize-1))/1000;
            }
        }
    }
    bestTotalJoule +=myObj.best1st_Joule;



    FILE *fid;
    fid = fopen(filename,"a");
    fprintf(fid,"%d\n",myObj.totalNodes);
    fprintf(fid,"%d\n",myObj.maxChNum);
    fprintf(fid,"%e\n",myObj.powerMax);
    fprintf(fid,"%e\n",myObj.indEntropy);
    fprintf(fid,"%d\n", myObj.bestFeasibleSupNum-myObj.maxChNum);
    fprintf(fid,"%d\n", myObj.SAIter);
    //fprintf(fid,"%f\n", myObj.computingTimes);
    fprintf(fid,"%5e\n", myObj.best1st_ms);
    fprintf(fid,"%5e\n", myObj.best2nd_ms);
    fprintf(fid,"%5e\n", myObj.best1st_Joule);
    fprintf(fid,"%5e\n", best2ndJoule);


    for (int i= 0; i<myObj.maxChNum; i++)
    {
        for( int j=0; j<myObj.totalNodes; j++)
        {
            fprintf(fid,"%d ",myObj.bestMaClusterStru[i][j]);
        }
        fprintf(fid,"\n");
    }

    for (int i= 0; i<myObj.maxChNum; i++)
    {
        fprintf(fid,"%d ",myObj.vecHeadNameBest[i]+1);
    }
    fprintf(fid,"\n");
    for (int i= 0; i<myObj.totalNodes; i++)
    {
        fprintf(fid,"%E ", myObj.powerBest[i]);
    }
    fprintf(fid,"\n");
    it1=myObj.listCluMemBest->begin();
    for (; it1!=myObj.listCluMemBest->end(); it1++ )
    {
        list <int>::iterator it2 = it1->begin();
        for (; it2!=it1->end(); it2++)
            fprintf(fid,"%d ",(*it2)+1);
        fprintf(fid,"\n");
    }
    fprintf(fid,"\n");
    fprintf(fid,"============Find Best @ Round %d================\n",myObj.roundBest);
    list <list <int> >::iterator it11 = myObj.listCluMemBest->begin();
    for(int headIndex=0; it11!=myObj.listCluMemBest->end(); headIndex++,it11++)
    {
        //Compute the Interference headIndex received.
        fprintf(fid,"==========Cluster %d, RcvInterf=%.5e============\n",headIndex,myObj.vecBestReceivedInterference[headIndex]);
        list <int>::iterator it22 =it11->begin();
        for(int i=0; it22!=it11->end(); i++,it22++)
        {
            fprintf(fid,"Node %d, SINR:%.5e, BpsHz:%.5e\n", *it22,myObj.vecBestSINR_forVerification[*it22],myObj.vecBestBpshz_forVerification[*it22]);
        }
        fprintf(fid,"---------------\n");
    }






    fclose(fid);
}


template<class  T> void ULSAOutputToolSet<T>::writeAll(char* filename,double density, T& myobj)
{
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=myobj.best2nd_ms+myobj.best1st_ms;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myobj.listCluMemBest->begin();
    for(; it1!=myobj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        if (tempSize==1)continue;
        for(; it2!=it1->end(); it2++)
        {
            if(myobj.powerBest[*it2]>0) {
                bestTotalJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
            }
        }
    }
    bestTotalJoule +=myobj.best1st_Joule;

    FILE *fid;
    cout<< filename<<endl;
    fid = fopen(filename,"a");

    // cout format:
    // first two rows are given parameters are constant
    // totalNodes, maxChNum, wholeSystemEntropy, CompressionRatiom
    // fidelity Ratio = -1 , in ULSA3/2g ; Given, Otherwise

    //            1  2    3    4   5  6   7   8           9 10  11  12  13  14       15  16
    fprintf(fid,"%d %d %.2f %.2f %2f %d %.1f %d        %.3f %d %5e %5e %5e %5e      %5e %5e\n",  \
            myobj.totalNodes,myobj.maxChNum, myobj.wholeSystemEntopy, myobj.returnComprRatio(),\
            myobj.fidelityRatio, myobj.quantizationBits, density, myobj.roundBest, \
             myobj.bestFeasibleJEntropy, myobj.bestFeasibleSupNum, \
            myobj.best1st_Joule, best2ndJoule, myobj.best1st_ms, myobj.best2nd_ms, \
             bestTotalJoule, bestTotal_ms);
    fclose(fid);

}



template<class  T> void ULSAOutputToolSet<T>::writeClusterInfo(char* filename, T& myobj,char *time)
{
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=myobj.best2nd_ms+myobj.best1st_ms;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myobj.listCluMemBest->begin();

    vector<double> clusterSecondTier;
    for(; it1!=myobj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        if (tempSize==1)continue;
        double tempEnergy=0;
        for(; it2!=it1->end(); it2++)
        {
            if(myobj.powerBest[*it2]>0) {
                bestTotalJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
                tempEnergy+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
            }
        }
        clusterSecondTier.push_back(tempEnergy);
    }
    bestTotalJoule +=myobj.best1st_Joule;

    FILE *fid;
    cout<< filename<<endl;
    fid = fopen(filename,"a");

    for(unsigned int i=0;i<myobj.vecHeadNameBest.size();i++){

    //                1  2  3     4    5    6       7   8     9      9 10  11  12  13  14       15  16
        fprintf(fid,"%s  %d %d %.2e    %.1f %.2e %.1e %.1e %.2f \n", \
            time, myobj.vecHeadNameBest[i]+1, \
            myobj.vecBestClusterSize[i], myobj.vecBestClusterBits[i]/(myobj.vecBestClusterSize[i]*myobj.indEntropy), \
            myobj.vecBestClusterBits[i],myobj.vecBestClusterHeadWatt[i]*myobj.vecBestClusterHeadMS[i]/1000, \
            clusterSecondTier[i],  \
            myobj.vecBestClusterHeadWatt[i], myobj.vecBestClusterHeadMS[i]);

    }
    fprintf(fid,"\n");
    fclose(fid);

}



template<class T> void ULSAOutputToolSet<T>::summaryNwrite2tiers(char* filename, T& myobj) {
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=myobj.best2nd_ms+myobj.best1st_ms;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myobj.listCluMemBest->begin();
    for(; it1!=myobj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        for(; it2!=it1->end(); it2++)
        {
            if(myobj.powerBest[*it2]>0) {
                bestTotalJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
            }
        }
    }
    bestTotalJoule +=myobj.best1st_Joule;
    FILE *fid=fopen(filename,"a");
    fprintf(fid,"------------------\n");
    fprintf(fid,"SupNum %3d       Sup Ratio: %f\n" , myobj.bestFeasibleSupNum, myobj.bestFeasibleSupNum/static_cast<double>(myobj.totalNodes));
    fprintf(fid,"Info   %5f     Sup Info Ratio: %5e CR:%.5e\n" , myobj.bestFeasibleJEntropy, myobj.bestFeasibleJEntropy/myobj.wholeSystemEntopy,myobj.returnComprRatio());
    fprintf(fid,"total: %f(ms),     1st:%5e(ms)     2nd:%5e(ms)\n" ,bestTotal_ms,myobj.best1st_ms,myobj.best2nd_ms);
    fprintf(fid,"total: %e(joule),  1st:%5e(joule)  2nd:%5e(joule)\n" ,  bestTotalJoule, myobj.best1st_Joule,best2ndJoule);


    for(int i=0; i<myobj.maxChNum; i++)
    {
        fprintf(fid,"%d-th Head %d, DataLoad:%f(bits), Power: %5e(Watt), Time: %5e(ms) Energy: %5e(Joule) \n", \
                i ,myobj.vecHeadNameBest[i] , myobj.vecBestClusterBits[i],myobj.vecBestClusterHeadWatt[i], myobj.vecBestClusterHeadMS[i], \
                myobj.vecBestClusterHeadWatt[i]*myobj.vecBestClusterHeadMS[i]/1000);
    }
    fclose(fid);


}

template<class T> void ULSAOutputToolSet<T>::writeEnergy(char *filename,T &myobj) {
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=myobj.best2nd_ms+myobj.best1st_ms;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myobj.listCluMemBest->begin();
    for(; it1!=myobj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        if (tempSize==1)continue;
        for(; it2!=it1->end(); it2++)
        {
            if(myobj.powerBest[*it2]>0) {
                bestTotalJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
            }
        }
    }
    bestTotalJoule +=myobj.best1st_Joule;
    FILE *fid=fopen(filename,"a");
    fprintf(fid,"%f %f %f %f %d %f\n",bestTotalJoule, myobj.best1st_Joule,best2ndJoule,bestTotal_ms, myobj.maxChNum,myobj.returnComprRatio());

    fclose(fid);

}
template<class T> void ULSAOutputToolSet<T>::summaryNwrite2tiers_MinResors_with2ndPowerControl(char* filename, T& myobj, double best2nd_ms) {
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=best2nd_ms+myobj.best1st_ms;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myobj.listCluMemBest->begin();
    for(; it1!=myobj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        if (tempSize==1)continue;
        for(; it2!=it1->end(); it2++)
        {
            if(myobj.powerBest[*it2]>0) {
                bestTotalJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myobj.powerBest[*it2]*myobj.best2nd_ms/(tempSize-1))/1000;
            }
        }
    }
    bestTotalJoule +=myobj.best1st_Joule;
    FILE *fid=fopen(filename,"a");
    fprintf(fid,"------------------\n");
    fprintf(fid,"SupNum %3d       Sup Ratio: %f\n" , myobj.bestFeasibleSupNum, myobj.bestFeasibleSupNum/static_cast<double>(myobj.totalNodes));
    fprintf(fid,"Info   %5f     Sup Info Ratio: %5e CR:%.5e\n" , myobj.bestFeasibleJEntropy, myobj.bestFeasibleJEntropy/myobj.wholeSystemEntopy,myobj.returnComprRatio());
    fprintf(fid,"total: %f(ms),     1st:%5e(ms)     2nd:%5e(ms)\n" ,bestTotal_ms,myobj.best1st_ms,best2nd_ms);
    fprintf(fid,"total: %e(joule),  1st:%5e(joule)  2nd:%5e(joule)\n" ,  bestTotalJoule, myobj.best1st_Joule,best2ndJoule);


    for(int i=0; i<myobj.maxChNum; i++)
    {
        fprintf(fid,"%d-th Head %d, DataLoad:%f(bits), Power: %5e(Watt), Time: %5e(ms) Energy: %5e(Joule) \n", \
                i ,myobj.vecHeadNameBest[i] , myobj.vecBestClusterBits[i],myobj.vecBestClusterHeadWatt[i], myobj.vecBestClusterHeadMS[i], \
                myobj.vecBestClusterHeadWatt[i]*myobj.vecBestClusterHeadMS[i]/1000);
    }
    fclose(fid);


}



template<class T> void ULSAOutputToolSet<T>::writePeformance_MinResors_with2ndPowerControl(char *filename,T &myobj, double best2nd_ms, double fidelityRatio) {
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=best2nd_ms+myobj.best1st_ms;
    double density = -1;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myobj.listCluMemBest->begin();
    for(; it1!=myobj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        if (tempSize==1)continue;
        for(; it2!=it1->end(); it2++)
        {
            if(myobj.powerBest[*it2]>0) {
                bestTotalJoule+=(myobj.powerBest[*it2]*best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myobj.powerBest[*it2]*best2nd_ms/(tempSize-1))/1000;
            }
        }
    }
    bestTotalJoule +=myobj.best1st_Joule;
    FILE *fid=fopen(filename,"a");
    fprintf(fid,"%d %d %.2f %.2f %2f %d %.1f %d        %.3f %d %5e %5e %5e %5e      %5e %5e\n",  \
            myobj.totalNodes,myobj.maxChNum, myobj.wholeSystemEntopy, myobj.returnComprRatio(),\
            myobj.fidelityRatio, myobj.quantizationBits, density, myobj.roundBest, \
             myobj.bestFeasibleJEntropy, myobj.bestFeasibleSupNum, \
            myobj.best1st_Joule, best2ndJoule, myobj.best1st_ms, myobj.best2nd_ms, \
             bestTotalJoule, bestTotal_ms);

    fclose(fid);

}



template<class T> void ULSAOutputToolSet<T>::writePeformance_MinResors_with2ndPowerControl_4b(char *filename,T &myobj, double best2nd_ms, double fidelityRatio,int bestChNum) {
    bestTotalJoule=0;
    best2ndJoule=0;
    bestTotal_ms=best2nd_ms+myobj.best1st_ms;
    double density = -1;
    //Calculate 2nd-tier first
    list<list<int> >::iterator it1=myobj.listCluMemBest->begin();
    for(; it1!=myobj.listCluMemBest->end(); it1++)
    {
        list<int>::iterator it2=it1->begin();
        double tempSize = static_cast<double>(it1->size());
        if (tempSize==1)continue;
        for(; it2!=it1->end(); it2++)
        {
            if(myobj.powerBest[*it2]>0) {
                bestTotalJoule+=(myobj.powerBest[*it2]*best2nd_ms/(tempSize-1))/1000;
                best2ndJoule+=(myobj.powerBest[*it2]*best2nd_ms/(tempSize-1))/1000;
            }
        }
    }
    bestTotalJoule +=myobj.best1st_Joule;
    FILE *fid=fopen(filename,"a");
    fprintf(fid,"%d %d %.2f %.2f %2f %d %.1f %d        %.3f %d %5e %5e %5e %5e %d     %5e %5e\n",  \
            myobj.totalNodes,myobj.maxChNum, myobj.wholeSystemEntopy, myobj.returnComprRatio(),\
            myobj.fidelityRatio, myobj.quantizationBits, density, myobj.roundBest, \
             myobj.bestFeasibleJEntropy, myobj.bestFeasibleSupNum, \
            myobj.best1st_Joule, best2ndJoule, myobj.best1st_ms, myobj.best2nd_ms, \
            bestChNum, \
             bestTotalJoule, bestTotal_ms);

    fclose(fid);

}


template<class T> void ULSAOutputToolSet<T>::showIndEntropy(T &myobj){
    cout<<"In ULSAOutput indEntropy="<<myobj.indEntropy<<endl;
}

template<class T> std::string ULSAOutputToolSet<T>::getIpAddr() {

  int s;
  struct ifconf ifconf;
  struct ifreq ifr[50];
  int ifs;
  int i;
  int domain =  AF_INET;
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
  close(s);
  return strIpAddr.substr(dlim.length(),strIpAddr.length()) ;

}
//template class ULSAOutputToolSet<class ULSA2k_MC>;
//template class ULSAOutputToolSet<class ULSA2j_MC>;
//template class ULSAOutputToolSet<class T2_Change_Developing>;
//template class ULSAOutputToolSet<class ULSAKmeans>;
//template class ULSAOutputToolSet<class ULSA2g_MC>;
//template class ULSAOutputToolSet<class ULSA2gx_MC>;
//template class ULSAOutputToolSet<class ULSA2i_MC>;
//template class ULSAOutputToolSet<class ULSA2i2_MC>;
//template class ULSAOutputToolSet<class ULSA2i3_MC>;
//template class ULSAOutputToolSet<class ULSAkmeans2i_MC>;


//template class ULSAOutputToolSet<class ULSA3g_DC>;
//template class ULSAOutputToolSet<class ULSA3gx_DC>;
//template class ULSAOutputToolSet<class ULSA3i_DC>;
//template class ULSAOutputToolSet<class ULSA3i2_DC>;
//template class ULSAOutputToolSet<class ULSA4a_DC>;
//template class ULSAOutputToolSet<class ULSA4b_DC>;
//template class ULSAOutputToolSet<class ULSA4b2_DC>;
//template class ULSAOutputToolSet<class ULSA4b3_DC>;
//template class ULSAOutputToolSet<class ULSAkmeans4b2_DC>;
