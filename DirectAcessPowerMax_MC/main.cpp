#include <iostream>
#include <cstdio>
#include <cstdlib>

#include<iomanip>
#include <algorithm>    // std::sort
#include "../commonLibrary/ULAGENT.h"
#include "../commonLibrary/CORRE_MA_OPE.h"
#include <cmath>
#include<algorithm>
#include<utility>
#include<cassert>
using namespace std;
const static float pathLoss0 = 131.1;// no unit
const static float pathLossAlpha =42.68; //no unit
const static float noise = -100;// no unit
bool comparator ( const pair<int,double>&l, const pair<int,double> &r);
int main(int argc, char* argv[]){
    if(argc!=7){
        cout<<"Please Enter the Map Path"<<endl;
        return 0;
    }
    FILE* fid=fopen(argv[4],"r");
    if(fid==NULL){
        cout<<"Read Map Error"<<endl;
    }
    double powerDbm=atof(argv[1]);
    double powerMax_Watt=pow(10,powerDbm/10)/1000;

    double bandWidthKhz=atof(argv[2]);

    double dataBits=atof(argv[3]);
    double fidratio=atof(argv[5]);
    //double correlationFactor=atof(argv[6]);
    double compressionRatio=atof(argv[6]);
    double n0 = 1*1e-19;


    double topologyRadius=0;
    int totalNodes=0;
    fscanf(fid,"%d %lf",&totalNodes,&topologyRadius);
    double Gib[totalNodes];
    double rateIb_bps[totalNodes];
    double Energy_Joule[totalNodes];
    float** distanceOf2Nodes;

    ULAGENT* nodes = new ULAGENT [totalNodes];
    for (int i=0;i<totalNodes;i++){
            double tempX=0,tempY=0;
            fscanf(fid,"%lf %lf", &tempX, &tempY);
            nodes[i].aryConstructor(i,tempX,tempY);
    }
    distanceOf2Nodes = new float* [totalNodes];

    for(int i=0; i<totalNodes; i++)
    {
        distanceOf2Nodes[i] = new float [totalNodes];

        for(int j = 0; j<totalNodes; j++)
        {
            if (i==j)
            {
                distanceOf2Nodes[i][j]=0.0;
            }
            else if (i>j)
            {
                distanceOf2Nodes[i][j]=distanceOf2Nodes[j][i];
            }
            else
            {
                double tempDib = pow(nodes[i].locX-nodes[j].locX,2) +pow(nodes[i].locY - nodes[j].locY,2);
                distanceOf2Nodes[i][j] = tempDib;
            }
        }
    }




    CORRE_MA_OPE *matrixComputer = new CORRE_MA_OPE(totalNodes, 0.0, distanceOf2Nodes);
    //Set correlation factor by compression ratio
    double tmpCompR = matrixComputer->returnNSetCorrelationFactorByCompressionRatio \
    (compressionRatio,dataBits,static_cast<double>(totalNodes));

    assert(0<fidratio&&fidratio<1.1);
    bool *allsup;
    allsup=new bool [totalNodes];
    for (int i=0;i<totalNodes;i++)allsup[i]=1;
    double wholeSysInfo = totalNodes*dataBits+matrixComputer->computeLog2Det(1.0,allsup);

    double comprRatio=1-wholeSysInfo/(totalNodes*dataBits);
    for (int i=0;i<totalNodes;i++)allsup[i]=0;

    //Calculate Basic Parameters
    pair<int,double>mypair[totalNodes];
    double allInfoMS=0;
    for(int i=0;i<totalNodes;i++){
        Gib[i]=pow(10,-1*(pathLoss0 + pathLossAlpha * (0.5*log10(nodes[i].dibSQ) - 3))/10);
        rateIb_bps[i]=bandWidthKhz*1e3*log2(1+powerMax_Watt*Gib[i]/(n0*bandWidthKhz*1e3));
        mypair[i]=make_pair(i,rateIb_bps[i]);
        Energy_Joule[i] =(dataBits/rateIb_bps[i]) *powerMax_Watt;
        allInfoMS+=(dataBits/rateIb_bps[i]);
    }
    sort(mypair,mypair+totalNodes,comparator);


    double sofarInfo=0;
    int sofarNode=0;

    double tmptotalResource = 0;
    double compensateMS=0;
    for (int i=0;i<totalNodes;i++){
        int index1= mypair[i].first;
        allsup[i]=1;
        sofarNode++;
        sofarInfo=sofarNode*dataBits+matrixComputer->computeLog2Det(1.0,allsup);
        tmptotalResource+=(dataBits/rateIb_bps[index1]);
        if(sofarInfo>wholeSysInfo)break;
    }
    compensateMS=allInfoMS-tmptotalResource;
    compensateMS*=1000;
    for (int i=0;i<totalNodes;i++)allsup[i]=0;
    char str[500];
    sprintf(str,"N%d_R%.1f_CR%.2f_FR%.2f_DirectAcess_PowerMaxMC.txt",totalNodes,topologyRadius,compressionRatio,fidratio);
    FILE *fid2;
    fid2=fopen(str,"w");

    //Caculate transmission resource usage for target scenario
    sofarNode=0;
    sofarInfo=0;
    double totalResource = 0;
    double totalEnergy = 0;
    for (int i=0;i<totalNodes;i++){
        int index= mypair[i].first;
        fprintf(fid2,"%d\n",index+1);

        allsup[i]=1;
        sofarNode++;
        sofarInfo=sofarNode*dataBits+matrixComputer->computeLog2Det(1.0,allsup);
        totalResource+=(dataBits/rateIb_bps[index]);
        Energy_Joule[index] =(dataBits/rateIb_bps[index]) *powerMax_Watt;
        totalEnergy+=Energy_Joule[index];
        if((sofarInfo>fidratio*wholeSysInfo)&&(static_cast<double>(sofarNode)>static_cast<double>(totalNodes)*fidratio))break;
    }
    fclose(fid2);
    totalResource *=1000;
    cout<<totalNodes<<" "<<dataBits<<" "<<tmpCompR<<" "<<wholeSysInfo<<" "<<sofarInfo<<" "<<sofarNode<<" "<<fidratio \
    <<" "<<comprRatio<<" "<<totalResource<<" "<<totalEnergy<<endl;

    /*for (int i=0;i<totalNodes;i++){
      cout<<setw(3)<<i<<setw(3)<<"  "<<dataBits<<" Gib "<<setw(12)<<Gib[i]<<" "<<setw(12)<<rateIb_bps[i]<<"(bps)  "<<setw(12)<<\
      dataBits/rateIb_bps[i]*1000<<"(ms) "<< \
      setw(4)<<powerMax_Watt<<"(Watt) "<<setw(12)<<Energy_Joule[i]<<"(Joule)"<<endl;
    }*/
    for(int i=0; i<totalNodes; i++) delete [] distanceOf2Nodes[i];
    delete [] distanceOf2Nodes;
    delete matrixComputer;
    delete [] nodes;
    delete [] allsup;
}
    bool comparator ( const pair<int,double>&l, const pair<int,double> &r)
   { return l.second > r.second; }
