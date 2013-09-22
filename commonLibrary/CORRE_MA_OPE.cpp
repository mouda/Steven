#include<iostream>
#include<cv.h>
#include<cmath>
#include<cfloat>
#include<cassert>
#include "CORRE_MA_OPE.h"
using namespace std;
CORRE_MA_OPE::CORRE_MA_OPE(int inTotalNodes, double inCorrelationFactor, float **inDijSQ)
{
  DijSQ = inDijSQ;
  totalNodes = inTotalNodes;
  correlationFac = inCorrelationFactor;
}

double CORRE_MA_OPE::computeLog2Det( double inVariance, bool * inClusterStru)
{
  variance = inVariance;
  //Read information from cluster structure array
  int covMaSize = 0;
  for(int i=0;i<totalNodes;i++)
  {if(inClusterStru[i] == true)covMaSize++;}
  int supSet[covMaSize];
  int cursor = 0;
  for(int i=0;i<totalNodes;i++)
  {
    if (inClusterStru[i]== true )
    {
      supSet[cursor]=i;
      cursor++;
    }
  }
  //cout<<covMaSize<<endl;
  //-----------------------------------------------
  int matrixLength = covMaSize * covMaSize;
  double covAry[matrixLength];
  computeCovMa(covAry,covMaSize ,supSet);
  CvMat covarianceMatrix;
  cvInitMatHeader(&covarianceMatrix, covMaSize,covMaSize,CV_64F,covAry);


  double det = cvDet(&covarianceMatrix);
  //printf("Correlation = %.15f\n",correlationFac);
  //printf("Determinant = %.15f\n",det);
  //printf("log2det = %.15f\n",log2(det));
  assert ((log2(det) <= DBL_MAX && log2(det)>= -DBL_MAX));

  return log2(det);
}

double CORRE_MA_OPE::returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int totalNodes){
    double step =10;
    double start=0;
    bool inClu[totalNodes];
    for(int i=0; i<totalNodes; i++)inClu[i]=true;
    while(1){
        correlationFac=start;

        double tmpCompR=1-(totalNodes*indEntropy+computeLog2Det(1.0, inClu))/(totalNodes*indEntropy);
        //cout<<"Correlation Factor="<<correlationFac<<"; Compression Ratio="<<tmpCompR<<endl;
        //cout<<"total Entropy = "<<(totalNodes*indEntropy)<<";redundancy="<<computeLog2Det(1.0, inClu)<<endl;
        if(tmpCompR>compressionRatio)
            return tmpCompR;
        else
            start+=step;

        if(!(tmpCompR>-10000)){
            cout<<"tmpCOmpr="<<tmpCompR<<endl;
            assert(0);
        }
    }

}
//Purpose: Compute covaraince Matrix
//Input: covarianve matrix diminsion(covMaSize), support set(supSet)
//Output: covAry is covariance array which will be used in later CvMat
void CORRE_MA_OPE::computeCovMa(double* covAry,int covMaSize, int* supSet)
{
  for(int i=0;i<covMaSize;i++)
  {
    for(int j=0;j<covMaSize;j++)
    {
      if(i==j)covAry[i*covMaSize+j] = variance;
      else if(i>j)covAry[i*covMaSize+j] = covAry[j*covMaSize+i];
      else
      {
        covAry[i*covMaSize+j] = variance * exp(-1*((double) DijSQ[supSet[i]][supSet[j]])/correlationFac);
      }
    }
  }
}
