#include <iostream>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <armadillo>
#include <eigen3/Eigen/Cholesky>

#include "CORRE_MA_OPE.h"
#include "../lib/cholesky.hpp"
using namespace std;
//using namespace boost::numeric;

CORRE_MA_OPE::CORRE_MA_OPE(int inTotalNodes, double inCorrelationFactor, double **inDijSQ)
{
  DijSQ = inDijSQ;
  totalNodes = inTotalNodes;
  correlationFac = inCorrelationFactor;
}

CORRE_MA_OPE::~CORRE_MA_OPE()
{

}

double CORRE_MA_OPE::computeLog2Det( double inVariance, bool * inClusterStru) const
{
  //m_variance = inVariance;
  //Read information from cluster structure array
  int covMaSize = 0;
  for(int i=0;i<totalNodes;i++)
  {if(inClusterStru[i] == true)covMaSize++;}
  int* supSet= new int [covMaSize];
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
  double *covAry = new double [matrixLength];
  vector<vector<double> > covMat(covMaSize,vector<double>(covMaSize));
  //computeCovMa(covAry,covMaSize ,supSet);
  matConstComputeCovMa(covMat, covMaSize ,supSet, inVariance);

//  cout << "new   : " << choleskyLogDet(covAry,covMaSize) << endl;
//  cout << "arma  : " << armaLogDet(covAry, covMaSize) << endl;
//  cout << "Eigen : " << eigenCholeskyLogDet(covAry, covMaSize) << endl;

//  return armaLogDet(covAry, covMaSize);
    delete [] covAry;
    delete [] supSet;
    
    return matEigenCholeskyLogDet(covMat, covMaSize);
    //return choleskyLogDet(covAry,covMaSize);
}

double CORRE_MA_OPE::computeLog2Det( double inVariance, const vector<int>& vecClusterStru) const
{
  int covMaSize = 0;
  for(int i = 0; i < totalNodes; ++i) {
    if(vecClusterStru[i] == true) ++covMaSize;
  }
  int* supSet= new int [covMaSize];
  int cursor = 0;
  for(int i=0;i<totalNodes;i++)
  {
    if (vecClusterStru[i]== true )
    {
      supSet[cursor]=i;
      cursor++;
    }
  }
  int matrixLength = covMaSize * covMaSize;
  vector<vector<double> > covMat(covMaSize,vector<double>(covMaSize));
  matConstComputeCovMa(covMat, covMaSize ,supSet, inVariance);
  delete [] supSet;

  return matEigenCholeskyLogDet(covMat, covMaSize);
}

double CORRE_MA_OPE::returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int totalNodes)
{
    double step =100;
    double start=0;
    bool* inClu = new bool [totalNodes];
    for(int i = 0; i < totalNodes; ++i )
      inClu[i]=true;
    while(1){
        correlationFac=start;
        double redundancy = computeLog2Det(1.0, inClu);
        redundancy = computeLog2Det(1.0, inClu);
        double tmpCompR=1-(totalNodes*indEntropy+redundancy)/(totalNodes*indEntropy);
//        cout << "Compression Ration: " << compressionRatio << endl;
//        cout << "Correlation Factor: " << correlationFac <<  ";Compression Ratio=" << tmpCompR <<endl;
//        cout << "total Entropy: " << (totalNodes * indEntropy) << ";redundancy=" << redundancy <<endl;
//        cout << "Ratio: " << endl;
        if( tmpCompR > compressionRatio ){
//          cerr << tmpCompR << endl;
          return tmpCompR;
        }
        else
            start+=step;

        if(!(tmpCompR>-10000)){
//            cerr<<"tmpCOmpr="<<tmpCompR<<endl;
            assert(0);
        }
    }
//    cout << correlationFac << endl;

    delete [] inClu;

}
//Purpose: Compute covaraince Matrix
//Input: covarianve matrix diminsion(covMaSize), support set(supSet)
//Output: covAry is covariance array which will be used in later CvMat
void CORRE_MA_OPE::computeCovMa(double* covAry,int covMaSize, int* supSet) const 
{
  for(int i=0;i<covMaSize;i++)
  {
    for(int j=0;j<covMaSize;j++)
    {
      if(i==j)covAry[i*covMaSize+j] = m_variance;
      else if(i>j)covAry[i*covMaSize+j] = covAry[j*covMaSize+i];
      else
      {
        covAry[i*covMaSize+j] = m_variance * exp(-1*((double) DijSQ[supSet[i]][supSet[j]])/correlationFac);
      }
    }
  }
}

void CORRE_MA_OPE::constComputeCovMa(double* covAry,int covMaSize, int* supSet, const double variance ) const 
{
  for(int i=0;i<covMaSize;i++)
  {
    for(int j=0;j<covMaSize;j++)
    {
      if( i == j ) {
        covAry[i*covMaSize+j] = variance;
      }
      else if( i > j ) {
        covAry[i*covMaSize+j] = covAry[j*covMaSize+i];
      }
      else {
        covAry[i*covMaSize+j] = variance * exp(-1*((double) DijSQ[supSet[i]][supSet[j]])/correlationFac);
      }
    }
  }
}

void CORRE_MA_OPE::matConstComputeCovMa(vector<vector<double> >& covMat, int covMaSize ,int* supSet, const double inVariance) const
{
  for(int i=0;i<covMaSize;i++)
  {
    for(int j=0;j<covMaSize;j++)
    {
      if( i == j ) {
        covMat[i][j] = inVariance;
      }
      else if( i > j ) {
        covMat[i][j] = covMat[j][i];
      }
      else {
        covMat[i][j] = inVariance * exp(-1*((double) DijSQ[supSet[i]][supSet[j]])/correlationFac);
      }
    }
  }

}
double CORRE_MA_OPE::choleskyLogDet( double const * const aryCovariance, const int& dimSize) const  
{
  boost::numeric::ublas::matrix<double> covMatrix(dimSize,dimSize);
  for (int i = 0; i < dimSize; ++i) {
    for (int j = 0; j < dimSize; ++j) {
      covMatrix(i,j) = aryCovariance[ i * dimSize + j ];  
    }
  }
  boost::numeric::ublas::matrix<double> TRM (dimSize, dimSize);
  cholesky_decompose(covMatrix, TRM);
  double logDet = 0.0;
  for (int i = 0; i < dimSize; ++i) {
    for (int j = 0; j < dimSize; ++j) {
      if (i == j) {
        logDet += log2(TRM(i,j));
      }
    }
  }
  return 2*logDet;

}

double CORRE_MA_OPE::armaLogDet( double const * const aryCovariance, const int& dimSize)
{
  arma::Mat<double> covMatrix(dimSize, dimSize);
  for (int i = 0; i < dimSize; ++i) {
    for (int j = 0; j < dimSize; ++j) {
      covMatrix(i,j) = aryCovariance[ i * dimSize + j ];  
    }
  }
  return log2(arma::det(covMatrix));
  

}

double CORRE_MA_OPE::eigenCholeskyLogDet( double const * const aryCovariance, const int& dimSize) const
{
  Eigen::MatrixXd covMatrix(dimSize,dimSize);
  for (int i = 0; i < dimSize; ++i) {
    for (int j = 0; j < dimSize; ++j) {
      covMatrix(i,j) = aryCovariance[ i * dimSize + j ];  
  //    cout << aryCovariance[i*dimSize +j] << ' ';
    }
  //  cout << endl;
  }
  Eigen::MatrixXd TRM( covMatrix.llt().matrixL() );
  double logDet = 0.0;
  for (int i = 0; i < dimSize; ++i) {
    for (int j = 0; j < dimSize; ++j) {
      if (i == j) {
        logDet += log2(TRM(i,j));
      }
    }
  }

  return 2*logDet;
}

double CORRE_MA_OPE::matEigenCholeskyLogDet( const vector<vector<double> >& covMat, const int& dimSize) const
{
  Eigen::MatrixXd covMatrix(dimSize,dimSize);
  for (int i = 0; i < dimSize; ++i) {
    for (int j = 0; j < dimSize; ++j) {
      covMatrix(i,j) = covMat[ i][ j ];  
  //    cout << aryCovariance[i*dimSize +j] << ' ';
    }
  //  cout << endl;
  }
  Eigen::MatrixXd TRM( covMatrix.llt().matrixL() );
  double logDet = 0.0;
  for (int i = 0; i < dimSize; ++i) {
    for (int j = 0; j < dimSize; ++j) {
      if (i == j) {
        logDet += log2(TRM(i,j));
      }
    }
  }

  return 2*logDet;
}
