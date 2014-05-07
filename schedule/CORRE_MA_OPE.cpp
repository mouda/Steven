#include <iostream>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <armadillo>

#include "CORRE_MA_OPE.h"
#include "../lib/cholesky.hpp"
using namespace std;
//using namespace boost::numeric;

CORRE_MA_OPE::CORRE_MA_OPE(int inTotalNodes, double spatialCorrFactor, 
    double temporalCorrFactor, 
    double **inDijSQ, double qBits):
  m_qBits(qBits),
  m_temporalCorrFac(temporalCorrFactor),
  m_numNodes(inTotalNodes),
  DijSQ(inDijSQ),
  m_totalEntropyPerSlot(0)
{
  double idtEntropy = 0.5*log2(2*3.1415*exp(1)) + qBits;
  returnNSetCorrelationFactorByCompressionRatio(spatialCorrFactor, idtEntropy, static_cast<double>(m_numNodes) );
  Eigen::MatrixXd covMat(m_numNodes, m_numNodes);
  GetCovMaVariance(covMat);
  m_totalEntropyPerSlot = covMat.determinant();
}

CORRE_MA_OPE::CORRE_MA_OPE(int inTotalNodes, double spatialCorrFactor, 
    double **inDijSQ, double qBits):
  m_qBits(qBits),
  m_temporalCorrFac(0),
  m_numNodes(inTotalNodes),
  DijSQ(inDijSQ),
  m_totalEntropyPerSlot(0)
{
  double idtEntropy = 0.5*log2(2*3.1415*exp(1)) + qBits;
  returnNSetCorrelationFactorByCompressionRatio(spatialCorrFactor, idtEntropy, static_cast<double>(m_numNodes) );
  Eigen::MatrixXd covMat(m_numNodes, m_numNodes);
  GetCovMaVariance(covMat);
  m_totalEntropyPerSlot = covMat.determinant();
}

CORRE_MA_OPE::~CORRE_MA_OPE()
{

}

double CORRE_MA_OPE::computeLog2Det( double inVariance, bool * inClusterStru) const
{
  //m_variance = inVariance;
  //Read information from cluster structure array
  int covMaSize = 0;
  for(int i=0;i<m_numNodes;i++)
  {if(inClusterStru[i] == true)covMaSize++;}
  int* supSet= new int [covMaSize];
  int cursor = 0;
  for(int i=0;i<m_numNodes;i++)
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

//  cout << "arma  : " << armaLogDet(covAry, covMaSize) << endl;
//  cout << "Eigen : " << eigenCholeskyLogDet(covAry, covMaSize) << endl;

//  return armaLogDet(covAry, covMaSize);
    delete [] covAry;
    delete [] supSet;
    
    return matEigenCholeskyLogDet(covMat, covMaSize);
}

double CORRE_MA_OPE::computeLog2Det( double inVariance, bool * inClusterStru) 
{
  //m_variance = inVariance;
  //Read information from cluster structure array
  int covMaSize = 0;
  for(int i=0;i<m_numNodes;i++)
  {if(inClusterStru[i] == true)covMaSize++;}
  int* supSet= new int [covMaSize];
  int cursor = 0;
  for(int i=0;i<m_numNodes;i++)
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

//  cout << "arma  : " << armaLogDet(covAry, covMaSize) << endl;
//  cout << "Eigen : " << eigenCholeskyLogDet(covAry, covMaSize) << endl;

//  return armaLogDet(covAry, covMaSize);
    delete [] covAry;
    delete [] supSet;
    
    return matEigenCholeskyLogDet(covMat, covMaSize);
}

double CORRE_MA_OPE::computeLog2Det( double inVariance, const vector<int>& vecClusterStru) const
{
  int covMaSize = 0;
  for(int i = 0; i < m_numNodes; ++i) {
    if(vecClusterStru[i] == 1) ++covMaSize;
  }
  int* supSet= new int [covMaSize];
  int cursor = 0;
  for(int i=0;i<m_numNodes;i++)
  {
    if (vecClusterStru[i]== 1 )
    {
      supSet[cursor]=i;
      cursor++;
    }
  }
  vector<vector<double> > covMat(covMaSize,vector<double>(covMaSize));
  matConstComputeCovMa(covMat, covMaSize ,supSet, inVariance);
  delete [] supSet;

  return matEigenCholeskyLogDet(covMat, covMaSize);
}

void
CORRE_MA_OPE::UpdateVariance(const vector<double>& curVecVariance, vector<double>& nextVecVariance, const vector<int>& vecSupport, vector<int>& vecSlots, const double timeDiff) const
{
  for (int i = 0; i < m_numNodes; ++i) {
    if (vecSupport[i] == 1) {
      vecSlots.at(i) = 1;
      ///nextVecVariance.at(i) = 1.0 * (1 - pow(exp(-1*timeDiff/m_temporalCorrFac),2)) ;
      nextVecVariance.at(i) = curVecVariance.at(i);
    }
    else{
      if (vecSlots.at(i) == 0) {
        nextVecVariance.at(i) = 1.0 ;
      }
      else {
        vecSlots.at(i) += 1;
        //nextVecVariance.at(i) = 1.0 * (1 - pow(exp(-1*timeDiff*static_cast<double>(vecSlots[i])/m_temporalCorrFac),2));
        nextVecVariance.at(i) = curVecVariance.at(i);
      }
    }
  }
}

double 
CORRE_MA_OPE::GetJointEntropy(const vector<int>& vecClusterStru, const vector<double>& vecVariance, const double currTime, const double qBits) const
{
  double idtEntropy = 0.0;
  for (int i = 0; i < m_numNodes; ++i) {
    if (vecClusterStru[i] == 1) {
      idtEntropy += 0.5*log2(2*3.1415*exp(1)) + log2(vecVariance[i]) + qBits;
    }
  }

  int covMaSize = 0;
  for(int i = 0; i < m_numNodes; ++i) {
    if(vecClusterStru[i] == 1) ++covMaSize;
  }
  int* supSet= new int [covMaSize];
  int cursor = 0;
  for(int i=0;i<m_numNodes;i++)
  {
    if (vecClusterStru[i]== 1 )
    {
      supSet[cursor]=i;
      cursor++;
    }
  }
  vector<vector<double> > covMat(covMaSize,vector<double>(covMaSize));
  ComputeCovMaDiffVariance(covMat, covMaSize, supSet, vecVariance);
  double redundancy = matEigenCholeskyLogDet(covMat, covMaSize);

  delete [] supSet;
  return idtEntropy + redundancy;
}

double
CORRE_MA_OPE::GetRateDistortion(const vector<int>& vecClusterStru, const vector<double>& vecVariance, const double currTime, const double qBits) const
{
  double entropy = GetJointEntropy(vecClusterStru, vecVariance, currTime, qBits);
  int* supSet = new int [m_numNodes];
  for (int i = 0; i < m_numNodes; ++i) {
    supSet[i] = i;
  }
  delete [] supSet;
  return pow( m_totalEntropyPerSlot / pow(2,2*entropy), 1.0/static_cast<double>(m_numNodes));
}

double CORRE_MA_OPE::returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int numNodes)
{
    double step =100;
    double start=0;
    bool* inClu = new bool [m_numNodes];
    for(int i = 0; i < m_numNodes; ++i )
      inClu[i]=true;
    while(1){
        m_spatialCorrFac=start;
        double redundancy = computeLog2Det(1.0, inClu);
        redundancy = computeLog2Det(1.0, inClu);
        double tmpCompR=1-(m_numNodes*indEntropy+redundancy)/(m_numNodes*indEntropy);
//        cout << "Compression Ration: " << compressionRatio << endl;
//        cout << "Correlation Factor: " << m_spatialCorrFac <<  ";Compression Ratio=" << tmpCompR <<endl;
//        cout << "total Entropy: " << (m_numNodes * indEntropy) << ";redundancy=" << redundancy <<endl;
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
//    cout << m_spatialCorrFac << endl;

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
        covAry[i*covMaSize+j] = m_variance * exp(-1*((double) DijSQ[supSet[i]][supSet[j]])/m_spatialCorrFac);
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
        covAry[i*covMaSize+j] = variance * exp(-1*((double) DijSQ[supSet[i]][supSet[j]])/m_spatialCorrFac);
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
        covMat[i][j] = inVariance * exp(-1*((double) DijSQ[supSet[i]][supSet[j]])/m_spatialCorrFac);
      }
    }
  }
}

void 
CORRE_MA_OPE::ComputeCovMaDiffVariance(vector<vector<double> >& covMat, int covMaSize ,int* supSet, const vector<double>& vecVariance) const
{
  for(int i=0;i<covMaSize;i++)
  {
    for(int j=0;j<covMaSize;j++)
    {
      if( i == j ) {
        covMat[i][j] = vecVariance.at(supSet[i]) * vecVariance.at(supSet[i]);
      }
      else if( i > j ) {
        covMat[i][j] = covMat[j][i];
      }
      else {
        covMat[i][j] = vecVariance.at(supSet[i])*vecVariance.at(supSet[j])* exp(-1*((double) DijSQ[supSet[i]][supSet[j]])/m_spatialCorrFac);
      }
    }
  }
}

void
CORRE_MA_OPE::GetCovMaVariance(Eigen::MatrixXd& covMat) const
{
  for(int i=0;i<m_numNodes;i++)
  {
    for(int j=0;j<m_numNodes;j++)
    {
      if( i == j ) {
        covMat(i,j) = 1.0;
      }
      else if( i > j ) {
        covMat(i,j) = covMat(i,j);
      }
      else {
        covMat(i,j) = 1.0 * exp(-1*((double) DijSQ[i][j])/m_spatialCorrFac);
      }
    }
  }
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
