#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "imageSource.h"
using namespace std;
//using namespace boost::numeric;

ImageSource::ImageSource(int inTotalNodes, double spatialCorrFactor, 
    double temporalCorrFactor, 
    double **inDijSQ, double qBits):
  m_qBits(qBits),
  m_temporalCorrFac(temporalCorrFactor),
  m_numNodes(inTotalNodes),
  DijSQ(inDijSQ),
  m_totalEntropyPerSlot(0)
{
  m_matImageCovariance = Eigen::MatrixXd::Zero(inTotalNodes, inTotalNodes); 
  m_vecImageIndependent = Eigen::MatrixXd::Zero(inTotalNodes, 1);
  fstream myImageCovFile("paper720_30cam_corrMatrix.txt", std::ios::in) ; 
  fstream myImageIdtFile("paper720_30cam_indepByte.txt",std::ios::in);
  double tmp = 0;
  if (myImageCovFile.good()) {
    for (int i = 0; i < inTotalNodes; ++i) {
      for (int j = 0; j < inTotalNodes; ++j) {
        myImageCovFile >> tmp;
        m_matImageCovariance(i,j) = tmp;
      }
    }
    for (int i = 0; i < inTotalNodes; ++i) {
      myImageIdtFile >> tmp;
      m_vecImageIndependent(i) = tmp;
    }
  }
  else {
    cerr << "Error: Cannot read image source" << endl;
    myImageCovFile.close();
    myImageIdtFile.close();
    assert(0);
  }
  myImageIdtFile.close();
  myImageCovFile.close();
}

ImageSource::ImageSource(int inTotalNodes, const std::string& idtFName, const std::string& corrFName, double ** inDijSQ):
  m_numNodes(inTotalNodes),
  m_qBits(0.0),
  m_temporalCorrFac(0.0),
  DijSQ(inDijSQ),
  m_totalEntropyPerSlot(0)
{
  m_matImageCovariance = Eigen::MatrixXd::Zero(inTotalNodes, inTotalNodes); 
  m_vecImageIndependent = Eigen::MatrixXd::Zero(inTotalNodes, 1);
  fstream myImageCovFile(corrFName.c_str(), std::ios::in) ; 
  fstream myImageIdtFile(idtFName.c_str(),std::ios::in);
  double tmp = 0;
  if (myImageCovFile.good() && myImageIdtFile.good()) {
    for (int i = 0; i < inTotalNodes; ++i) {
      for (int j = 0; j < inTotalNodes; ++j) {
        myImageCovFile >> tmp;
        m_matImageCovariance(i,j) = tmp;
      }
    }
    for (int i = 0; i < inTotalNodes; ++i) {
      myImageIdtFile >> tmp;
      m_vecImageIndependent(i) = tmp;
    }
  }
  else {
    cerr << "Error: Cannot read image source" << endl;
    myImageCovFile.close();
    myImageIdtFile.close();
    assert(0);
  }
  myImageIdtFile.close();
  myImageCovFile.close();

}

ImageSource::~ImageSource()
{

}

double ImageSource::computeLog2Det( double inVariance, bool * inClusterStru) const
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
  std::vector<std::vector<double> > covMat(covMaSize,std::vector<double>(covMaSize));
  //computeCovMa(covAry,covMaSize ,supSet);
  matConstComputeCovMa(covMat, covMaSize ,supSet, inVariance);

//  cout << "arma  : " << armaLogDet(covAry, covMaSize) << endl;
//  cout << "Eigen : " << eigenCholeskyLogDet(covAry, covMaSize) << endl;

//  return armaLogDet(covAry, covMaSize);
    delete [] covAry;
    delete [] supSet;
    
    return matEigenCholeskyLogDet(covMat, covMaSize);
}

double ImageSource::computeLog2Det( double inVariance, bool * inClusterStru) 
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
  std::vector<std::vector<double> > covMat(covMaSize,std::vector<double>(covMaSize));
  //computeCovMa(covAry,covMaSize ,supSet);
  matConstComputeCovMa(covMat, covMaSize ,supSet, inVariance);

//  cout << "arma  : " << armaLogDet(covAry, covMaSize) << endl;
//  cout << "Eigen : " << eigenCholeskyLogDet(covAry, covMaSize) << endl;

//  return armaLogDet(covAry, covMaSize);
    delete [] covAry;
    delete [] supSet;
    
    return matEigenCholeskyLogDet(covMat, covMaSize);
}

double ImageSource::computeLog2Det( double inVariance, const std::vector<int>& vecClusterStru) const
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
  std::vector<std::vector<double> > covMat(covMaSize,std::vector<double>(covMaSize));
  matConstComputeCovMa(covMat, covMaSize ,supSet, inVariance);
  delete [] supSet;

  return matEigenCholeskyLogDet(covMat, covMaSize);
}

/* @brief    given the subset of machine to report the joint entropy
 * @param   support set vecClusterStru 
 * @retval  entropy 
 */


double 
ImageSource::GetJointEntropy(const std::vector<int>& vecClusterStru, const std::vector<double>& vecVariance, const double currTime, const double qBits) const
{
  int     minIframeIdx = -1;
  int     mySupportCount = 0;
  double  minIframeBytes = DBL_MAX;
  for (int i = 0; i < vecClusterStru.size(); ++i) {
    if (vecClusterStru.at(i) && m_vecImageIndependent(i) < minIframeBytes ) {
      minIframeIdx = i;
      minIframeBytes = m_vecImageIndependent(i);
    }
    if (vecClusterStru.at(i)) ++mySupportCount;
  }
  std::vector<int> myDeterminedNode;
  myDeterminedNode.push_back(minIframeIdx);
  double totalBytes = minIframeBytes;
  while (myDeterminedNode.size() < mySupportCount) {
    int minPframeIdx = -1;
    int minPframeBytes = DBL_MAX;
    for (int i = 0; i < myDeterminedNode.size(); ++i) {
      for (int j = 0; j < m_numNodes; ++j) {
        if (vecClusterStru.at(j) && 
            std::find(myDeterminedNode.begin(), myDeterminedNode.end(), j) == myDeterminedNode.end() &&
            m_matImageCovariance(myDeterminedNode.at(i),j) < minPframeBytes
            ) {
          minPframeIdx = j;
          minPframeBytes = m_matImageCovariance(myDeterminedNode.at(i),j);
        }
        if (vecClusterStru.at(j)) {
        }
      }
    }
    totalBytes += minPframeBytes;
    myDeterminedNode.push_back(minPframeIdx);
  }
  
  
  return totalBytes*8;
}

double ImageSource::returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int numNodes)
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
void ImageSource::computeCovMa(double* covAry,int covMaSize, int* supSet) const 
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

void ImageSource::constComputeCovMa(double* covAry,int covMaSize, int* supSet, const double variance ) const 
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

void ImageSource::matConstComputeCovMa(std::vector<std::vector<double> >& covMat, int covMaSize ,int* supSet, const double inVariance) const
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
ImageSource::ComputeCovMaDiffVariance(std::vector<std::vector<double> >& covMat, int covMaSize ,int* supSet, const std::vector<double>& vecVariance) const
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
ImageSource::GetCovMaVariance(Eigen::MatrixXd& covMat) const
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


double ImageSource::eigenCholeskyLogDet( double const * const aryCovariance, const int& dimSize) const
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

double ImageSource::matEigenCholeskyLogDet( const std::vector<std::vector<double> >& covMat, const int& dimSize) const
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
