#include "gaussianField.h"

GaussianField::GaussianField(int totalNodes, double inCorrelationFactor, double ** inDijSQ)
{
  m_DijSQ = inDijSQ;
  totalNodes = inTotalNodes;
  correlationFac = inCorrelationFactor;
}

GaussianField::~GaussianField()
{

}

double
GaussianField::GetEntropy( const vector<double>& vecCurEntropy, double time )
{

} 

double
GaussianField::computeLog2Det(double inVariance, const vector<int>& vecClusterStru) const
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

double 
GaussianField::returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int totalNodes)
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
        if( tmpCompR > compressionRatio ){
          return tmpCompR;
        }
        else
            start+=step;

        if(!(tmpCompR>-10000)){
            assert(0);
        }
    }

    delete [] inClu;
    return correlationFac; 
}

double 
GaussianField::matEigenCholeskyLogDet( const vector<vector<double> >& covMat, const int& dimSize) const
{
  Eigen::MatrixXd covMatrix(dimSize,dimSize);
  for (int i = 0; i < dimSize; ++i) {
    for (int j = 0; j < dimSize; ++j) {
      covMatrix(i,j) = covMat[ i][ j ];  
    }
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

void 
GaussianField::matConstComputeCovMa(vector<vector<double> >& covMat, int covMaSize ,int* supSet, const double inVariance) const
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
        covMat[i][j] = inVariance * exp(-1*((double) m_DijSQ[supSet[i]][supSet[j]])/correlationFac);
      }
    }
  }

}
