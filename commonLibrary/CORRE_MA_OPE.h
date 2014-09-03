#ifndef CORRE_MAOPERATION_H
#define CORRE_MAOPERATION_H
#include "../lib/cholesky.hpp"
class CORRE_MA_OPE
{
public:
  CORRE_MA_OPE(int intotalNodes, double inCorrelationFactor, float ** inDijSQ);
  int totalNodes;
  float **DijSQ;
  double correlationFac;
  double variance;

  double computeLog2Det(double inVariance, bool* inClusterStru );
  double returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int totalNodes);

private:
  void computeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet);//inCovAry is output of function
  double choleskyLogDet( double const * const aryCovariance, const int& dimSize);
  double armaLogDet( double const * const aryCovariance, const int& dimSize);
  double eigenCholeskyLogDet( double const * const aryCovariance, const int& dimSize);

};

#endif // CORRE_MAOPERATION_H
