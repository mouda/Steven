#ifndef CORRE_MAOPERATION_H
#define CORRE_MAOPERATION_H
#include "../lib/cholesky.hpp"
class CORRE_MA_OPE
{
public:
  CORRE_MA_OPE(int intotalNodes, double inCorrelationFactor, double ** inDijSQ);
  int totalNodes;
  double **DijSQ;
  double correlationFac;
  double m_variance;

  double GetVariance() const { return  m_variance;}
  double GetCorrationFactor() const { return correlationFac; }
  double computeLog2Det(double inVariance, bool* inClusterStru ) const;
  double returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int totalNodes);

private:
  void computeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet) const ;//inCovAry is output of function
  void constComputeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet, const double variance) const;
  double choleskyLogDet( double const * const aryCovariance, const int& dimSize);
  double armaLogDet( double const * const aryCovariance, const int& dimSize);
  double eigenCholeskyLogDet( double const * const aryCovariance, const int& dimSize) const;

};

#endif // CORRE_MAOPERATION_H
