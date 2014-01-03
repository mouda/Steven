#ifndef CORRE_MAOPERATION_H
#define CORRE_MAOPERATION_H
#include "../lib/cholesky.hpp"
#include <vector>
using std::vector;
class CORRE_MA_OPE
{
public:
  CORRE_MA_OPE(int intotalNodes, double inCorrelationFactor, double ** inDijSQ);
  ~CORRE_MA_OPE();
  int totalNodes;
  double **DijSQ;
  double correlationFac;
  double m_variance;

  double GetVariance() const { return  m_variance;}
  double GetCorrationFactor() const { return correlationFac; }
  double computeLog2Det(double inVariance, bool* inClusterStru ) const;
  double returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int totalNodes);

private:
  double matEigenCholeskyLogDet(const vector<vector<double> >&  covMat , const int& dimSize) const;
  void matConstComputeCovMa(vector<vector<double> >& covMat, int covMaSize ,int* supSet, const double inVariance) const;
  void computeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet) const ;//inCovAry is output of function
  void constComputeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet, const double variance) const;
  double choleskyLogDet( double const * const aryCovariance, const int& dimSize) const;
  double armaLogDet( double const * const aryCovariance, const int& dimSize);
  double eigenCholeskyLogDet( double const * const aryCovariance, const int& dimSize) const;

};

#endif // CORRE_MAOPERATION_H
