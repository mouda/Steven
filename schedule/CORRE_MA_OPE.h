#ifndef CORRE_MAOPERATION_H
#define CORRE_MAOPERATION_H
#include "../lib/cholesky.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

#include <vector>
using std::vector;
class CORRE_MA_OPE
{
public:
  CORRE_MA_OPE(int inm_numNodes, double spatialCorrFactor, double temporalCorrFactor, double ** inDijSQ, double qBits);
  ~CORRE_MA_OPE();
  int m_numNodes;
  double m_qBits;
  double **DijSQ;
  double m_spatialCorrFac;
  double m_temporalCorrFac;
  double m_startTime;
  double m_currTime;
  double m_variance;
  double m_totalEntropyPerSlot;

  double SetCurrTime( const double myTime ) { m_currTime = myTime;}
  double GetCurrTime() const { return m_currTime; }

  double GetVariance() const { return  m_variance;}
  double GetCorrationFactor() const { return m_spatialCorrFac; }
  double GetSpatialCorrelationFactor() const { return m_spatialCorrFac; }
  double GetTemporalCorrelationFactor() const { return m_temporalCorrFac; }
  double GetDijSQByPair( const int i, const int j) const { return DijSQ[i][j]; }

  void   UpdateVariance(const vector<double>& curVecVariance, vector<double>& nextVecVariance, const vector<int>& vecSupport, const double timeDiff) const;
  double GetJointEntropy(const vector<int>& vecClusterStru, const vector<double>& vecVariance, const double currTime, const double qBits) const;
  double GetRateDistortion(const vector<int>& vecClusterStru, const vector<double>& vecVariance, const double currTime, const double qBits) const;
  double computeLog2Det(double inVariance, bool* inClusterStru ) const;
  double computeLog2Det(double inVariance, const vector<int>& vecClusterStru) const;
  double returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int numNodes);

private:
  double matEigenCholeskyLogDet(const vector<vector<double> >&  covMat , const int& dimSize) const;
  void matConstComputeCovMa(vector<vector<double> >& covMat, int covMaSize ,int* supSet, const double inVariance) const;
  void ComputeCovMaDiffVariance(vector<vector<double> >& covMat, int covMaSize ,int* supSet, const vector<double>& vecVariance) const;
  void GetCovMaVariance(Eigen::MatrixXd& covMat) const; 
  void computeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet) const ;//inCovAry is output of function
  void constComputeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet, const double variance) const;
  double armaLogDet( double const * const aryCovariance, const int& dimSize);
  double eigenCholeskyLogDet( double const * const aryCovariance, const int& dimSize) const;

};

#endif // CORRE_MAOPERATION_H
