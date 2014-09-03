#ifndef IMAGESOURCE_H
#define IMAGESOURCE_H
#include "../lib/cholesky.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

#include <vector>
#include <string>

class ImageSource
{
public:
  ImageSource(int inm_numNodes, double spatialCorrFactor, double temporalCorrFactor, double ** inDijSQ, double qBits);
  ImageSource(int inm_numNodes, const std::string& idtFName, const std::string& corrFName, double ** inDijSQ);
  ~ImageSource();
  int m_numNodes;
  double m_qBits;
  double **DijSQ;
  double m_spatialCorrFac;
  double m_temporalCorrFac;
  double m_startTime;
  double m_currTime;
  double m_variance;
  double m_totalEntropyPerSlot;


  double GetCorrationFactor() const { return m_spatialCorrFac; }
  double GetSpatialCorrelationFactor() const { return m_spatialCorrFac; }

  double GetJointEntropy(const std::vector<int>& vecClusterStru, const std::vector<double>& vecVariance, const double currTime, const double qBits) const;
  double computeLog2Det(double inVariance, bool* inClusterStru ) const;
  double computeLog2Det(double inVariance, bool* inClusterStru ) ;
  double computeLog2Det(double inVariance, const std::vector<int>& vecClusterStru) const;
  double returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int numNodes);

private:
  double matEigenCholeskyLogDet(const std::vector<std::vector<double> >&  covMat , const int& dimSize) const;
  void matConstComputeCovMa(std::vector<std::vector<double> >& covMat, int covMaSize ,int* supSet, const double inVariance) const;
  void ComputeCovMaDiffVariance(std::vector<std::vector<double> >& covMat, int covMaSize ,int* supSet, const std::vector<double>& vecVariance) const;
  void GetCovMaVariance(Eigen::MatrixXd& covMat) const; 
  void computeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet) const ;//inCovAry is output of function
  void constComputeCovMa(double* inCovAry, int inCovMaSize, int* inSupSet, const double variance) const;
  double eigenCholeskyLogDet( double const * const aryCovariance, const int& dimSize) const;

  Eigen::MatrixXd       m_matImageCovariance;
  Eigen::MatrixXd       m_vecImageIndependent;

};

#endif // CORRE_MAOPERATION_H
