#ifndef _GAUSSIANFIELD_
#define _GAUSSIANFIELD_

#include <vector>

class GaussianField
{
  public:
    GaussianField(int totalNodes, double inCorrelationFactor, double ** inDijSQ);
    ~GaussianField();

    double GetEntropy( const std::vector<double>& vecCurEntropy, double time);
    double GetVariance() const { return  m_variance;}
    double GetCorrationFactor() const { return correlationFac; }
    double GetDijSQByPair( const int i, const int j) const { return m_DijSQ[i][j]; }

    double computeLog2Det(double inVariance, const std::vector<int>& vecClusterStru) const;
    double returnNSetCorrelationFactorByCompressionRatio(double compressionRatio,double indEntropy, int totalNodes);

  private:

    double  matEigenCholeskyLogDet(
        const std::vector<std::vector<double> >&  covMat, 
        const int& dimSize
        ) const;
    void    matConstComputeCovMa(
        std::vector<std::vector<double> >& covMat, 
        int covMaSize,
        int* supSet, 
        const double inVariance) const;

    int       m_numNodes;
    double**  m_DijSQ;
    double    correlationFac;
    double    m_variance;

};

#endif
