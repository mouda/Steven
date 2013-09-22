#ifndef CORRE_MAOPERATION_H
#define CORRE_MAOPERATION_H
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

};

#endif // CORRE_MAOPERATION_H
