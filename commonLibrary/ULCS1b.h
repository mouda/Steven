/*
  File: ULCS1b1b.h
  Brief: Record all the coalition structure information of SA
        1b: keep one more new structure
        allSupStru[totalNodes]

  Author: Steven
  Date: 2012/09/27.
  Latest Progress:

  Abbreviation:

  Quote:
  Index means the array index, so the data store in the vecHeadName is head name. However we need index to access it.
  Name means the number of each machine

*/
#ifndef ULCS1b_H
#define ULCS1b_H
#include<vector>
#include<list>

using namespace std;
class ULCS1b
{
public:
  ULCS1b();
  ULCS1b( int inputTotalNodes, int inputMaxChNum);
  ~ULCS1b();
  void addNewHeadCs(int inputHead);
  void addMemberCs(int headIndex, int memberName,bool iniDone );//headIndex := headName
  void discardMemberCs(int headIndex2, int memberName2 );//headIndex := headName
  void join_FromHead(int JoiningHeadIndex,int targetH);
  void isolate_FromHead(int isoNodeName,int isoCluI,int targetH);


  void reverseAdd(int inHeadIndex, int inMember);
  void reverseDiscard(int inHeadIndex1, int inMember1);
  void reverseJoin(int lastJoinHeadIndex,int lastJoinHead, int lastTargetIndex, std::vector<int>&lastjoinall);
  void reverseisolate(int isoNodeName,int isoCluI,int targetH);
  void resetSystem();

  bool returnIfClusterSmall(int thershold, int &numOfClu);



  int calSupNodes();
  void writeCsFile(char* fileName);
  int* returnHeadPtr(int inputHeadIndex);
  int maxChNum;
  int totalNodes;
  //private:
  bool systemon;
  vector <int> vecHeadName;// vexHeadIndex store the name of the machine which is CH
  vector <int> vecClusterSize;// include Head itself

  list<list<int> >::iterator iteLastDiscard;
  list<list<int> > *listCluMember;
  list <int> *listUnSupport;
  bool ** clusterStru;
  bool * allSupStru;
//        ULCS1b& operator= (const ULCS1b &rhs);
};
#endif
