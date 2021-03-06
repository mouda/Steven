#ifndef _CLUSTERSTRUCTURE_
#define _CLUSTERSTRUCTURE_ 
#include <iostream>
#include <vector>
#include <list>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::list;

class ClusterStructure 
{
  public:
    ClusterStructure( const int numNodes, const int maxNumHeads);
    ~ClusterStructure();

    void                          SetRecord( const vector<int>& invecheadname, const list<list<int> >& inlistclumember, const vector<int>& vecAllSupStru);
    void                          SetNotSupport( const int idx){ m_allSupStru.at(idx) = 0; }
    void                          SetSupport( const int idx){ m_allSupStru.at(idx) = 1; }
    void                          SetVecHeadNameByIdx( const int idx, const int name){ m_vecHeadName.at(idx) = name; }
    void                          SetTier1TotalPower( const double tier1Power){m_tier1TotalPower = tier1Power; }

    int                           GetChIdxByName( const int& nodeName ) const;
    int                           GetChNameByName( const int& nodeName ) const;
    int                           GetNumNodes() const {return m_numNodes;}
    int                           GetNumHeads() const {return m_maxNumHeads;}
    double                        GetTier1TotalPower() const { return m_tier1TotalPower;}

    const vector<int>&            GetVecHeadName() const { return m_vecHeadName; }
    const list<list<int> >&       GetListCluMemeber() const{return m_listCluMember; }
    const vector<int>&            GetVecClusterSize() const { return m_vecClusterSize;}
    const vector<int>&            GetVecSupport() const { return m_vecSupport; }
    const vector<int>&            GetAllSupStru() const { return m_allSupStru; }
    const vector<double>&         GetVecPower() const { return m_vecPower; }
    const vector<vector<int> >&   GetMatClusterStru() const { return m_matClusterStru; }
    const list<int>&              GetListUnSupport() const { return m_listUnSupport; }

    void Print() const;

    void                          addNewHeadCs(int inputHead);
    void                          addMemberCs(int headIndex, int memberName,bool iniDone );//headIndex := headName
    void                          discardMemberCs(int headIndex2, int memberName2 );//headIndex := headName
    void                          join_FromHead(int JoiningHeadIndex,int targetH);
    void                          isolate_FromHead(int isoNodeName,int isoCluI,int targetH);

    void                          reverseAdd(int inHeadIndex, int inMember);
    void                          reverseDiscard(int inHeadIndex1, int inMember1);
    void                          reverseJoin(int lastJoinHeadIndex,int lastJoinHead, int lastTargetIndex, std::vector<int>&lastjoinall);
    void                          reverseisolate(int isoNodeName,int isoCluI,int targetH);
    void                          resetSystem();

    bool                          returnIfClusterSmall(int thershold, int &numOfClu);
    int*                          returnHeadPtr(int inputHeadIndex);
    int                           calSupNodes();
  private:
    const int                     m_numNodes;
    const int                     m_maxNumHeads;
    int                           m_numHeads;
    bool                          m_systemon;
    double                        m_tier1TotalPower;

    list<list<int> >              m_listCluMember;
    vector<int>                   m_vecSupport;
    vector<int>                   m_vecHeadName;        // size = # of heads
    vector<int>                   m_vecCHIdxForNodes;   // size = # of nodes
    vector<int>                   m_vecCHNameForNodes;  // size = # of nodes
    vector<double>                m_vecPower;

    vector<int>                   m_vecClusterSize;// include Head itself
    list<list<int> >::iterator    m_iteLastDiscard;
    list<int>                     m_listUnSupport;
    vector<vector<int> >          m_matClusterStru;
    vector<int>                   m_allSupStru;
};
#endif
