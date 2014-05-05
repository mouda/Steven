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

    void  SetRecord( const vector<int>& invecheadname, 
          const list<list<int> >& inlistclumember);

    int                     GetChIdxByName( const int& nodeName ) const;
    int                     GetChNameByName( const int& nodeName ) const;
    int                     GetNumNodes() const {return m_numNodes;}
    int                     GetNumHeads() const {return m_maxNumHeads;}

    const vector<int>&      GetVecHeadName() const { return m_vecHeadName; }
    const list<list<int> >& GetListCluMemeber() const{return m_listCluMember; }
    const vector<int>&      GetVecSupport() const { return m_vecSupport; }
    const vector<double>&   GetVecPower() const { return m_vecPower; }

    void Print() const;

  private:
    const int   m_numNodes;
    const int   m_maxNumHeads;
    int         m_numHeads;

    list<list<int> >  m_listCluMember;
    vector<int>       m_vecSupport;
    vector<int>       m_vecHeadName;        // size = # of heads
    vector<int>       m_vecCHIdxForNodes;   // size = # of nodes
    vector<int>       m_vecCHNameForNodes;  // size = # of nodes
    vector<double>    m_vecPower;
};
#endif
