#ifndef _CLUSTERSTRUCTURE_
#define _CLUSTERSTRUCTURE_ 
#include <vector>
#include <list>

using std::vector;
using std::list;

class ClusterStructure 
{
  public:
    ClusterStructure( const int numNodes, const int maxNumHeads);
    ~ClusterStructure();

    int GetChIdxByName( const int& nodeName );
    int GetChNameByName( const int& nodeName );

  private:
    const int   m_numNodes;
    const int   m_maxNumHeads;
    int         m_numHeads;

    list<list<int> >  m_listCluMember;
    vector<int>       m_vecSupport;
    vector<int>       m_vecHeadName;
    vector<int>       m_vecCHIdxForNodes;
    vector<int>       m_vecCHNameForNodes;
};
#endif
