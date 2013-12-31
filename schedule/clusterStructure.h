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

    void SetRecord( const vector<int>& invecheadname, 
        const list<list<int> >& inlistclumember);

    int GetChIdxByName( const int& nodeName ) const;
    int GetChNameByName( const int& nodeName ) const;

    const vector<int>&  GetVecHeadName() const { return m_vecHeadName; }
    const list<list<int> >& GetListCluMemeber() const{return m_listCluMember; }

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
