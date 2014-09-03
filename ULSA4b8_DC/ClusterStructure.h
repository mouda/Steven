#include <vector>
#include <list>

using std::list;
using std::vector;
class ClusterStructure 
{
  public:
    ClusterStructure():m_numNode(0),m_maxNumCh(0){};
    ClusterStructure(const int& totalNodes, const int& maxChNum):
      m_numNode(totalNodes), m_maxNumCh(maxChNum) {};
    ~ClusterStructure();

    void SetVecHeadName( const vector<int>& inVecHeadName)
      {m_vecHeadName = inVecHeadName;}
    void SetListCluMember( const list<list<int> >& inListCluMember)
      {m_listCluMember = inListCluMember;}

    void SetRecord( const vector<int>& invecheadname, const list<list<int> >& inlistclumember);
      
    int GetChIdxByName( const int& nodeName );
    int GetChNameByName( const int& nodeName );
    const vector<int>&  GetVecHeadName(){ return m_vecHeadName; }
    const list<list<int> >& GetListCluMemeber(){return m_listCluMember; }

  private:
    const int         m_numNode;
    const int         m_maxNumCh;
    list<list<int> >  m_listCluMember;
    vector<int>       m_vecHeadName;
    vector<int>       m_vecCHIdxForNodes;
    vector<int>       m_vecCHNameForNodes;
};
