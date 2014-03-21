#ifndef _GREEDYPHYSICAL_
#define _GREEDYPHYSICAL_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/directed_graph.hpp> 
#include <boost/bind.hpp>
#include <queue>


#include "scheduler.h"

using namespace boost;
typedef adjacency_list<vecS, vecS, directedS,
        /* Vertex properties */
        property<vertex_index_t, int>,
        /*  Edge properties */
        property<edge_weight_t, double> 
        > BglGraph;
typedef graph_traits<BglGraph>::vertex_descriptor     BglVertex;
typedef graph_traits<BglGraph>::edge_descriptor       BglEdge;
typedef graph_traits<BglGraph>::out_edge_iterator     OutEdgeIter;
typedef property_map<BglGraph, vertex_index_t>::type  BglVertexMap;
typedef property_map<BglGraph, edge_weight_t>::type   BglEdgeMap;
typedef BglGraph::edge_property_type                  BglEdgeWeight;
typedef graph_traits<BglGraph>::edge_iterator         EdgeIter;
typedef graph_traits<BglGraph>::vertex_iterator       VertexIter;
typedef std::pair<EdgeIter,EdgeIter>                  EdgePair; 
typedef std::pair<BglVertex,BglVertex>                BglVertexPair;




class GreedyPhysical: public Scheduler
{
  public:
    GreedyPhysical(
        const double txTime, 
        const double bandwidthKhz, 
        Map const * const, 
        CORRE_MA_OPE* , 
        ClusterStructure const * const);
    ~GreedyPhysical();
    void SetGaussianField(CORRE_MA_OPE* myGField) { m_ptrMatComputer = myGField;}
    double ScheduleOneSlot( std::vector<int>& );
    bool ScheduleOneSlot( std::vector<int>& , const std::vector<double>& );
    string PrintSelf(){ return m_type; }
    double GetTxTimePerSlot() const { return m_txTimePerSlot; }
    
  private:
    double GetInterferenceNumber(const int source, const int target );
    double GetInterferenceNumber(const int source, const int target, const std::vector<BglEdge>& );
    int  GreedySelectOneNode( const std::vector<BglEdge>& );
    double GetRxPower(const int source, const int target);

    int                             m_numNodes;
    int                             m_numMaxHeads;
    const double                    m_txTimePerSlot;
    const double                    m_bandwidthKhz;
    double                          m_maxPower;
    Map const * const               m_ptrMap;
    ClusterStructure const * const  m_ptrCS;
    CORRE_MA_OPE*                   m_ptrMatComputer;
    std::vector<double>             m_vecNodePower;
    std::vector<int>                m_vecSched; /* record the scheduled node */
    const string                    m_type;

    /* Graph */

    BglGraph                        m_commGraph;
    BglGraph                        m_conflictGraph;
}; 
#endif
