#ifndef _SIMULATOR_
#define _SIMULATOR_ 
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "imageMap.h"
#include "slot.h"
#include "clusterStructure.h"
#include "scheduler.h"
#include "imageSource.h"
#include "../lib/fileHandler.h"

using std::string;
using std::pair;
using std::cout;
using std::endl;
using std::cerr;
using std::fstream;
using std::ios;
using std::setfill;
using std::setw;
using std::stringstream;
using std::make_pair;

class Simulator
{
  public:
    Simulator(ImageMap* myMap,ClusterStructure* myCS, Scheduler* myScheduler, ImageSource* myField, 
        const string&,
        const string&,
        const string&,
        const string&, 
        const string&,
        const string&);
    ~Simulator();

    void SequentialRun(int tier2NumSlot);
    void SequentialFree();
    void Run();

    void WriteCS( const string& );
    void WriteEntropy();
    void WriteMSE();
    void WriteSolution();
    void WriteSupport();
    void WriteTotalEntropy();
    void WriteVecPower();


    bool SelfCheck();

  private:
    double GetTotalEntropy(const std::vector<int>& vecSupport) const;
    double Get1stSlotEntropy(const std::vector<int>& vecSupport) const;
    void Print(const std::vector<int>& );
    template < class T>
    string VecToString( const std::vector<T>& );
    Slot* GetNextSlot(Slot*, std::vector<int>& );

    ImageMap*                m_ptrImageMap;
    ClusterStructure*   m_ptrCS;
    Scheduler*          m_ptrSched;
    ImageSource*       m_ptrImageSource;
    std::vector<int>*   m_vecSupport;

    std::vector<int>    m_vecTotal;

    /* event driven simulation data structure */
    std::list<Slot*>    m_listSlot;
    Slot*               m_ptrSlotHead;

    /* output  */
    FileHandler         m_entropyFHandler;
    FileHandler         m_MSEFHandler;
    FileHandler         m_solutionFHandler;
    FileHandler         m_supportNumFHandler;
    FileHandler         m_totalEntropyFHandler;
    FileHandler         m_vecPowerFHandler;
};
#endif
