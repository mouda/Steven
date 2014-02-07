#ifndef _EVENT_
#define _EVENT_

class Event
{
  public:
    Event();
    ~Event();
    double GetSpatialCorrFac() const { return m_spatialCorrFac;}
    double GetTimeCorrFac() const { return m_timeCorrFac;}

  private:
    double  m_spatialCorrFac;
    double  m_timeCorrFac;

};

#endif
