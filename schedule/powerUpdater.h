#ifndef _POWERUPDATER_
#define _POWERUPDATER_

#include <vector>

class PowerUpdater
{
  public:
    PowerUpdater();
    ~PowerUpdater();
  private:

    void Init();
    void UpdatePowerVector( std::vector<int>& );

};
#endif
