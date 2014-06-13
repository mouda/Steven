#include"TimeStamp.h"
#include <cassert>
void TimeStamp::returnRealWordTime(char timeBuf[], int sizeofBuffer){

    assert (sizeofBuffer>=32);
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (timeBuf,32,"%Y-%m-%d_%H-%M",timeinfo);
    return;
}
