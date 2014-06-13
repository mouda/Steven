#include "timesearch.h"
#include <cmath>
timesearch::timesearch(){
	lastResors=0;
	lastFidR=0;
	lastDelta=0;
	upperResors=0;
	lowerResors=0;
	firstAnsFound=false;
}

double timesearch::update2ndResors_GivenFidelity(double curResors, double curFidR, double targetFidR){

   double nextResors=0;
   double delta=0;
   double unit=0.2;
   double dFidR=-1;
   if(lastDelta!=0){
	dFidR=curFidR-lastFidR;
	double dResource=curResors-lastResors;
	unit=dResource/dFidR;
	delta=(targetFidR-curFidR)/dFidR*unit;
	nextResors=((curResors-delta>0)?curResors+delta:curResors/2);
   }
   else {
	lastDelta=((curFidR<targetFidR)?unit:-unit);;
	lastFidR=curFidR;
	lastResors=curResors;
	nextResors=((curFidR<targetFidR)?curResors+unit:curResors-unit);
   }
   cout<<"update unit "<<unit<<"per "<<dFidR<<endl;
   cout<<"Next Resource(ms)="<<nextResors<<endl;

   return nextResors;
}
double timesearch::update2ndResors_binerySearch_GivenFidelity(double curResors, double curFidR, double targetFidR,bool &curSolFound){

    double nextResors;

    if(curFidR>targetFidR){
        nextResors=(lowerResors+curResors)/2;
        upperResors=curResors;
        curSolFound=false;
        firstAnsFound = true;
        return nextResors;
    }

   else  if (!firstAnsFound){
        lowerResors=curResors;
        return 3*curResors;
    }
    else if (curFidR<targetFidR&&curSolFound){
        return curResors;

    }
    else{
        nextResors=(upperResors+curResors)/2;
        lowerResors=curResors;
        return nextResors;
    }

}


