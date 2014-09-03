#include<iostream>
using namespace std;
class timesearch{
	public:
		timesearch();
		double update2ndResors_GivenFidelity(double curResors, double curFidR, double targetFidR);
  		double update2ndResors_binerySearch_GivenFidelity(double curResors, double curFidR, double targetFidR,bool & cursolFound);

  		double lastResors;
		double lastFidR;
		double lastDelta;

		bool firstAnsFound;

		double upperResors;
		double lowerResors;
};
