#pragma once
#include "../variables.hpp"

//------------------------------------------------------------------------
class Potential {
	private:
	public:
		string potName="potential";
		std::vector<Pair> pairList;
		double ML2;
		virtual void printName(void) {cout<<potName<<endl;}
		virtual void compute(Variables *vars, FLAG *flags){};
		virtual void makePair(Variables *vars){};
		virtual void init(Variables *vars){};
		Potential(){};
		~Potential(){};
};
