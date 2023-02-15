#pragma once
#include "../../variables.hpp"
#include "../basePotential.hpp"

class PotentialBorn : public Potential {
	private:
		string potName="Ion BMH";
	public:
		double bornCoeff[2][2][5];
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void makePair(){};
		PotentialBorn(){
			//difine Born-Mayer-Huggins conefficients for "only" NaCl
			//https://doi.org/10.1063/1.1522375
			//https://doi.org/10.1016/0022-3697(64)90160-X
			//0:A/rho, 1:6C, 2:8D, 3:sigma, 4:1/rho
			//NaNa, A=25.4435kJ/mol, C=101.1719kJ/mol, D=48.1771kJ/mol, sigma=2.340A, 1/rho=3.1546A-1
			//NaCl, A=20.3548kJ/mol, C=674.4793kJ/mol, D=837.077kJ/mol, sigma=2.755A, 1/rho=3.1546A-1
			//ClCl, A=15.2661kJ/mol, C=6985.6786kJ/mol, D=14031.5785kJ/mol, sigma=3.170A, 1/rho=3.1546A-1
			bornCoeff[0][0][0]=25.4435*3.1546/4.184;
			bornCoeff[0][1][0]=bornCoeff[1][0][0]=20.3548*3.1546/4.184;
			bornCoeff[1][1][0]=15.2661*3.1546/4.184;

			bornCoeff[0][0][1]=6*101.1719/4.184;
			bornCoeff[0][1][1]=bornCoeff[1][0][1]=6*674.4793/4.184;
			bornCoeff[1][1][1]=6*6985.6786/4.184;

			bornCoeff[0][0][2]=8*48.1771/4.184;
			bornCoeff[0][1][2]=bornCoeff[1][0][2]=8*837.077/4.184;
			bornCoeff[1][1][2]=8*14031.5785/4.184;

			bornCoeff[0][0][3]=2.340;
			bornCoeff[0][1][3]=bornCoeff[1][0][3]=2.755;
			bornCoeff[1][1][3]=3.170;

			bornCoeff[0][0][4]=3.1546;
			bornCoeff[0][1][4]=bornCoeff[1][0][4]=3.1546;
			bornCoeff[1][1][4]=3.1546;

		};
		~PotentialBorn(){};
};