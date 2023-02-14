#pragma once
#include "../variables.hpp"
#include "basePotential.hpp"
#include "StillingerWeber/potentialSW.hpp"
#include "Tersoff/potentialTersoff.hpp"


class PotentialGasIon : public Potential {
	private:
		string potName="LJ gas-ion";
	public:
		double ML2;
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void makePair(Variables *vars);
		PotentialGasIon(double ml2){ML2=ml2;};
		~PotentialGasIon(){};
};

class PotentialVaporVapor : public Potential {
	private:
		string potName="LJ vapor-vapor";
	public:
		double ML2;
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void makePair(Variables *vars);
		PotentialVaporVapor(double ml2){ML2=ml2;};
		~PotentialVaporVapor(){};
};

class PotentialVaporGas : public Potential {
	private:
		string potName="LJ vapor-gas";
	public:
		double ML2;
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void makePair(Variables *vars);
		PotentialVaporGas(double ml2){ML2=ml2;};
		~PotentialVaporGas(){};
};

class PotentialVaporIon : public Potential {
	private:
		string potName="LJ vapor-ion";
	public:
		double ML2;
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void makePair(Variables *vars){};
		PotentialVaporIon(double ml2){ML2=ml2;};
		~PotentialVaporIon(){};
};

class PotentialGasGas : public Potential {
	private:
		string potName="LJ gas-gas";
	public:
		double ML2;
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void makePair(Variables *vars);
		PotentialGasGas(double ml2){ML2=ml2;};
		~PotentialGasGas(){};
};


class PotentialAMBER : public Potential {
	private:
		string potName="AMBER ion";
		std::vector<Pair> longPair;
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		void computeLong(Variables *vars, FLAG *flags);
		void computeBond(Variables *vars, FLAG *flags);
		void computeAngle(Variables *vars, FLAG *flags);
		void computeDihedral(Variables *vars, FLAG *flags);
		void initialAMBER(Variables *vars, FLAG *flags);
		PotentialAMBER(Variables *vars, FLAG *flags){
			initialAMBER(vars,flags);
		};
		~PotentialAMBER(){};
};


class PotentialGasIntra : public Potential {
	private:
		string potName="Gas intra";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialGasIntra(){};
		~PotentialGasIntra(){};
};

class PotentialBorn : public Potential {
	private:
		string potName="Ion BMH";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialBorn(){};
		~PotentialBorn(){};
};

class PotentialEfield : public Potential {
	private:
		string potName="Ion E-filed";
	public:
		void printName(void) {cout<<potName<<endl;}
		double Ecoeff[3];
		void compute(Variables *vars, FLAG *flags);
		void makePair(Variables *vars){};
		PotentialEfield(double Ex, double Ey, double Ez){
			Ecoeff[0]=Ex;
			Ecoeff[1]=Ey;
			Ecoeff[2]=Ez;
		};
		~PotentialEfield(){};
};

class PotentialIonDipole : public Potential {
	private:
		string potName="Induced dipole ion-gas";
	public:
		void printName(void) {cout<<potName<<endl;}
		double alphagas;
		double zion;
		void compute(Variables *vars, FLAG *flags);
		PotentialIonDipole(){};
		~PotentialIonDipole(){};
};

class PotentialVaporIntra : public Potential {
	private:
		string potName="AMBER vapor";
	public:
		void printName(void) {cout<<potName<<endl;}
		void compute(Variables *vars, FLAG *flags);
		PotentialVaporIntra(){};
		~PotentialVaporIntra(){};
};


//------------------------------------------------------------------------
