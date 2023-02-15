#pragma once
#include "variables.hpp"
#include "output/observer.hpp"
#include "potential/potential.hpp"
#include "PhysicalProp.hpp"
#include "MB.hpp"
//------------------------------------------------------------------------

class MD {
	private:


  public:

	double startTime;
	int Nth;
	int calculation_number;

	int gastype;	/*1:He, 2:Ar, 3:N2*/
	long int step_relax;
	long int Noftimestep;
	double p;
	double T;
	int Nof_around_gas;
	int Nof_around_vapor;
	int OBSERVE;
	double Ecoeff[3];

	double dt;
	double CUTOFF;
	double MARGIN;
	double ML2;
	double CL2;
	double d_size;
	double V;
	double margin_length;

	long int itime;

//	Pointers
	Variables *vars;
	Observer *obs;
	Physical *pp;
	FLAG *flags;
	MBdist *mbdist;
	MBdist *mbdistV;

//	Interactions
	std::vector<Potential*> InterInter;
	std::vector<Potential*> IntraInter;
	void setPotential(FLAG *flags,int mode);

//	velocity verlet
	void run(char** argv);
	void verlet(void);
	void update_position(void);
	void velocity_calculation(void);
	void updatePositionSHAKE(void);
	void shake(int ID);

//	pair list
	void check_pairlist(void);
	void make_pair(void);
	void update_vapor_in(void);
	void update_gas_in(void);
	void updateInCenters(void);
  	void makeDiatomicProp_in(int i);
	void makePolyatomicProp_in(int i);
	void periodic(void);	/*	periodic condition	*/
	void boundary_scaling_vapor_move(void);
	void boundary_scaling_gas_move(void);
	void boundary_scaling_ion_move(void);
	double preCenter[3];
	int loopPair, loop_update;	/*	current fixing time(loop) and update fixing time(loop) of out_gas for multi-timestep	*/

//	initialization
	void initialization_gas(void);
  	void initialization_vapor(void);
	void readCondFile(char* condfile);
	void setDefaultVariables(void);
	char atomFile[100];
	char vaporFile[100];

// related vapor sticking position
	int positionLogStep;
	std::vector<int> stickPositionList;

//	Thermostats
	void velocity_scaling(void);
	void nosehoover_ion(void);
	void nosehoover_zeta(void);
	double zeta;
	void setNVE(void);

/*other*/
	void export_dump(void);
	void export_dump_close(void);
	void output(void);
	void output_initial(void);
	void Ovin(int i);
	void Ovout(int i);
	void display(int output_ONOFF);
	char filepath[100];
	char vaporStickFile[100];
	char filepath_gyration[100];
	double gyration;
	void gyration_out(MD *md2);
	string gyration_path;
	double crsq;

	double totalPotential;
	MD(char* condfile,int calcNumber);
	~MD(void);
};




//------------------------------------------------------------------------
