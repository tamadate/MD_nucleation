#pragma once
#include "variables.hpp"
#include "output/observer.hpp"
#include "potential/potential.hpp"
#include "PhysicalProp.hpp"
#include "pairList/boundary/MBdist.hpp"
//------------------------------------------------------------------------

class MD {
	private:


  public:

	double startTime;
	int Nth;
	int calculation_number;

	int gastype;	/*1:He, 2:Ar, 3:N2*/
	int vaportype;	/*1:MeOH, 2:H2O, 3:EtOH*/
	long int step_relax;
	long int step_repre;
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
	void setCondition(char* condfile);
	void readCondFile(char* condfile);

	long int itime;
	std::vector<long int> collisionFlagGas;
	std::vector<long int> collisionFlagVapor;
	std::vector<Potential*> InterInter;
	std::vector<Potential*> IntraInter;
	void setPotential(FLAG *flags,int mode);

	Variables *vars;
	Observer *obs;
	Physical *pp;
	FLAG *flags;
	MBdist *mbdist;
	MBdist *mbdistV;

//	vectors for pairlist
	double margin_length;

//	General functions
	void updateInCenters(void);

//	velocity verlet
	void run_diff(char** argv);
	void verlet(void);
	void update_position(void);
	void velocity_calculation(void);
	void update_position_constrained(void);
	void update_velocity_constrained(void);

	void updatePositionSHAKE(void);
	void recordOldPosition(int ID);
	void shake(int ID);

//	pair list
	void update_vapor_in(void);
	void update_gas_in(void);
	void make_pair(void);
	void check_pairlist(void);
  	void makeDiatomicProp_in(int i);
	void makePolyatomicProp_in(int i);

//	initialization
	void initialization_gas(void);
  	void initialization_vapor(void);

//	periodic
	void periodic(void);	/*	periodic condition for gas_in	*/
	void boundary_scaling_gas_move(void);
	void boundary_scaling_ion_move(void);
	void boundary_scaling_vapor_move(void);
	int loopPair, loop_update;	/*	current fixing time(loop) and update fixing time(loop) of out_gas for multi-timestep	*/
	double pre_ion[3];

//	analysis (calculating position and velocity of center of mass)
	double gyration;

//	export
	void export_dump(void);
	void export_dump_close(void);

// related vapor sticking position
	void positionLog(void);
	int positionLogStep;
	std::vector<int> stickPositionList;

/*other*/
	void output(void);
	void output_gas(void);
	void output_temp(double gastemp, double iontemp);
	void output_initial(void);
	void Ovin(int i);
	void Ovout(int i);
	void display(int output_ONOFF);
	char filepath[100];
	char atomFile[100];
	char vaporFile[100];
	char vaporStickFile[100];
	char filepath_gyration[100];
	void gyration_out(MD *md2);
	string gyration_path;
	double crsq;

	void velocity_scaling(void);
	void nosehoover_ion(void);
	void nosehoover_zeta(void);
	double zeta;
	void setNVE(void);


	double del2,CD2,rmin2;
	double totalPotential;


	MD(char* condfile,int calcNumber);
	~MD(void);
	void run(char** argv);
	int yesno;	/*	flag for collision or not collision	*/
};




//------------------------------------------------------------------------
