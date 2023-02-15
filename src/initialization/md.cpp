#include "../md.hpp"

/////////////////////////////////////////////////////////////////////
/*
	constructor
*/
/////////////////////////////////////////////////////////////////////
MD::MD(char* condfile, int calcNumber) {
	startTime=omp_get_wtime();
	calculation_number =  calcNumber;
	vars = new Variables();
	obs = new Observer();
	pp = new Physical();
	flags = new FLAG();
	mbdist = new MBdist();
	mbdistV = new MBdist();

	setDefaultVariables();
	readCondFile(condfile);
	pp->readIonProp(atomFile);
	pp->readVaporProp(vaporFile);
	pp->setPhysicalProp(gastype,T,p);
	output_initial();

	vars->readIonFile(atomFile);
	vars->ionInitial(T);
	initialization_gas();	//Set initial positions & velocities for gas
	vars->readVaporFile(vaporFile);
	initialization_vapor();	//Set initial positions & velocities for vapor

	vars->setCrossPotentials();
	setPotential(flags,1);

	make_pair();
	margin_length = MARGIN;
	vars->tzero();

	mbdist -> makeWeightedMB(pp->cgas,pp->mgas,T);
	mbdistV -> makeWeightedMB(pp->cvapor,pp->mvapor,T);
}

/////////////////////////////////////////////////////////////////////
/*
	destructor
*/
/////////////////////////////////////////////////////////////////////
MD::~MD(void) {
	delete vars;
	delete obs;
	delete pp;
	delete flags;
}
