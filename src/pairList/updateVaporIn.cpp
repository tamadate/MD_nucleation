#include "../md.hpp"


/////////////////////////////////////////////////////////////////////
/*
	make vapor-ion pair list
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_vapor_in(void){
// clear vars->gas_in, vars->gas_out and pair list of gas-ion

	vars->vapor_in.clear();
	vars->vapor_out.clear();
    double HL=d_size*0.5;

    random_device seed;
	mt19937 mt(seed());
    normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
    std::array<double,3> *x=vars->position.data();
    std::array<double,3> *v=vars->velocity.data();

	for (auto i : vars->group[2]){
		double dx = x[i][0] - x[0][0];
		double dy = x[i][1] - x[0][1];
		double dz = x[i][2] - x[0][2];
		double r2 = (dx * dx + dy * dy + dz * dz);
		if (r2 < RI2){
            //if inter-gas interaction flag is ON, stand flag_in ON
            //if flag_in ON, molecule push back to vars->gas_in
            //if flag_in OFF, molecule push back to vars->gas_out
            vars->vapor_in.push_back(i);
            makePolyatomicProp_in(i);
            vars->region[i]=true;
        }
        else{
            vars->vapor_out.push_back(i);
            if(vars->region[i]) Ovout(i);
            vars->region[i]=false;
            if (dx<-HL) x[i][0] += d_size;
            if (dy<-HL) x[i][1] += d_size;
            if (dz<-HL) x[i][2] += d_size;
            if (dx>HL) x[i][0] -= d_size;
            if (dy>HL) x[i][1] -= d_size;
            if (dz>HL) x[i][2] -= d_size;
            v[i][0] = distvapor(mt) *1e-5;
			v[i][1] = distvapor(mt) *1e-5;
			v[i][2] = distvapor(mt) *1e-5;
        }
	}
}

void
MD::makePolyatomicProp_in(int i){
	if(!vars->region[i]){
        random_device seed;
        mt19937 mt(seed());
        normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mvapor));
        uniform_real_distribution<double> r(0,1);

        vars->positionAtom[i]=vars->positionVapor;
        for (auto &atom : vars->velocityAtom[i]) atom[0]=0, atom[1]=0, atom[2]=0;
        for (auto &atom : vars->forceAtom[i]) atom[0]=0, atom[1]=0, atom[2]=0;

        double a=r(mt)*2*M_PI;
        double b=r(mt)*2*M_PI;
        double c=r(mt)*2*M_PI;
        for (auto &ag : vars->positionAtom[i]){
            vars->ROTATION(ag[0],ag[1],ag[2],a,b,c,ag[0],ag[1],ag[2]);
            ag[0]+=vars->position[i][0];
            ag[1]+=vars->position[i][1];
            ag[2]+=vars->position[i][2];
        }
        for (auto &ag : vars->velocityAtom[i]){
            vars->ROTATION(ag[0],ag[1],ag[2],a,b,c,ag[0],ag[1],ag[2]);
            ag[0]+=vars->velocity[i][0];
            ag[1]+=vars->velocity[i][1];
            ag[2]+=vars->velocity[i][2];
        }
    }
}

