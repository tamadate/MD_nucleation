#include "../md.hpp"


/////////////////////////////////////////////////////////////////////
/*
	make gas-ion pair list
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_gas_in(void){
// clear vars->gas_in, vars->gas_out and pair list of gas-ion
	vars->gas_in.clear();
	vars->gas_out.clear();
    double HL=d_size*0.5;

    random_device seed;
	mt19937 mt(seed());
    normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
    std::array<double,3> *x=vars->position.data();
    std::array<double,3> *v=vars->velocity.data();
	for (auto i : vars->group[1]){
		double dx = x[i][0] - x[0][0];
		double dy = x[i][1] - x[0][1];
		double dz = x[i][2] - x[0][2];
		double r2 = (dx * dx + dy * dy + dz * dz);
		if (r2 < RI2){
            //if inter-gas interaction flag is ON, stand flag_in ON
            //if flag_in ON, molecule push back to vars->gas_in
            //if flag_in OFF, molecule push back to vars->gas_out
            vars->gas_in.push_back(i);
            makeDiatomicProp_in(i);
            vars->region[i]=true;
        }
        else{
            vars->gas_out.push_back(i);
            vars->region[i]=false;
            if (dx<-HL) x[i][0] += d_size;
            if (dy<-HL) x[i][1] += d_size;
            if (dz<-HL) x[i][2] += d_size;
            if (dx>HL) x[i][0] -= d_size;
            if (dy>HL) x[i][1] -= d_size;
            if (dz>HL) x[i][2] -= d_size;
            v[i][0] = distgas(mt) *1e-5;
			v[i][1] = distgas(mt) *1e-5;
			v[i][2] = distgas(mt) *1e-5;
        }
	}
}

void
MD::makeDiatomicProp_in(int i){
    // if previous region is out
  	if(!vars->region[i]){
        random_device seed;
        mt19937 mt(seed());
        normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
		uniform_real_distribution<double> r(0,1);

        vars->positionAtom[i] = vars->positionGas;      // reset gas position to origin
        // reset gas velocity = 0
		std::array<double,3> *vi = vars->velocityAtom[i].data();
        int Natom=vars->massAtom[i].size();
        for(int ii=0; ii<Natom; ii++){
            vi[ii][0]=0;
            vi[ii][1]=0;
            vi[ii][2]=0;
        }
        // give internal velocity (rotation)
        if(gastype==2){
            double vy_rot = distgas(mt) *1e-5;
		    double vz_rot = distgas(mt) *1e-5;
            vi[0][0]=0;
            vi[0][1]=vy_rot;
            vi[0][2]=vy_rot;
            vi[1][0]=0;
            vi[1][1]=-vy_rot;
            vi[1][2]=-vy_rot;
		}

        // randomly rotate gas
		double a=r(mt)*2*M_PI;
		double b=r(mt)*2*M_PI;
		double c=r(mt)*2*M_PI;
		for (auto &ag : vars->positionAtom[i]){
	        vars->ROTATION(ag[0],ag[1],ag[2],a,b,c,ag[0],ag[1],ag[1]);
            // set center of mass position
            ag[0]+=vars->position[i][0];
			ag[1]+=vars->position[i][1];
			ag[2]+=vars->position[i][2];
		}
        for (auto &ag : vars->velocityAtom[i]){
	        vars->ROTATION(ag[0],ag[1],ag[2],a,b,c,ag[0],ag[1],ag[1]);
            // set translational velocity
            ag[0]+=vars->velocity[i][0];
			ag[1]+=vars->velocity[i][1];
			ag[2]+=vars->velocity[i][2];
		}
  }
}