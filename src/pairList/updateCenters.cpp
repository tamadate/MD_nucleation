//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------

void
MD:: updateInCenters(void){
	std::array<double,3> *mol=vars->position.data();
	std::array<double,3> *molv=vars->velocity.data();
	for (auto i : vars->gas_in){
		std::array<double,3> *xi = vars->positionAtom[i].data();
		std::array<double,3> *vi = vars->velocityAtom[i].data();
		double *mi = vars->massAtom[i].data();
		mol[i][0]=mol[i][1]=mol[i][2]=molv[i][0]=molv[i][1]=molv[i][2]=0;
		int Natom=vars->positionAtom[i].size();
		for (int ii=0; ii<Natom; ii++){
			mol[i][0]+=xi[ii][0]*mi[ii];
			mol[i][1]+=xi[ii][1]*mi[ii];
			mol[i][2]+=xi[ii][2]*mi[ii];
			molv[i][0]+=vi[ii][0]*mi[ii];
			molv[i][1]+=vi[ii][1]*mi[ii];
			molv[i][2]+=vi[ii][2]*mi[ii];
		}
		double massInv=1/vars->mass[i];
		//  averaged position
		mol[i][0]*=massInv;
		mol[i][1]*=massInv;
		mol[i][2]*=massInv;
		//  averaged velocity
		molv[i][0]*=massInv;
		molv[i][1]*=massInv;
		molv[i][2]*=massInv;
	}

	for (auto i : vars->vapor_in){
		std::array<double,3> *xi = vars->positionAtom[i].data();
		std::array<double,3> *vi = vars->velocityAtom[i].data();
		double *mi = vars->massAtom[i].data();
		mol[i][0]=mol[i][1]=mol[i][2]=molv[i][0]=molv[i][1]=molv[i][2]=0;
		int Natom=vars->positionAtom[i].size();
		for (int ii=0; ii<Natom; ii++){
			mol[i][0]+=xi[ii][0]*mi[ii];
			mol[i][1]+=xi[ii][1]*mi[ii];
			mol[i][2]+=xi[ii][2]*mi[ii];
			molv[i][0]+=vi[ii][0]*mi[ii];
			molv[i][1]+=vi[ii][1]*mi[ii];
			molv[i][2]+=vi[ii][2]*mi[ii];
		}
		double massInv=1/vars->mass[i];
		//  averaged position
		mol[i][0]*=massInv;
		mol[i][1]*=massInv;
		mol[i][2]*=massInv;
		//  averaged velocity
		molv[i][0]*=massInv;
		molv[i][1]*=massInv;
		molv[i][2]*=massInv;
	}
}