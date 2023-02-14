//------------------------------------------------------------------------
#include "potentialSW.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Calculate the intra-atomic interaction (Stillinger-Weber)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialSW::compute(Variables *vars, FLAG *flags) {
	std::array<double,3> *x = vars->positionAtom[0].data();
	std::array<double,3> *f = vars->forceAtom[0].data();
	const int is = vars->positionAtom[0].size();
	U=0;
	Pressure=0;
	virial=0;

	for(auto &a : pairs){
		int i=a.i;
		neighshort.clear();
//	repulsive forces d(fCfR)/dr
		for(auto j : a.j){
			delr[0]=x[i][0]-x[j][0];
			delr[1]=x[i][1]-x[j][1];
			delr[2]=x[i][2]-x[j][2];
			//adjust_periodic(delr[0], delr[1], delr[2]);
			double rsq = (delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2]);
			if(rsq>cut2) continue;
			double force_pair=twobody(rsq)*0.5;
			f[i][0]+=force_pair*delr[0];
			f[i][1]+=force_pair*delr[1];
			f[i][2]+=force_pair*delr[2];
			f[j][0]-=force_pair*delr[0];
			f[j][1]-=force_pair*delr[1];
			f[j][2]-=force_pair*delr[2];

			neighshort.push_back(j);

//			virial+=force_pair*r*r;
		}

//	three body interaction
		for(auto j : neighshort){
			delr[0]=x[i][0]-x[j][0];
			delr[1]=x[i][1]-x[j][1];
			delr[2]=x[i][2]-x[j][2];
			//adjust_periodic(delr[0], delr[1], delr[2]);
			double rsq = (delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2]);

//			virial+=force_pair*r*r;

			for (auto k : neighshort){
				if(k==j) continue;
				delr2[0] = x[i][0]-x[k][0];
				delr2[1] = x[i][1]-x[k][1];
				delr2[2] = x[i][2]-x[k][2];
			//	adjust_periodic(delr2[0], delr2[1], delr2[2]);
				double rsq2 = (delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2]);
				if(rsq2>cut2) continue;
				double fi[3],fj[3],fk[3];
				threebody(rsq,rsq2,delr,delr2,fj,fk);
				f[i][0]+=(fj[0]+fk[0])*0.5;
				f[i][1]+=(fj[1]+fk[1])*0.5;
				f[i][2]+=(fj[2]+fk[2])*0.5;
				f[j][0]-=fj[0]*0.5;
				f[j][1]-=fj[1]*0.5;
				f[j][2]-=fj[2]*0.5;
				f[k][0]-=fk[0]*0.5;
				f[k][1]-=fk[1]*0.5;
				f[k][2]-=fk[2]*0.5;
			}
		}
	}

	for(int i=0;i<is;i++){
/*	FILE*f=fopen("force.dat", "a");
	fprintf(f, "%f\t%f\t%f\t%e\t%e\t%e\t%f\n", ions[i].qx, ions[i].qy, ions[i].qz, ions[i].fx, ions[i].fy, ions[i].fz, energy);
	fclose(f);*/
		f[i][0]*=23.06;
		f[i][1]*=23.06;
		f[i][2]*=23.06;
	}
	U*=0.5;
//	vars->Potential+=energy*23.06;
//	vars->virial+=virial*23.06;

}

double
PotentialSW::computeVirial(Variables *vars) {
/*	double virial=0;
	for(auto &a : vars->ions){
		virial+=a.fx*a.qx;
		virial+=a.fy*a.qy;
		virial+=a.fz*a.qz;
	}*/
	return virial*23.061;
}



double
PotentialSW::twobody(double rsq){
	double r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;

	r = sqrt(rsq);
	rinvsq = 1.0/rsq;
	rp = pow(r,-powerp);
	rq = pow(r,-powerq);
	rainv = 1.0 / (r - cut);
	rainvsq = rainv*rainv*r;
	expsrainv = exp(sigma * rainv);
	U+=(c5*rp-c6*rq)*expsrainv;
	return (c1*rp - c2*rq + (c3*rp -c4*rq) * rainvsq) * expsrainv * rinvsq;
}


void
PotentialSW::threebody(double rsq1, double rsq2,double *rij, double *rik,double *fj, double *fk)
{
	double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
	double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
	double rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2;
	double facang,facang12,csfacang,csfac1,csfac2;

	r1 = sqrt(rsq1);
	rinvsq1 = 1.0/rsq1;
	rainv1 = 1.0/(r1 - cut);
	gsrainv1 = sigma_gamma * rainv1;
	gsrainvsq1 = gsrainv1*rainv1/r1;
	expgsrainv1 = exp(gsrainv1);

	r2 = sqrt(rsq2);
	rinvsq2 = 1.0/rsq2;
	rainv2 = 1.0/(r2 - cut);
	gsrainv2 = sigma_gamma * rainv2;
	gsrainvsq2 = gsrainv2*rainv2/r2;
	expgsrainv2 = exp(gsrainv2);

	rinv12 = 1.0/(r1*r2);
	cs = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) * rinv12;
	delcs = cs - costheta0;
	delcssq = delcs*delcs;

	facexp = expgsrainv1*expgsrainv2;

	facrad = lambda_epsilon * facexp*delcssq;
	frad1 = facrad*gsrainvsq1;
	frad2 = facrad*gsrainvsq2;
	facang = lambda_epsilon2 * facexp*delcs;
	facang12 = rinv12*facang;
	csfacang = cs*facang;

	csfac1 = rinvsq1*csfacang;

	fj[0] = rij[0]*(frad1+csfac1)-rik[0]*facang12;
	fj[1] = rij[1]*(frad1+csfac1)-rik[1]*facang12;
	fj[2] = rij[2]*(frad1+csfac1)-rik[2]*facang12;

	csfac2 = rinvsq2*csfacang;

	fk[0] = rik[0]*(frad2+csfac2)-rij[0]*facang12;
	fk[1] = rik[1]*(frad2+csfac2)-rij[1]*facang12;
	fk[2] = rik[2]*(frad2+csfac2)-rij[2]*facang12;

	U+=facrad;
}

void
PotentialSW::check_pairlist(Variables *vars){
	loop_t++;
	if(loop_t>loop_update_t) make_pair(vars);
}

void
PotentialSW::make_pair(Variables *vars){
	pairs.clear();
	std::array<double,3> *x = vars->positionAtom[0].data();
	int is=vars->positionAtom[0].size();
	for (int i=0; i<is; i++){
		Pair_many p;
		p.i=i;
		vector<int> js;
		for (int j=0; j<is; j++){
			if (i==j) continue;
			double dx = x[i][0] - x[j][0];
			double dy = x[i][1] - x[j][1];
			double dz = x[i][2] - x[j][2];
	//		adjust_periodic(dx, dy, dz);
			double rsq = (dx * dx + dy * dy + dz * dz);
			if (rsq < sw_ML2) js.push_back(j);
		}
		p.j=js;
		pairs.push_back(p);
	}
	loop_t=0;
}
