#include "potential.hpp"
/*########################################################################################

-----compute intramolecular interaction-----

#######################################################################################*/
#define TOLERANCE 0.05
#define SMALL     0.001
#define SMALLER   0.00001
/**********************************Force calculation******************************************/


void
PotentialGasIntra::compute(Variables *vars, FLAG *flags){
	vars->times.tgas-=omp_get_wtime();
	for(auto &I : vars->gas_in){
		std::array<double,3> *x = vars->positionAtom[I].data();
		std::array<double,3> *f = vars->forceAtom[I].data();
		for(auto &p : vars->bonds[I]){
			int i=p[0];
			int j=p[1];
			int type=p[2];
			double dx = x[i][0] - x[j][0];
			double dy = x[i][1] - x[j][1];
			double dz = x[i][2] - x[j][2];
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r = sqrt(rsq);
			double dr= r - vars->btypes[type].coeff[1];
			double rk = vars->btypes[type].coeff[0] * dr;
			double force_bond_harmonic = -2.0*rk/r;
			f[i][0] += force_bond_harmonic * dx;
			f[i][1] += force_bond_harmonic * dy;
			f[i][2] += force_bond_harmonic * dz;
			f[j][0] -= force_bond_harmonic * dx;
			f[j][1] -= force_bond_harmonic * dy;
			f[j][2] -= force_bond_harmonic * dz;
			if(flags->eflag) vars->Utotal.Ugas+=rk*dr;
		}
	}
	vars->times.tgas+=omp_get_wtime();
}


void
PotentialVaporIntra::compute(Variables *vars, FLAG *flags) {
	double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2, C, Cs, dtheta, tk, a, a11, a12, a22, f1[3], f3[3];

	double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
	double edihedral,ff2[3],ff4[3],ff1[3],ff3[3];
	double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
	double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
	double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
	double c,s,p_,sx2,sy2,sz2, m;

	double sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
	double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
	double c2mag,sc1,sc2,s1,s12,p,pd;
	double a33,a13,a23;
	double s2,cx,cy,cz,cmag,dx,phi,si,siinv,sin2;
	vars->times.tvap-=omp_get_wtime();
	for(auto &I : vars->vapor_in){
		std::array<double,3> *x = vars->positionAtom[I].data();
		std::array<double,3> *f = vars->forceAtom[I].data();
		for(auto &p : vars->bonds[I]){
			int i=p[0];
			int j=p[1];
			int type=p[2];
			double dx = x[i][0] - x[j][0];
			double dy = x[i][1] - x[j][1];
			double dz = x[i][2] - x[j][2];
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r = sqrt(rsq);
			double dr= r - vars->btypes[type].coeff[1];
			double rk = vars->btypes[type].coeff[0] * dr;
			double force_bond_harmonic = -2.0*rk/r;
			f[i][0] += force_bond_harmonic * dx;
			f[i][1] += force_bond_harmonic * dy;
			f[i][2] += force_bond_harmonic * dz;
			f[j][0] -= force_bond_harmonic * dx;
			f[j][1] -= force_bond_harmonic * dy;
			f[j][2] -= force_bond_harmonic * dz;
			if(flags->eflag) vars->Utotal.Ugas+=rk*dr;
		}

		for(auto &p : vars->angles[I]){
			int i=p[0];
			int j=p[1];
			int k=p[2];
			int type=p[3];
			double dx1 = x[i][0] - x[j][0];
			double dy1 = x[i][1] - x[j][1];
			double dz1 = x[i][2] - x[j][2];
			double rsq1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
			double r1 = sqrt(rsq1);
			double dx2 = x[k][0] - x[j][0];
			double dy2 = x[k][1] - x[j][1];
			double dz2 = x[k][2] - x[j][2];
			double rsq2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
			double r2 = sqrt(rsq2);
			double C = dx1*dx2 + dy1*dy2 + dz1*dz2;
			C /= r1*r2;
			double Cs = 1/(sqrt(1.0-C*C));
			double dtheta = acos(C) - vars->ctypes[type].coeff[1];
			double tk = vars->ctypes[type].coeff[0] * dtheta;
			double a = -2.0 * tk * Cs;
			double a11 = a*C / rsq1;
			double a12 = -a / (r1*r2);
			double a22 = a*C / rsq2;
			f1[0] = a11*dx1 + a12*dx2;
			f1[1] = a11*dy1 + a12*dy2;
			f1[2] = a11*dz1 + a12*dz2;
			f3[0] = a22*dx2 + a12*dx1;
			f3[1] = a22*dy2 + a12*dy1;
			f3[2] = a22*dz2 + a12*dz1;
			f[i][0] += f1[0];
			f[i][1] += f1[1];
			f[i][2] += f1[2];
			f[j][0] -= f1[0] + f3[0];
			f[j][1] -= f1[1] + f3[1];
			f[j][2] -= f1[2] + f3[2];
			f[k][0] += f3[0];
			f[k][1] += f3[1];
			f[k][2] += f3[2];
			if(flags->eflag) vars->Utotal.Ugas+=tk*dtheta;
		}


		for (auto &p : vars->dihedrals[I]) {
			int i=p[0];
			int j=p[1];
			int k=p[2];
			int l=p[3];
			int type=p[4];
			// 1st bond
			double vb1x = x[i][0] - x[j][0];
			double vb1y = x[i][1] - x[j][1];
			double vb1z = x[i][2] - x[j][2];
			// 2nd bond
			double vb2x = x[k][0] - x[j][0];
			double vb2y = x[k][1] - x[j][1];
			double vb2z = x[k][2] - x[j][2];
			double vb2xm = -vb2x;
			double vb2ym = -vb2y;
			double vb2zm = -vb2z;
			// 3rd bond
			double vb3x = x[l][0] - x[k][0];
			double vb3y = x[l][1] - x[k][1];
			double vb3z = x[l][2] - x[k][2];
			
			double ax = vb1y*vb2zm - vb1z*vb2ym;
			double ay = vb1z*vb2xm - vb1x*vb2zm;
			double az = vb1x*vb2ym - vb1y*vb2xm;
			double bx = vb3y*vb2zm - vb3z*vb2ym;
			double by = vb3z*vb2xm - vb3x*vb2zm;
			double bz = vb3x*vb2ym - vb3y*vb2xm;
			double rasq = ax*ax + ay*ay + az*az;
			double rbsq = bx*bx + by*by + bz*bz;
			double rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
			double rg = sqrt(rgsq);
			double rginv=0.0;
			double ra2inv=0.0;
			double rb2inv=0.0;
			if (rg > 0) rginv = 1.0/rg;
			if (rasq > 0) ra2inv = 1.0/rasq;
			if (rbsq > 0) rb2inv = 1.0/rbsq;
			double rabinv = sqrt(ra2inv*rb2inv);
			double c = (ax*bx + ay*by + az*bz)*rabinv;
			double s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

			double df = 0.0;

			double p_,df1,ddf1;
			for(int JJ=0;JJ<vars->dtypes[type].multi;JJ++){
				int JJ5=JJ*5;
				int JJ5_1=JJ5+1;
				int JJ5_3=JJ5+3;
				int JJ5_4=JJ5+4;
				p_=1.0;
				ddf1=df1=0.0;
				for (int loop=0; loop<vars->dtypes[type].coeff[JJ5_1]; loop++){
					ddf1 = p_*c - df1*s;
					df1 = p_*s + df1*c;
					p_ = ddf1;
				}
				p_ = p_*vars->dtypes[type].coeff[JJ5_3] + df1*vars->dtypes[type].coeff[JJ5_4];
				df1 = df1*vars->dtypes[type].coeff[JJ5_3] - ddf1*vars->dtypes[type].coeff[JJ5_4];
				df1 *= -vars->dtypes[type].coeff[JJ5_1];
				p_ += 1.0;
				if(vars->dtypes[type].coeff[1]==0){
					p_=1.0+vars->dtypes[type].coeff[JJ5_3];
					df1=0.0;
				}
				df += (-vars->dtypes[type].coeff[JJ5] * df1);
				if (flags->eflag) vars->Utotal.Uion+= vars->dtypes[type].coeff[JJ5] * p_;
			}

			// cout<<df<<endl;
			double fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
			double hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
			double fga = fg*ra2inv*rginv;
			double hgb = hg*rb2inv*rginv;
			double gaa = -ra2inv*rg;
			double gbb = rb2inv*rg;
			double dtfx = gaa*ax;
			double dtfy = gaa*ay;
			double dtfz = gaa*az;
			double dtgx = fga*ax - hgb*bx;
			double dtgy = fga*ay - hgb*by;
			double dtgz = fga*az - hgb*bz;
			double dthx = gbb*bx;
			double dthy = gbb*by;
			double dthz = gbb*bz;
			double sx2 = df*dtgx;
			double sy2 = df*dtgy;
			double sz2 = df*dtgz;
			ff1[0] = df*dtfx;
			ff1[1] = df*dtfy;
			ff1[2] = df*dtfz;
			ff2[0] = sx2 - ff1[0];
			ff2[1] = sy2 - ff1[1];
			ff2[2] = sz2 - ff1[2];
			ff4[0] = df*dthx;
			ff4[1] = df*dthy;
			ff4[2] = df*dthz;
			ff3[0] = -sx2 - ff4[0];
			ff3[1] = -sy2 - ff4[1];
			ff3[2] = -sz2 - ff4[2];

			f[i][0] += ff1[0];
			f[i][1] += ff1[1];
			f[i][2] += ff1[2];
			f[j][0] += ff2[0];
			f[j][1] += ff2[1];
			f[j][2] += ff2[2];
			f[k][0] += ff3[0];
			f[k][1] += ff3[1];
			f[k][2] += ff3[2];
			f[l][0] += ff4[0];
			f[l][1] += ff4[1];
			f[l][2] += ff4[2];

		}
	}
	vars->times.tvap+=omp_get_wtime();
}
