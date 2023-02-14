#include "../md.hpp"
/*########################################################################################

-----Output-----

#######################################################################################*/

/**********************************initialization******************************************/
void
MD::output_initial(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number));
	FILE*f=fopen(filepath, "w");
	fclose(f);
	/*sprintf(filepath, "%d_%d_relax.dat", int(T), int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "gas_%d_%d.dat", int(T), int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);*/
	sprintf(filepath, "vapor_collision_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	/*sprintf(filepath, "gas_collision_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);*/
	sprintf(filepath, "K_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "U_%d.dat", int(calculation_number));
	f=fopen(filepath, "w");
	fclose(f);
}


void
MD::output(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number));
	FILE*f=fopen(filepath, "a");
	std::array<double,3> *x = vars->position.data();
	std::array<double,3> *v = vars->velocity.data();
	fprintf(f, "%e %e %e %e %e %e %e\n", vars->time/1e6, x[0][0], x[0][1], x[0][2], v[0][0], v[0][1], v[0][2]);
	fclose(f);
}

void
MD::output_gas(void){
	sprintf(filepath, "gas_%d_%d.dat", int(T), int(calculation_number));
	std::array<double,3> *x = vars->position.data();
	std::array<double,3> *v = vars->velocity.data();
	double X[3], V[3], Mass;
	X[0]=X[1]=X[2]=V[0]=V[1]=V[2]=Mass=0;
	for (auto &i : vars->group[1]){
		X[0]+=x[i][0]*vars->mass[i];
		X[1]+=x[i][1]*vars->mass[i];
		X[2]+=x[i][2]*vars->mass[i];
		V[0]+=v[i][0]*vars->mass[i];
		V[1]+=v[i][1]*vars->mass[i];
		V[2]+=v[i][2]*vars->mass[i];
		Mass+=vars->mass[i];
	}
	X[0]/=Mass;
	X[1]/=Mass;
	X[2]/=Mass;
	V[0]/=Mass;
	V[1]/=Mass;
	V[2]/=Mass;

	FILE*f=fopen(filepath, "a");
	fprintf(f, "%e %e %e %e %e %e %e\n", vars->time/1e6, X[0], X[1], X[2], V[0], V[1], V[2]);
	fclose(f);
}

void
MD::Ovin(int i){
	sprintf(filepath, "vapor_in_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	std::array<double,3> &x = vars->position[0];
	std::array<double,3> &v = vars->velocity[0];
	std::array<double,3> &xi = vars->position[i];
	std::array<double,3> &vi = vars->velocity[i];
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,xi[0]-x[0],xi[1]-x[1],xi[2]-x[2],vi[0]-v[0],vi[1]-v[1],vi[2]-v[2]);
	fclose(f);
}

void
MD::Ovout(int i){
	sprintf(filepath, "vapor_out_%d.dat", int(calculation_number));
	FILE*f=fopen(filepath, "a");
	std::array<double,3> &x = vars->position[0];
	std::array<double,3> &v = vars->velocity[0];
	std::array<double,3> &xi = vars->position[i];
	std::array<double,3> &vi = vars->velocity[i];
	fprintf(f, "%d %e %e %e %e %e %e %e\n", i,itime*dt,xi[0]-x[0],xi[1]-x[1],xi[2]-x[2],vi[0]-v[0],vi[1]-v[1],vi[2]-v[2]);
	fclose(f);
}

void
MD::output_temp(double gastemp, double iontemp){
	sprintf(filepath, "%d_%d_relax.dat", int(T), int(calculation_number));
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%f\t%f\t%f\n", vars->time/1e6, gastemp, iontemp);
	fclose(f);
}


void
MD::display(int output_ONOFF){
	obs->computeProps(vars);
    double virial=0;//ters->compute_tersoff_virial(vars)/3.0/V*Cpress;
    double gaspress=(kb*Nof_around_gas*obs->T_g + vars->totalVirial/3.0*6.95e-21)/(V*1e-30);
    double U = vars->Usum();
	double K = obs->Kion+obs->K_g+obs->K_v;
    std::cout << "----------------------TIME = " << vars->time/1000.0 << " ps-------------------------" << endl;
		cout<<"System propeties"<<endl;
    	printf("  Kion = %1.2e  Kgas = %1.2e  Kvap = %1.2e\n  Tion = %1.2f  Tgas = %1.2f  Tvap = %1.2f\n", 
		obs->Kion, obs->K_g, obs->K_v, obs->Tion, obs->T_g, obs->T_v);
		printf("Uion = %1.2e  Ugas = %1.2e  Uvap = %1.2e\n  Ugi = %1.2e  Ugg = %1.2e  Uvi = %1.2e	\n  Uvg = %1.2e	Uvv= %1.2e  \n",
		vars->Utotal.Uion, vars->Utotal.Ugas, vars->Utotal.Uvap, vars->Utotal.Ugi, vars->Utotal.Ugg,	vars->Utotal.Uvi, vars->Utotal.Uvg, vars->Utotal.Uvv);
		cout<<"System propeties"<<endl;
		printf("  K = %1.2e	U = %1.2e	Press = %f\n", K, U, gaspress/101300.0);
		cout<<"Times"<<endl;
		printf("  tion  = %1.1f s	tgas = %1.1f s 	tvap = %1.1f s\n",vars->times.tion,vars->times.tgas,vars->times.tvap);
		printf("  tvi   = %1.1f s	tgi  = %1.1f s	tvg  = %1.1f s	tvv  = %1.1f s\n",vars->times.tvi,vars->times.tgi,vars->times.tvg,vars->times.tvv);
		printf("  tpair = %1.1f s	tpos = %1.1f s	tvel = %1.1f s	tetc = %1.1f s\n",vars->times.tpair,vars->times.tpos,vars->times.tvel,vars->times.tetc);
		printf("  tpot  = %1.1f s	ttot = %1.1f s\n",(vars->times.tvi+vars->times.tgi+vars->times.tvv+vars->times.tvg+vars->times.tion+vars->times.tgas+vars->times.tvap), omp_get_wtime()-startTime);
		printf("  NCPU = %d\n",Nth);

		cout <<endl;

		sprintf(filepath, "K_%d.dat", int(calculation_number));
		FILE*f=fopen(filepath, "a");
		fprintf(f,"%e %e %e %e %e\n",vars->time,obs->Kion,obs->K_g,obs->K_v,K);
		fclose(f);

		sprintf(filepath, "U_%d.dat", int(calculation_number));
		f=fopen(filepath, "a");
		fprintf(f,"%e %e %e %e %e %e %e %e %e\n",vars->time,vars->Utotal.Uion,vars->Utotal.Ugas,vars->Utotal.Uvap,vars->Utotal.Ugi,vars->Utotal.Ugg,vars->Utotal.Uvg,vars->Utotal.Uvi,vars->Utotal.Uvv);
		fclose(f);
}

void
MD::export_dump(void) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int Ntotal=vars->mass.size();

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n",
		count, Ntotal, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {
		X=vars->position[0][0];
		Y=vars->position[0][1];
		Z=vars->position[0][2];
	}
	int ID=0;
	for(int i=0; i<Ntotal; i++){
		int Natom=vars->massAtom[i].size();
		std::array<double,3> &xi = vars->position[i];
		std::array<double,3> &vi = vars->velocity[i];
		int ty=vars->type[i];
		fprintf(f,"%d %s %f %f %f %f %f %f\n",
			ID,vars->atypes[ty].name.c_str(),xi[0]-X,xi[1]-Y,xi[2]-Z,vi[0],vi[1],vi[2]);
		ID++;
	}
	fclose(f);

}

void
MD::export_dump_close(void) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int Ntotal=0;
	int Nmol=vars->mass.size();
	for(int i=0; i<Nmol; i++){
		if(vars->region[i]){
			int Natom=vars->massAtom[i].size();
			Ntotal+=Natom;
		}
	}

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z vx vy vz\n",
		count, Ntotal, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {
		X=vars->position[0][0];
		Y=vars->position[0][1];
		Z=vars->position[0][2];
	}
	int ID=0;
	for(int i=0; i<Nmol; i++){
		if(vars->region[i]){
			int Natom=vars->massAtom[i].size();
			std::array<double,3> *xi = vars->positionAtom[i].data();
			std::array<double,3> *vi = vars->velocityAtom[i].data();
			for(int ii=0; ii<Natom; ii++){
				int ty=vars->typeAtom[i][ii];
				fprintf(f,"%d %s %f %f %f %f %f %f\n",
					ID,vars->atypes[ty].name.c_str(),xi[ii][0]-X,xi[ii][1]-Y,xi[ii][2]-Z,vi[ii][0],vi[ii][1],vi[ii][2]);
				ID++;
			}
		}
	}
	fclose(f);
}
