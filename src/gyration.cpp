#include "md.hpp"

void
MD::gyration_out(MD *md2){
	FILE*f=fopen(pp->gyration_path, "a");
	double numerator=0;
	int Nion=vars->massAtom[0].size();
	std::array<double,3> *xi = vars->positionAtom[0].data();
	for(int i=0; i<Nion; i++){
		double dx=xi[i][0]-vars->position[0][0];
		double dy=xi[i][1]-vars->position[0][1];
		double dz=xi[i][2]-vars->position[0][2];
		double rsq=dx*dx+dy*dy+dz*dz;
		numerator+=vars->massAtom[0][i]*rsq;
	}
	fprintf(f, "%f\t%f\t",vars->time,sqrt(numerator/pp->Mion));
	numerator=0;
	std::array<double,3> *vi = vars->velocityAtom[0].data();
	for(int i=0; i<Nion; i++){
		double dx=vi[i][0]-vars->position[0][0];
		double dy=vi[i][1]-vars->position[0][1];
		double dz=vi[i][2]-vars->position[0][2];
		double rsq=dx*dx+dy*dy+dz*dz;
		numerator+=vars->massAtom[0][i]*rsq;
	}
	fprintf(f, "%f\n",sqrt(numerator/md2->pp->Mion));

	fclose(f);
}
