#include "../md.hpp"


void
MD::readCondFile(char* condfile){
	ifstream stream(condfile);
	string str;
	int iflag=0;
	cout<<"**************************************************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"************  Calculation conditions  ************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"**************************************************"<<endl;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		std::vector<string> readings;
		istringstream stream(str);
		string reading;
		while(getline(stream,reading,'\t')) {
			readings.push_back(reading);
		}
		if(readings[0]=="Input"){
			strcpy(atomFile,readings[1].c_str());
			cout<<"Atom file -->\t\t"<<readings[1]<<endl;
		}
		if(readings[0]=="Vapor"){
			strcpy(vaporFile,readings[1].c_str());
			Nof_around_vapor=stoi(readings[2]);
			cout<<"Vapor file -->\t\t"<<readings[1]<<endl;
			cout<<"Number of vapors\t"<<Nof_around_vapor<<endl;
		}
		if(readings[0]=="VaporStickPositions"){
			strcpy(vaporStickFile,readings[1].c_str());
			positionLogStep=stoi(readings[2]);
			ifstream stream2(vaporStickFile);
			string str2;
			while(getline(stream2,str2)) {
				if(str2.length()==0) continue;
				stickPositionList.push_back(stoi(str2)-1);
			}
			cout<<"Vapor stick position file -->\t\t"<<readings[1]<<endl;
		}
		if(readings[0]=="Gas"){
			std::array<double,3> x;
			std::array<double,3> newProp;
			newProp[0]=newProp[1]=newProp[2]=0;
			Atom_type at;
			if(readings[1]=="He"){
				gastype=1;
				at.mass=4.027;
				at.name="He(g)";
				at.coeff1=0.0203;
				at.coeff2=2.556;
				x[0]=0;
				x[1]=0;
				x[2]=0;
				vars->positionGas.push_back(x);
				vars->chargeGas.push_back(0);
				vars->typeGas.push_back(0);
			}
			if(readings[1]=="N2"){
				gastype=2;
				at.mass=14.01;
				at.name="N(g)";
				at.coeff1=0.1098;
				at.coeff2=3.27351824993;
				x[0]=-0.549;
				x[1]=0;
				x[2]=0;
				vars->positionGas.push_back(x);
				vars->chargeGas.push_back(0);
				vars->typeGas.push_back(0);
				vars->massGas.push_back(at.mass);
				x[0]=0.549;
				vars->positionGas.push_back(x);
				vars->chargeGas.push_back(0);
				vars->typeGas.push_back(0);
				vars->massGas.push_back(at.mass);
				Bond_type bt;
				bt.coeff[0]=1221.7;
				bt.coeff[1]=1.098;
				vars->btypes.push_back(bt);
				std::array<double,3> bon;
				bon[0]=0;
				bon[1]=1;
				bon[2]=0;
				//vars->bondsGasRigid.push_back(bon);
				vars->bondsGas.push_back(bon);
			}
			if(readings[1]=="N2monoatomic"){
				gastype=3;
				at.mass=28.02;
				at.name="N2(g)";
				at.coeff1=0.14397;
				at.coeff2=3.798;
				x[0]=0;
				x[1]=0;
				x[2]=0;
				vars->positionGas.push_back(x);
				vars->chargeGas.push_back(0);
			}
			if(readings[1]=="Ar"){
				gastype=4;
				at.mass=28.02;
				at.name="Ar(g)";
				at.coeff1=0.14397;
				at.coeff2=3.798;
				x[0]=0;
				x[1]=0;
				x[2]=0;
				vars->positionGas.push_back(x);
				vars->chargeGas.push_back(0);
			}
			vars->atypes.push_back(at);
			Nof_around_gas=stoi(readings[2]);
			cout<<"Gastype\t\t\t"<<readings[1]<<endl;
			cout<<"Number of gases\t\t"<<Nof_around_gas<<endl;
		}
		if(readings[0]=="Temperature"){
			T=stod(readings[1]);
			cout<<"Temperature\t\t"<<T<<" K"<<endl;
		}
		if(readings[0]=="Pressure"){
			p=stod(readings[1]);
			cout<<"Pressure\t\t"<<p<<" Pa"<<endl;
		}
		if (readings[0]=="dt") {
			dt=stod(readings[1]);
			cout<<"Time step\t\t"<<dt<<" fs"<<endl;
		}
		if (readings[0]=="TotalSteps") {
			Noftimestep=stod(readings[1]);
			cout<<"Total steps\t\t"<<float(Noftimestep)<<endl;
		}
		if (readings[0]=="RelaxSteps") {
			step_relax=stod(readings[1]);
			cout<<"Relax steps\t\t"<<float(step_relax)<<endl;
		}
		if (readings[0]=="CutOff") {
			CUTOFF=stod(readings[1]);
			cout<<"Cutoff\t\t\t"<<CUTOFF<<" ang."<<endl;
		}
		if (readings[0]=="Margin") {
			MARGIN=stod(readings[1]);
			cout<<"Margin size\t\t"<<MARGIN<<" ang."<<endl;
		}
		if (readings[0]=="Output") {
			ostringstream ss;
			ss<<readings[1]<<"_"<<calculation_number<<".dump";
			string tmp2=ss.str();
			pp->dump_path=new char[tmp2.length()+1];
			strcpy(pp->dump_path,tmp2.c_str());
			OBSERVE=stoi(readings[2]);
			FILE*f=fopen(pp->dump_path, "w");
			fclose(f);
			cout<<"Dump file -->\t\t"<<tmp2<<endl;
		}
		if (readings[0]=="NVTion") {
			if (readings[1]=="OFF") {
				flags->nose_hoover_ion=0;
				cout<<"Nose-Hoover for ion --> OFF"<<endl;
			}
			else if(readings[1]=="scale") {
				flags->nose_hoover_ion=0;
				flags->velocity_scaling=1;
				cout<<"Nose-Hoover for ion --> OFF\nVelocity scaling for ion --> ON"<<endl;
			}
			else {
				flags->nose_hoover_ion=1;
				pp->Tnh_ion=stod(readings[1]);
				cout<<"Nose-Hoover for ion --> ON --> "<<pp->Tnh_ion<<" K"<<endl;
			}
		}
		if (readings[0]=="NVTgas") {
			if (readings[1]=="OFF") {
				flags->nose_hoover_gas=0;
				cout<<"Nose-Hoover for gas --> OFF"<<endl;
			}
			else {
				flags->nose_hoover_gas=1;
				pp->Tnh_gas=stod(readings[1]);
				cout<<"Nose-Hoover for gas --> ON --> "<<pp->Tnh_gas<<" K"<<endl;
			}
		}

		if (readings[0]=="Interactions") {continue;}
		if (readings[1]=="gg") {
			if (readings[2]=="LJ") flags->inter_gg=1;
			else if (readings[2]=="OFF") flags->inter_gg=0;
			else printf("**************Uknown gas gas parameter was found**************\n");
		}
		if (readings[1]=="gi"||readings[1]=="ig") {
			if (readings[2]=="LJ") flags->force_lj=1;
			else if (readings[2]=="ion dipole") flags->force_ion_dipole=1;
			else if (readings[2]=="OFF") {
				flags->force_ion_dipole=0;
				flags->force_lj=0;
			}
			else printf("**************Uknown gas ion parameter was found**************\n");
		}
		if (readings[1]=="ion") {
			if (readings[2]=="AMBER") flags->intra_AMBER=1;
			else if (readings[2]=="Stilinger-Weber") flags->force_sw=1;
			else if (readings[2]=="Tersoff") flags->force_ters=1;
			else if (readings[2]=="Born-Mayer-Huggins-NaCl") flags->force_born=1;
			else printf("**************Uknown ion parameter was found**************\n");
		}
		if (readings[1]=="vi"||readings[1]=="iv") {
			if (readings[2]=="LJcoul") flags->inter_vi=1;
			else if (readings[2]=="OFF") flags->inter_vi=0;
			else printf("**************Uknown vapor ion parameter was found**************\n");
		}
		if (readings[1]=="vv"||readings[1]=="vv") {
			if (readings[2]=="LJcoul") flags->inter_vv=1;
			else if (readings[2]=="OFF") flags->inter_vv=0;
			else printf("**************Uknown vapor vapor parameter was found**************\n");
		}
		if (readings[1]=="gv"||readings[1]=="vg") {
			if (readings[2]=="LJ") flags->inter_vg=1;
			else if (readings[2]=="OFF") flags->inter_vg=0;
			else printf("**************Uknown vapor vapor parameter was found**************\n");
		}
		if (readings[1]=="Efield") {
			flags->efield=1;
			Ecoeff[0]=stod(readings[2]);
			Ecoeff[1]=stod(readings[3]);
			Ecoeff[2]=stod(readings[4]);
		}
		if (readings[0]=="Gyration") {
			ostringstream ss;
			ss<<readings[1]<<"_"<<calculation_number<<".dat";
			string tmp=ss.str();
			pp->gyration_path=new char[tmp.length()+10];
			strcpy(pp->gyration_path,tmp.c_str());
			flags->gyration=1;
			FILE*f=fopen(pp->gyration_path, "w");
			fclose(f);
			cout<<"Gyration --> ON -->\t"<<tmp<<endl;
		}
		if (readings[0]=="SHAKE") {
			ostringstream ss;
			if(readings[1]=="H") flags->shakeH=true;
			cout<<"SHAKE --> ON -->\tHydrogen"<<endl;
		}
	}
	vars->flags=flags;
	d_size=pow(Nof_around_gas*kb*T/p,1/3.0)*1e10;//pow(28.0855*8/6.02e23/(2.218e-24),1/3.0)*5;
	V=d_size*d_size*d_size;
	CL2 = (CUTOFF)*(CUTOFF);
	ML2 = (CUTOFF+MARGIN)*(CUTOFF+MARGIN);
	cout<<"Cut off length\t\t"<<CUTOFF<<" ang."<<endl;
	cout<<"Margin length\t\t"<<MARGIN<<" ang."<<endl;
	cout<<"Domain size\t\t"<<d_size<<" ang."<<endl;
}
