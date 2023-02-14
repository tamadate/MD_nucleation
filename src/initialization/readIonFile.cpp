#include "../variables.hpp"

void
Variables::readIonFile(char* infile){
	ifstream stream(infile);
	string str;
	int iflag=0;
	std::vector<std::array<double,3>> xs;
	std::vector<std::array<double,3>> vs;
	std::vector<std::array<double,3>> fs;
	std::vector<double> ms;
	std::vector<double> ps;
	std::vector<int> types;
	std::array<double,3> vec;
	std::vector<double> ws;
	vec[0]=0;
	vec[1]=0;
	vec[2]=0;

	std::vector<std::array<double,3>> bondsIon;
	std::vector<std::array<double,3>> bondsIonRigid;
	std::vector<std::array<double,4>> anglesIon;
	std::vector<std::array<double,5>> dihedralsIon;

	int atypesOffset=atypes.size()-1;
	int btypesOffset=btypes.size()-1;
	int ctypesOffset=ctypes.size()-1;
	int dtypesOffset=dtypes.size()-1;

	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="atom type name mass coeff1 coeff2") {iflag=1; continue;}
		if (str=="bond type name coeff1 coeff2") {iflag=2; continue;}
		if (str=="angle type name coeff1 coeff2") {iflag=3; continue;}
		if (str=="dihedral type name coeff1 coeff2 coeff3 coeff4") {iflag=4; continue;}
		if (str=="atoms") {iflag=5; continue;}
		if (str=="bonds") {iflag=6; continue;}
		if (str=="angles") {iflag=7; continue;}
		if (str=="dihedrals") {iflag=8; continue;}
	    string tmp;
    	istringstream stream(str);
		if (iflag==1) {
			Atom_type at;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==1) {
					if(tmp=="#c"||tmp=="#cs"||tmp=="#c1"||tmp=="#c2"||tmp=="#c3"||tmp=="#ca"
					||tmp=="#cp"||tmp=="#cq"||tmp=="#cc"||tmp=="#cd"||tmp=="#ce"||tmp=="#cf"
					||tmp=="#cg" ||tmp=="#ch"||tmp=="#cx"||tmp=="#cy"||tmp=="#cu"||tmp=="#cv"||tmp=="#cz") at.name="C";
					else if(tmp=="#h1"||tmp=="#h2"||tmp=="#h3"||tmp=="#h4"||tmp=="#h5"||tmp=="#ha"
					||tmp=="#hc"||tmp=="#hn"||tmp=="#ho"||tmp=="#hp"||tmp=="#hs"||tmp=="#hw"||tmp=="#hx") at.name="H";
					else if(tmp=="f") at.name="F";
					else if(tmp=="cl") at.name="Cl";
					else if(tmp=="br") at.name="Br";
					else if(tmp=="i") at.name="I";
					else if(tmp=="#n"||tmp=="#n1"||tmp=="#n2"||tmp=="#n3"||tmp=="#n4"||tmp=="#na"
					||tmp=="#nb"||tmp=="#nc"||tmp=="#nd"||tmp=="#ne"||tmp=="#nf"||tmp=="#nh"
					||tmp=="#no"||tmp=="#ns"||tmp=="#nt"||tmp=="#nx"||tmp=="#ny"||tmp=="#nz"
					||tmp=="#n+"||tmp=="#nu"||tmp=="#nv"||tmp=="#n7"||tmp=="#n8"||tmp=="#n9") at.name="N";
					else if(tmp=="#o"||tmp=="#oh"||tmp=="#os"||tmp=="#ow") at.name="O";
					else if(tmp=="#p2"||tmp=="#p3"||tmp=="#p4"||tmp=="#p5"||tmp=="#pb"||tmp=="#pc"
					||tmp=="#pd"||tmp=="#pe"||tmp=="#pf"||tmp=="#px"||tmp=="#py") at.name="P";
					else if(tmp=="#s"||tmp=="#s2"||tmp=="#s4"||tmp=="#s6"||tmp=="#sh"||tmp=="#ss"
					||tmp=="#sx"||tmp=="#sy") at.name="S";
					else at.name=tmp;
				}
				if (loop==2) at.mass=stod(tmp);
				if (loop==3) at.coeff1=stod(tmp);
				if (loop==4) at.coeff2=stod(tmp);
				loop++;
			}
			atypes.push_back(at);
		}
		if (iflag==2) {
			Bond_type bt;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) bt.coeff[0]=stod(tmp);
				if (loop==3) bt.coeff[1]=stod(tmp);
				loop++;
			}
			btypes.push_back(bt);
		}
		if (iflag==3) {
			Angle_type at;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) at.coeff[0]=stod(tmp);
				if (loop==3) at.coeff[1]=stod(tmp)/180.0*M_PI;
				loop++;
			}
			ctypes.push_back(at);
		}
		if (iflag==4) {
			Dihedral_type dit;
			int loop=0;
			double coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7, coeff8, coeff9, coeff10,  coeff11, coeff12, coeff13, coeff14, coeff15;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) dit.multi=stoi(tmp);
				if (loop==3) dit.coeff[0]=stod(tmp);
				if (loop==4) dit.coeff[1]=stod(tmp);
				if (loop==5) {
					dit.coeff[2]=stod(tmp);
					dit.coeff[3]=cos(stod(tmp)/180.0*M_PI);
					dit.coeff[4]=sin(stod(tmp)/180.0*M_PI);
				}
				if (loop==6) dit.coeff[5]=stod(tmp);
				if (loop==7) dit.coeff[6]=stod(tmp);
				if (loop==8) {
					dit.coeff[7]=stod(tmp);
					dit.coeff[8]=cos(stod(tmp)/180.0*M_PI);
					dit.coeff[9]=sin(stod(tmp)/180.0*M_PI);
				}
				if (loop==9) dit.coeff[10]=stod(tmp);
				if (loop==10) dit.coeff[11]=stod(tmp);
				if (loop==11) {
					dit.coeff[12]=stod(tmp);
					dit.coeff[13]=cos(stod(tmp)/180.0*M_PI);
					dit.coeff[14]=sin(stod(tmp)/180.0*M_PI);
				}
				loop++;
			}
			dtypes.push_back(dit);
		}
		if (iflag==5) {
			std::array<double,3> x;
			std::array<double,3> v;
			int loop=0;
			while(getline(stream,tmp,'\t')) {
				//if (loop==0) a.id=stoi(tmp);
				if (loop==1) {
					int ty=stoi(tmp)+atypesOffset;
					types.push_back(ty);
					ms.push_back(atypes[ty].mass);
				}
				if (loop==2) ps.push_back(stod(tmp));
				if (loop==3) x[0]=stod(tmp);
				if (loop==4) x[1]=stod(tmp);
				if (loop==5) x[2]=stod(tmp);
				if (loop==6) v[0]=stod(tmp);
				if (loop==7) v[1]=stod(tmp);
				if (loop==8) v[2]=stod(tmp);
				loop++;
			}
			xs.push_back(x);
			vs.push_back(v);
			fs.push_back(vec);
			ws.push_back(1);
		}
		if (iflag==6) {
			int loop=0;
			std::array<double,3> bondIon;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) bondIon[0]=stoi(tmp)-1;
				if (loop==1) bondIon[1]=stoi(tmp)-1;
				if (loop==2) bondIon[2]=stoi(tmp)+btypesOffset;
				loop++;
			}
			if(flags->shakeH){
				int type1=types[bondIon[0]];
				int type2=types[bondIon[1]];
				if(atypes[type1].name=="H"||atypes[type2].name=="H") bondsIonRigid.push_back(bondIon);
				else{bondsIon.push_back(bondIon);}
			}
			else bondsIon.push_back(bondIon);
		}
		if (iflag==7) {
			int loop=0;
			std::array<double,4> angleIon;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) angleIon[0]=stoi(tmp)-1;
				if (loop==1) angleIon[1]=stoi(tmp)-1;
				if (loop==2) angleIon[2]=stoi(tmp)-1;
				if (loop==3) angleIon[3]=stoi(tmp)+ctypesOffset;
				loop++;
			}
			anglesIon.push_back(angleIon);
		}
		if (iflag==8) {
			int loop=0;
			std::array<double,5> dihedralIon;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) dihedralIon[0]=stoi(tmp)-1;
				if (loop==1) dihedralIon[1]=stoi(tmp)-1;
				if (loop==2) dihedralIon[2]=stoi(tmp)-1;
				if (loop==3) dihedralIon[3]=stoi(tmp)-1;
				if (loop==4) dihedralIon[4]=stoi(tmp)+dtypesOffset;
				loop++;
			}
			dihedralsIon.push_back(dihedralIon);
		}
	}
	// Atoms push back
	positionAtom.push_back(xs);
	velocityAtom.push_back(vs);
	forceAtom.push_back(fs);
	massAtom.push_back(ms);
	chargeAtom.push_back(ps);
	weight.push_back(ws);
	typeAtom.push_back(types);
	bonds.push_back(bondsIon);
	angles.push_back(anglesIon);
	dihedrals.push_back(dihedralsIon);
	rigid_pairs.push_back(bondsIonRigid);
	
	// Molecule push back
	position.push_back(vec);
	velocity.push_back(vec);
	force.push_back(vec);
	type.push_back(0);
	double totalMass=0;
	for(auto &a : ms) totalMass+=a;
	mass.push_back(totalMass);
	region.push_back(true);
	group[0].push_back(0);
}
