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
		
		if (str=="atom type") {iflag=1; continue;}
		if (str=="bond type") {iflag=2; continue;}
		if (str=="angle type") {iflag=3; continue;}
		if (str=="dihedral type") {iflag=4; continue;}
		if (str=="atoms") {iflag=5; continue;}
		if (str=="bonds") {iflag=6; continue;}
		if (str=="angles") {iflag=7; continue;}
		if (str=="dihedrals") {iflag=8; continue;}
		std::vector<string> readings;
		istringstream stream(str);
		string reading;
		while(getline(stream,reading,'\t')) {
			readings.push_back(reading);
		}
		if (iflag==1) {
			Atom_type at;
			if(readings[1]=="#c"||readings[1]=="#cs"||readings[1]=="#c1"||readings[1]=="#c2"||readings[1]=="#c3"||readings[1]=="#ca"
			||readings[1]=="#cp"||readings[1]=="#cq"||readings[1]=="#cc"||readings[1]=="#cd"||readings[1]=="#ce"||readings[1]=="#cf"
			||readings[1]=="#cg" ||readings[1]=="#ch"||readings[1]=="#cx"||readings[1]=="#cy"||readings[1]=="#cu"||readings[1]=="#cv"||readings[1]=="#cz") at.name="C";
			else if(readings[1]=="#h1"||readings[1]=="#h2"||readings[1]=="#h3"||readings[1]=="#h4"||readings[1]=="#h5"||readings[1]=="#ha"
			||readings[1]=="#hc"||readings[1]=="#hn"||readings[1]=="#ho"||readings[1]=="#hp"||readings[1]=="#hs"||readings[1]=="#hw"||readings[1]=="#hx") at.name="H";
			else if(readings[1]=="f") at.name="F";
			else if(readings[1]=="cl") at.name="Cl";
			else if(readings[1]=="br") at.name="Br";
			else if(readings[1]=="i") at.name="I";
			else if(readings[1]=="#n"||readings[1]=="#n1"||readings[1]=="#n2"||readings[1]=="#n3"||readings[1]=="#n4"||readings[1]=="#na"
			||readings[1]=="#nb"||readings[1]=="#nc"||readings[1]=="#nd"||readings[1]=="#ne"||readings[1]=="#nf"||readings[1]=="#nh"
			||readings[1]=="#no"||readings[1]=="#ns"||readings[1]=="#nt"||readings[1]=="#nx"||readings[1]=="#ny"||readings[1]=="#nz"
			||readings[1]=="#n+"||readings[1]=="#nu"||readings[1]=="#nv"||readings[1]=="#n7"||readings[1]=="#n8"||readings[1]=="#n9") at.name="N";
			else if(readings[1]=="#o"||readings[1]=="#oh"||readings[1]=="#os"||readings[1]=="#ow") at.name="O";
			else if(readings[1]=="#p2"||readings[1]=="#p3"||readings[1]=="#p4"||readings[1]=="#p5"||readings[1]=="#pb"||readings[1]=="#pc"
			||readings[1]=="#pd"||readings[1]=="#pe"||readings[1]=="#pf"||readings[1]=="#px"||readings[1]=="#py") at.name="P";
			else if(readings[1]=="#s"||readings[1]=="#s2"||readings[1]=="#s4"||readings[1]=="#s6"||readings[1]=="#sh"||readings[1]=="#ss"
			||readings[1]=="#sx"||readings[1]=="#sy") at.name="S";
			else at.name=readings[1];
			at.mass=stod(readings[2]);
			at.coeff1=stod(readings[3]);
			at.coeff2=stod(readings[4]);
			atypes.push_back(at);
		}
		if (iflag==2) {
			Bond_type bt;
			bt.coeff[0]=stod(readings[2]);
			bt.coeff[1]=stod(readings[3]);
			btypes.push_back(bt);
		}
		if (iflag==3) {
			Angle_type at;
			at.coeff[0]=stod(readings[2]);
			at.coeff[1]=stod(readings[3])/180.0*M_PI;
			ctypes.push_back(at);
		}
		if (iflag==4) {
			Dihedral_type dit;
			dit.multi=stoi(readings[2]);
			dit.coeff[0]=stod(readings[3]);
			dit.coeff[1]=stod(readings[4]);
			dit.coeff[2]=stod(readings[5]);
			dit.coeff[3]=cos(stod(readings[5])/180.0*M_PI);
			dit.coeff[4]=sin(stod(readings[5])/180.0*M_PI);
			dit.coeff[5]=stod(readings[6]);
			dit.coeff[6]=stod(readings[7]);
			dit.coeff[7]=stod(readings[8]);
			dit.coeff[8]=cos(stod(readings[8])/180.0*M_PI);
			dit.coeff[9]=sin(stod(readings[8])/180.0*M_PI);
			dit.coeff[10]=stod(readings[9]);
			dit.coeff[11]=stod(readings[10]);
			dit.coeff[12]=stod(readings[11]);
			dit.coeff[13]=cos(stod(readings[11])/180.0*M_PI);
			dit.coeff[14]=sin(stod(readings[11])/180.0*M_PI);
			dtypes.push_back(dit);
		}
		if (iflag==5) {
			std::array<double,3> x;
			std::array<double,3> v;
			int ty=stoi(readings[1])+atypesOffset;
			types.push_back(ty);
			ms.push_back(atypes[ty].mass);
			ps.push_back(stod(readings[2]));
			x[0]=stod(readings[3]);
			x[1]=stod(readings[4]);
			x[2]=stod(readings[5]);
			if(readings.size()>6){
				v[0]=stod(readings[6]);
				v[1]=stod(readings[7]);
				v[2]=stod(readings[8]);
			}

			xs.push_back(x);
			vs.push_back(v);
			fs.push_back(vec);
			ws.push_back(1);
		}
		if (iflag==6) {
			std::array<double,3> bondIon;
			bondIon[0]=stoi(readings[0])-1;
			bondIon[1]=stoi(readings[1])-1;
			bondIon[2]=stoi(readings[2])+btypesOffset;
			if(flags->shakeH){
				int type1=types[bondIon[0]];
				int type2=types[bondIon[1]];
				if(atypes[type1].name=="H"||atypes[type2].name=="H") bondsIonRigid.push_back(bondIon);
				else{bondsIon.push_back(bondIon);}
			}
			else bondsIon.push_back(bondIon);
		}
		if (iflag==7) {
			std::array<double,4> angleIon;
			angleIon[0]=stoi(readings[0])-1;
			angleIon[1]=stoi(readings[1])-1;
			angleIon[2]=stoi(readings[2])-1;
			angleIon[3]=stoi(readings[3])+ctypesOffset;
			anglesIon.push_back(angleIon);
		}
		if (iflag==8) {
			std::array<double,5> dihedralIon;
			dihedralIon[0]=stoi(readings[0])-1;
			dihedralIon[1]=stoi(readings[1])-1;
			dihedralIon[2]=stoi(readings[2])-1;
			dihedralIon[3]=stoi(readings[3])-1;
			dihedralIon[4]=stoi(readings[4])+dtypesOffset;
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
