#include "../md.hpp"
/*########################################################################################

-----Periodic conditions-----

#######################################################################################*/


/************************periodec condition for in gas***************************/
void
MD::periodic(void) {
	std::array<double,3> *mol = vars->position.data();
	double HL=d_size*0.5;
	for (auto &a : vars->position) {
		if (a[0] < mol[0][0]-HL) a[0] += d_size;
		if (a[1] < mol[0][1]-HL) a[1] += d_size;
		if (a[2] < mol[0][2]-HL) a[2] += d_size;
		if (a[0] > mol[0][0]+HL) a[0] -= d_size;
		if (a[1] > mol[0][1]+HL) a[1] -= d_size;
		if (a[2] > mol[0][2]+HL) a[2] -= d_size;
	}
}

void
adjust_periodic(double &dx, double &dy, double &dz, double d_size) {
	const double LH = d_size * 0.5;
	if (dx < -LH)dx += d_size;
	if (dx > LH) dx -= d_size;
	if (dy < -LH)dy += d_size;
	if (dy > LH) dy -= d_size;
	if (dz < -LH)dz += d_size;
	if (dz > LH) dz -= d_size;
}