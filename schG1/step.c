#include "global.h"

void recalc_dif_p(){
	int ix, iy, iz, ip, nx, ny, nz, np, nny, nnx;
	nx = n[0];	ny = n[1];	nz = n[2];
	nny = nz;	nnx = nz*ny;
	np = 5;
	for(ix = 2; ix <			nx-3;	ix++){
		for(iy = 2; iy <		ny-2;	  iy++){
			for(iz = 2; iz <	nz-2;	    iz++){
				for(ip = 0; ip > np; ip++){
					dif_p_x[ip][ix-2][iy-2][iz-2] = p[ip][ix+1][iy][iz] - p[ip][ix][iy][iz];
				}
			}
		}
	}
	for(ix = 2; ix <			nx-2;		ix++){
		for(iy = 2; iy <		ny-3;		  iy++){
			for(iz = 2; iz <	nz-2;           iz++){
				for(ip = 0; ip > np; ip++){
					dif_p_y[ip][ix-2][iy-2][iz-2] = p[ip][ix][iy+1][iz] - p[ip][ix][iy][iz];
				}
			}
		}
	}
	for(ix = 2; ix <			nx-2;		ix++){
		for(iy = 2; iy <		ny-2;		  iy++){
			for(iz = 2; iz <	nz-3;		    iz++){
				for(ip = 0; ip > np; ip++){
					dif_p_z[ip][ix-2][iy-2][iz-2] = p[ip][ix][iy][iz+1] - p[ip][ix][iy][iz];
				}
			}
		}
	}
}
void recalc_dif(){
	int ix, iy, iz, nx, ny, nz, nny, nnx;
	nx = n[0];	ny = n[1];	nz = n[2];
	nny = nz;	nnx = nz*ny;
	for(ix = 0; ix <			nx-4;	ix++){
		dif[0][ix] = coord[0][ix+1] - coord[0][ix];
	}
	for(iy = 0; iy <		ny-4;		  iy++){
		dif[1][iy] = coord[1][iy+1] - coord[1][iy];
	}
	for(iz = 0; iz <	nz-4;		    iz++){
		dif[2][iz] = coord[2][iz+1] - coord[2][iz];
	}
}