void recalc_dif_p(){
	int ix, iy, iz, nx, ny, nz;
	for(ix = 2; ix <			nx-3;	ix++){
		for(iy = 2; iy <		ny-2;	  iy++){
			for(iz = 2; iz <	nz-2;	    iz++){
				dif_p_x[ix-2][iy-2][iz-2] = p[ix+1][iy][iz] - p[ix][iy][iz]
			}
		}
	}
	for(ix = 2; ix <			nx-2;		ix++){
		for(iy = 2; iy <		ny-3;		  iy++){
			for(iz = 2; iz <	nz-2;           iz++){
				dif_p_y[ix-2][iy-2][iz-2] = p[ix][iy+1][iz] - p[ix][iy][iz]
			}
		}
	}
	for(ix = 2; ix <			nx-2;		ix++){
		for(iy = 2; iy <		ny-2;		  iy++){
			for(iz = 2; iz <	nz-3;		    iz++){
				dif_p_z[ix-2][iy-2][iz-2] = p[ix][iy][iz+1] - p[ix][iy][iz]
			}
		}
	}
}
void recalc_dif(){
	int ix, iy, iz, nx, ny, nz;
	for(ix = 0; ix <			nx-4;	ix++){
		for(iy = 0; iy <		ny-3;	  iy++){
			for(iz = 0; iz <	nz-3;	    iz++){
				dif_x[ix][iy][iz] = coord[ix+1][iy][iz] - coord[ix][iy][iz]
			}
		}
	}
	for(ix = 0; ix <			nx-3;		ix++){
		for(iy = 0; iy <		ny-4;		  iy++){
			for(iz = 0; iz <	nz-3;           iz++){
				dif_y[ix][iy][iz] = coord[ix][iy+1][iz] - coord[ix][iy][iz]
			}
		}
	}
	for(ix = 0; ix <			nx-3;		ix++){
		for(iy = 0; iy <		ny-3;		  iy++){
			for(iz = 0; iz <	nz-4;		    iz++){
				dif_z[ix][iy][iz] = coord[ix][iy][iz+1] - coord[ix][iy][iz]
			}
		}
	}
}