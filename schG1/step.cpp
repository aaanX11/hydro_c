#pragma once
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include "global.h"
#include "step.h"
#include "riemi.h"

void recalc_dif_p(){
	int ix, iy, iz, ip, nx, ny, nz, np;
	double max, min;
	nx = n[0];	ny = n[1];	nz = n[2];
	fprintf(stderr, "recalc_dif_p: nx = %d ny = %d nz = %d\n", nx, ny, nz);
	np = 5;
	for(ip = 0; ip < np; ip++){
		max = -1000;
		min = 1000;
		for(ix = 2; ix <			nx+1;	ix++){
			for(iy = 2; iy <		ny+2;	  iy++){
				for(iz = 2; iz <	nz+2;	    iz++){
					dif_p_x[ip][ix-2][iy-2][iz-2] = p[ip][ix+1][iy][iz] - p[ip][ix][iy][iz];
					if(min > dif_p_x[ip][ix-2][iy-2][iz-2]){min = dif_p_x[ip][ix-2][iy-2][iz-2];}
					if(max < dif_p_x[ip][ix-2][iy-2][iz-2]){max = dif_p_x[ip][ix-2][iy-2][iz-2];}
				}
			}
		}
		fprintf(stderr, "max dif p = %lf min dif p = %lf\n", max, min);
	}
	for(ip = 0; ip < np; ip++){
		max = -1000;
		min = 1000;
		for(ix = 2; ix <			nx+2;		ix++){
			for(iy = 2; iy <		ny+1;		  iy++){
				for(iz = 2; iz <	nz+2;           iz++){
					dif_p_y[ip][ix-2][iy-2][iz-2] = p[ip][ix][iy+1][iz] - p[ip][ix][iy][iz];
					if(min > dif_p_y[ip][ix-2][iy-2][iz-2]){min = dif_p_y[ip][ix-2][iy-2][iz-2];}
					if(max < dif_p_y[ip][ix-2][iy-2][iz-2]){max = dif_p_y[ip][ix-2][iy-2][iz-2];}
				}
			}
		}
		fprintf(stderr, "max dif p = %lf min dif p = %lf\n", max, min);
	}
	for(ip = 0; ip < np; ip++){
		max = -1000;
		min =1000;
		for(ix = 2; ix <			nx+2;		ix++){
			for(iy = 2; iy <		ny+2;		  iy++){
				for(iz = 2; iz <	nz+1;		    iz++){
					dif_p_z[ip][ix-2][iy-2][iz-2] = p[ip][ix][iy][iz+1] - p[ip][ix][iy][iz];
					if(min > dif_p_z[ip][ix-2][iy-2][iz-2]){min = dif_p_z[ip][ix-2][iy-2][iz-2];}
					if(max < dif_p_z[ip][ix-2][iy-2][iz-2]){max = dif_p_z[ip][ix-2][iy-2][iz-2];}
				}
			}
		}
		fprintf(stderr, "max dif p = %lf min dif p = %lf\n", max, min);
	}
}
void calc_dif(){
	int ix, iy, iz, nx, ny, nz;
	double max, min;
	nx = n[0];	ny = n[1];	nz = n[2];
	fprintf(stderr, "%d %d %d\n", nx, ny, nz);
	for(ix = 0; ix <			nx;	ix++){
		dif[0][ix] = coord[0][ix+1] - coord[0][ix];
	}
	for(iy = 0; iy <		ny;		  iy++){
		dif[1][iy] = coord[1][iy+1] - coord[1][iy];
	}
	min = 20;
	max = -20;
	for(iz = 0; iz <	nz;		    iz++){
		dif[2][iz] = coord[2][iz+1] - coord[2][iz];
		if(dif[2][iz] < min){
			min = dif[2][iz];
		}
		if(dif[2][iz] > max){
			max = dif[2][iz];
		}
	}
	fprintf(stderr, "%lf %lf\n", min, max);
}
int tan1(int ax){
	switch(ax){
		case 'X':
			return 1;
		case 'Y':
			return 0;
		case 'Z':
			return 0;
	}
	std::cerr<<"\n";
	exit(1);
	return 4;
}
int tan2(int ax){
	switch(ax){
		case 'X':
			return 2;
		case 'Y':
			return 2;
		case 'Z':
			return 1;
	}
	std::cerr<<"\n";
	exit(1);
	return 4;
}
void tanvel(double s2, double& rtan1, double& rtan2, int ax, int ix, int iy, int iz){
	//forward
	if (s2 > 0 || obstacle(ix, iy, iz, ax)){
		rtan1 = p[tan1(ax)][ix][iy][iz];
		rtan2 = p[tan2(ax)][ix][iy][iz];
	}
	else{
		switch(ax){
			case 'X':
				rtan1 = p[tan1(ax)][ix+1][iy][iz];
				rtan2 = p[tan2(ax)][ix+1][iy][iz];
				break;
			case 'Y':
				rtan1 = p[tan1(ax)][ix][iy+1][iz];
				rtan2 = p[tan2(ax)][ix][iy+1][iz];
				break;
			case 'Z':
				rtan1 = p[tan1(ax)][ix][iy][iz+1];
				rtan2 = p[tan2(ax)][ix][iy][iz+1];
				break;
		}
	}
}
void init_aux(double* parx, double* pary, double* parz){
	//parx:		5*ny*nz
	//pary:		5*nz
	//parz:		5
	double* param, rtan1, rtan2;
	int ip, ax, nx, ny, nz, alarm;
	nx = n[0]; ny = n[1]; nz = n[2];
	param = new double[12];
	ip = 1;
	ax = 'X';
	for(int iy = 2; iy < ny+2; iy++){
		for(int iz = 2; iz < nz+2; iz++){
			riemparam1(param, 2, iy, iz, -ax);
			riemi(param, &alarm);
			if(alarm){
				std::cerr<<"function 'riemi' sends an alarm"<<'\n';}
			tanvel(param[10], rtan1, rtan2, ax, 2, iy, iz);
			parx[0+(iy-2)*nz+(iz-2)] = param[0];	parx[1+(iy-2)*nz+(iz-2)] = param[1];	parx[2+(iy-2)*nz+(iz-2)] = param[2];	parx[3+(iy-2)*nz+(iz-2)] = rtan1;	parx[4+(iy-2)*nz+(iz-2)] = rtan2;
		}
	}
	ip = 2;
	ax = 'Y';
	for(int iz = 2; iz < nz+2; iz++){
		riemparam1(param, 2, 2, iz, -ax);
		riemi(param, &alarm);
		if(alarm){
			std::cerr<<"function 'riemi' sends an alarm"<<'\n';}
		tanvel(param[10], rtan1, rtan2, ax, 2, 2, iz);
		pary[0+iz-2] = param[0];	pary[1+iz-2] = param[1];	pary[2+iz-2] = param[2];	pary[3+iz-2] = rtan1;	pary[4+iz-2] = rtan2;
	}
	ip = 3;
	int iz = 2;
	ax = 'Z';
	riemparam1(param, 2, 2, 2, -ax);	
	riemi(param, &alarm);
	if(alarm){
		std::cerr<<"function 'riemi' sends an alarm"<<'\n';}
	tanvel(param[10], rtan1, rtan2, ax, 2, 2, 2);
	parz[0] = param[0];	parz[1] = param[1];	parz[2] = param[2];	parz[3] = rtan1;	parz[4] = rtan2;
}
void updatetimestep(int ax, int ix, int iy, int iz, double s1, double s3, double& ddeltat){
	int i;
	switch(ax){
		case 'X':
			i = ix-2;
			break;
		case 'Y':
			i = iy-2;
			break;
		case 'Z':
			i = iz-2;
			break;
		case -'X':
			i = ix-2-1;
			break;
		case -'Y':
			i = iy-2-1;
			break;
		case -'Z':
			i = iz-2-1;
			break;
	}
	if(fabs(s1) < 1e-5 || fabs(s1) > 1e+2){
		std::cerr<<"function 'updatetimestep' s1 too big "<<char(ax)<<" "<<i<<"\n";	}
	if(fabs(s3) < 1e-5 || fabs(s3) > 1e+2){
		std::cerr<<"function 'updatetimestep' s3 too big "<<char(ax)<<" "<<i<<"\n";	}
	if(i < -1 || i > n[ax-'X']-1){				
		std::cerr<<"function 'updatetimestep' index out of range "<<char(ax)<<" "<<i<<"\n";	}
	if(i >= 0){				ddeltat = std::min(ddeltat, dif[ax-'X'][i]/fabs(s1));			}
	if(i <= n[ax-'X']-2){				ddeltat = std::min(ddeltat, dif[ax-'X'][i+1]/fabs(s3));			}

	if(ddeltat < 1e-5 || ddeltat > 2){
		std::cerr<<char(ax)<<"function 'updatetimestep' time step error dt = "<<ddeltat<<'\n';}
	return;	
}
int obstacle(int ix, int iy, int iz, int axis){
	if(axis == 0){
		if(lays[ix][iy][iz] == -1){return 0;}
		else return 3;
	}
	switch(axis){
		case 'X':
			if (lays[ix+2][iy][iz] == -1){return 2;}
			else{if(lays[ix+1][iy][iz] == -1) return 1;
			else{return 3;}
			}
		case 'Y':
			if (lays[ix][iy+1][iz] == -1){return 1;}
			else{if(lays[ix][iy+2][iz] == -1){ return 2;}
			else return 3;}
		case 'Z':
			if (lays[ix][iy][iz+1] == -1){return 1;}
			else{if(lays[ix][iy][iz+2] == -1){ return 2;}
			else return 3;}
		case -('X'):
			if (lays[ix-1][iy][iz] == -1){return 1;}
			else{if( lays[ix-2][iy][iz] == -1){ return 2;}
			else return 3;}
		case -('Y'):
			if (lays[ix][iy-1][iz] == -1){return 1;}
			else{if(lays[ix][iy-2][iz] == -1){ return 2;}
			else return 3;}
		case -('Z'):
			if (lays[ix][iy][iz-1] == -1){return 1;}
			else{if(lays[ix][iy][iz-2] == -1){ return 2;}
			else return 3;}
	}
	return -1;
}
double mmod(double x, double y){
	double t;
	if(x*y < 0){	return 0;	}
	if (fabs(x) < fabs(y)){t = fabs(x);}
	else{ t = fabs(y);}
	if(x < 0){	t *= -1; 	}
	return t;
}
double inenergy(double rho, double pressure){
	double kappa = 1.4;
	return pressure/((kappa - 1)*rho);
}
double energyconverter(double rho, double u, double v, double w, double e){
	double kappa = 1.4;	
	return (kappa - 1)*(e - rho*0.5*(u*u + v*v +w*w));
}
void riemparam1(double* param, int ix, int iy, int iz, int ax){
	int ip;
	switch(abs(ax)){
		case 'X':
			ip = 1;
			break;
		case 'Y':
			ip = 2;
			break;
		case 'Z':
			ip = 3;
			break;
	}
	if(obstacle(ix,iy,iz,ax)  == 1){
		param[0] = p[0][ix][iy][iz];	param[3] = -p[0][ix][iy][iz];
		param[1] = p[ip][ix][iy][iz];	param[4] = -p[ip][ix][iy][iz];
		param[2] = p[4][ix][iy][iz];	param[5] = -p[4][ix][iy][iz];
	}
	else{
		param[0] = param[3] = p[0][ix][iy][iz];
		param[1] = param[4] = p[ip][ix][iy][iz];
		param[2] = param[5] = p[4][ix][iy][iz];
	}
}
void riemparam2(double* param, int ix, int iy, int iz, int axis){
	if(obstacle(ix,iy,iz,axis) < 3){
		riemparam1(param,ix,iy,iz,axis);
	}
	else{
	int ip;
	switch(axis){
		case 'X':
			if(ix <  3 || ix > n[0]){riemparam1(param,ix,iy,iz,axis);}
			else{
				ip = 1;
				param[0] = p[0][ix][iy][iz]   + 0.5*dif[0][ix-2]*mmod(dif_p_x[0][ix][iy][iz]/dif[0][ix-2],		dif_p_x[0][ix-1][iy][iz]/dif[0][ix-3]);
				param[3] = p[0][ix+1][iy][iz] - 0.5*dif[0][ix-1]*mmod(dif_p_x[0][ix+1][iy][iz]/dif[0][ix-1], dif_p_x[0][ix][iy][iz]/dif[0][ix-2]);
				param[1] = p[ip][ix][iy][iz]  + 0.5*dif[0][ix-2]*mmod(dif_p_x[ip][ix][iy][iz]/dif[0][ix-2],	dif_p_x[ip][ix-1][iy][iz]/dif[0][ix-3]);
				param[4] = p[ip][ix+1][iy][iz]- 0.5*dif[0][ix-1]*mmod(dif_p_x[ip][ix+1][iy][iz]/dif[0][ix-1],dif_p_x[ip][ix][iy][iz]/dif[0][ix-2]);
				param[2] = p[4][ix][iy][iz]   + 0.5*dif[0][ix-2]*mmod(dif_p_x[4][ix][iy][iz]/dif[0][ix-2],		dif_p_x[4][ix-1][iy][iz]/dif[0][ix-3]);
				param[5] = p[4][ix+1][iy][iz] - 0.5*dif[0][ix-1]*mmod(dif_p_x[4][ix+1][iy][iz]/dif[0][ix-1], dif_p_x[4][ix][iy][iz]/dif[0][ix-2]);
			}
			break;
		case 'Y':
			if(iy <  3 || iy > n[1]){riemparam1(param,ix,iy,iz,axis);}
			else{
				ip = 2;
				param[0] = p[0][ix][iy][iz]   + 0.5*dif[1][iy-2]*mmod(dif_p_y[0][ix][iy][iz]/dif[1][iy-2],		dif_p_y[0][ix][iy-1][iz]/dif[1][iy-3]);
				param[3] = p[0][ix][iy+1][iz] - 0.5*dif[1][iy-1]*mmod(dif_p_y[0][ix][iy+1][iz]/dif[1][iy-1], dif_p_y[0][ix][iy][iz]/dif[1][iy-2]);
				param[1] = p[ip][ix][iy][iz]  + 0.5*dif[1][iy-2]*mmod(dif_p_y[ip][ix][iy][iz]/dif[1][iy-2],	dif_p_y[ip][ix][iy-1][iz]/dif[1][iy-3]);
				param[4] = p[ip][ix][iy+1][iz]- 0.5*dif[1][iy-1]*mmod(dif_p_y[ip][ix][iy+1][iz]/dif[1][iy-1],dif_p_y[ip][ix][iy][iz]/dif[1][iy-2]);
				param[2] = p[4][ix][iy][iz]   + 0.5*dif[1][iy-2]*mmod(dif_p_y[4][ix][iy][iz]/dif[1][iy-2],		dif_p_y[4][ix][iy-1][iz]/dif[1][iy-3]);
				param[5] = p[4][ix][iy+1][iz] - 0.5*dif[1][iy-1]*mmod(dif_p_y[4][ix][iy+1][iz]/dif[1][iy-1], dif_p_y[4][ix][iy][iz]/dif[1][iy-1]);
			}
			break;
		case 'Z':
			if(iz <  3 || iz > n[2]){riemparam1(param,ix,iy,iz,axis);}
			else{
				ip = 3;
				param[0] = p[0][ix][iy][iz]   + 0.5*dif[2][iz-2]*mmod(dif_p_x[0][ix][iy][iz]/dif[2][iz-2],		dif_p_z[0][ix][iy][iz-1]/dif[2][iz-3]);
				param[3] = p[0][ix][iy][iz+1] - 0.5*dif[2][iz-1]*mmod(dif_p_z[0][ix][iy][iz+1]/dif[2][iz-1], dif_p_z[0][ix][iy][iz]/dif[2][iz-2]);
				param[1] = p[ip][ix][iy][iz]  + 0.5*dif[2][iz-2]*mmod(dif_p_z[ip][ix][iy][iz]/dif[2][iz-1],	dif_p_z[ip][ix][iy][iz-1]/dif[2][iz-3]);
				param[4] = p[ip][ix][iy][iz+1]- 0.5*dif[2][iz-1]*mmod(dif_p_z[ip][ix][iy][iz+1]/dif[2][iz-1],dif_p_z[ip][ix][iy][iz]/dif[2][iz-2]);
				param[2] = p[4][ix][iy][iz]   + 0.5*dif[2][iz-2]*mmod(dif_p_z[4][ix][iy][iz]/dif[2][iz-2],		dif_p_z[4][ix][iy][iz-1]/dif[2][iz-3]);
				param[5] = p[4][ix][iy][iz+1] - 0.5*dif[2][iz-1]*mmod(dif_p_z[4][ix][iy][iz+1]/dif[2][iz-1], dif_p_z[4][ix][iy][iz]/dif[2][iz-2]);
			}
			break;
	}
	}
}
void f(int ix, int iy, int iz, double* A){
	double kappa = 1.4;	
	A[0] = p[0][ix][iy][iz];
	A[1] = p[0][ix][iy][iz]*p[1][ix][iy][iz];
	A[2] = p[0][ix][iy][iz]*p[2][ix][iy][iz];
	A[3] = p[0][ix][iy][iz]*p[3][ix][iy][iz];
	double ie = p[4][ix][iy][iz]/((kappa - 1)*p[0][ix][iy][iz]);
	A[4] = p[0][ix][iy][iz]*(ie + 0.5*(p[1][ix][iy][iz]*p[1][ix][iy][iz] + p[2][ix][iy][iz]*p[2][ix][iy][iz] + p[3][ix][iy][iz]*p[3][ix][iy][iz]));
}
void f_inv(int ix, int iy, int iz, double* A){
	double kappa = 1.4;	
	if(fabs(A[0])>1e+2||fabs(A[1])>1e+2||fabs(A[2])>1e+2||fabs(A[3])>1e+2||fabs(A[4])>1e+2){
		std::cerr<<"Collapse\n";
	}
	p[0][ix][iy][iz] = A[0];
	p[1][ix][iy][iz] = A[1]/A[0];
	p[2][ix][iy][iz] = A[2]/A[0];
	p[3][ix][iy][iz] = A[3]/A[0];
	p[4][ix][iy][iz] = (kappa - 1)*(A[4] - 0.5*(A[1]*A[1] + A[1]*A[1] + A[1]*A[1])/A[0]);
}

void show(std::string datafname){
	std::ofstream datafile((datafname.append(".dat")).c_str());
	int nx, ny, nz, i[3];
	nx = n[0];	ny = n[1];	nz = n[2];
	datafile<<"VARIABLES = \"X\", \"Y\", \"Z\", \"R\", \"U\", \"V\",\"W\", \"P\", \"L\"\n";
	datafile<<"ZONE F=BLOCK I= "<<nx+1<<" J= "<<ny+1<<" K = "<<nz+1<<'\n';
	datafile<<"VARLOCATION=([4-9]=CELLCENTERED)\n";
	for(int ax = 0; ax < 3; ax++){
		for( i[2] = 0; i[2] <= nz; ++i[2]){	
			for( i[1] = 0; i[1] <= nz; ++i[1]){	
				for( i[0] = 0; i[0] <= nz; ++i[0]){	datafile<<coord[ax][i[ax]]<<'\t';	}
				datafile<<'\n';
			}
		}
	}
	for(int ip = 0; ip < 5; ip++){
		for(int iz = 2; iz < nz+2; iz++){
			for(int iy = 2; iy < ny+2; iy++){
				for(int ix = 2; ix < nx+2; ix++){	datafile<<p[ip][ix][iy][iz]<<'\t';}
				datafile<<'\n';
			}
		}	
	}
	for(int iz = 2; iz < nz+2; iz++){
		for(int iy = 2; iy < ny+2; iy++){
			for(int ix = 2; ix < nx+2; ix++){	datafile<<lays[ix][iy][iz]<<'\t';}
			datafile<<'\n';
		}
	}
	datafile.close();
}