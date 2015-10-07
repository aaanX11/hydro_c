#define _CRTDBG_MAP_ALLOC
#include <string>
#include <sstream>
#include <stdlib.h>
#include <crtdbg.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "global.h"
#include "riemi.h"
#include "step.h"
void read_grid(const std::string& grd_fn){
	std::ifstream fp;
	std::string s;
	std::stringstream ss;
	fp.open(grd_fn.c_str());
	int ix, iy, iz, ip, nx, ny, nz, i;
	if(!fp){		std::cerr<<"Cannot open file with grid: "<<grd_fn<<"\n";	free_all();	exit(1);	}
	coord = new double*[3];
	dif = new double*[3];
	for(int ax = 0; ax < 3; ++ax){		
		while(getline(fp, s) != NULL && s[0] != 'X'+ax){	}
		getline(fp, s);	ss.str(s);	ss>>n[ax]; 	ss.clear();	//s.erase( remove(s.begin(), s.end(), '\r'), s.end() );
		if(n[ax] < 0){std::cerr<<"Grid file "<<grd_fn<<": strange number of points in "<<char(ax+'X')<<"\n"; fp.close();	free_all();		exit(1);		}
		coord[ax] = new double[n[ax]+1];
		dif[ax] = new double[ n[ax]];
		getline(fp, s); ss.clear();	ss.str(s); std::string s1;	i = 0; double tmp; while(ss>>tmp) {coord[ax][i] = tmp; i++;}ss.clear();
	}
	calc_dif();
	nx = n[0];	ny = n[1];	nz = n[2];
	if(nx*ny*nz > 1e+7){std::cerr<<"Grid file "<<grd_fn<<" has too many points in area: "<<n[0]*n[1]*n[2]<<"\n";free_all();	exit(1);}	
	p = new double***[5];
	if(p == NULL)std::cerr<<"NULL\n";
	for(int ip = 0; ip < 5; ip++){
		p[ip] = new double**[nx+4];
		if(p[ip] == NULL)std::cerr<<"NULL\n";
		for(ix = 0; ix < nx+4; ix++){
			p[ip][ix] = new double*[ny+4];
			if(p[ip][ix] == NULL)std::cerr<<"NULL\n";
			for(iy = 0; iy < ny+4; iy++){
				p[ip][ix][iy] = new double[nz+4];
				if(p[ip][ix][iy] == NULL){std::cerr<<"NULL\n";}
				for(iz = 0; iz < nz+4; iz++){
					p[ip][ix][iy][iz] = 0.0;
				}
			}
		}
	}
	lays = new int**[nx+4];
	if(lays == NULL){std::cerr<<"NULL\n";}
	for(ix = 0; ix < nx+4; ix++){
		lays[ix] = new int*[ny+4];
		if(lays[ix] == NULL){std::cerr<<"NULL\n";}
		for(iy = 0; iy < ny+4; iy++){
			lays[ix][iy] = new int[nz+4];
			if(lays[ix][iy] == NULL){std::cerr<<"NULL\n";}
			for(iz = 0; iz < nz+4; iz++){
				lays[ix][iy][iz] = -1;
			}
		}
	}
	dif_p_x = new double***[5];
	if(dif_p_x == NULL){std::cerr<<"NULL\n";}
	for(ip = 0; ip < 5; ip++){
		dif_p_x[ip] = new double**[nx+3];
		if(dif_p_x[ip] == NULL){std::cerr<<"NULL\n";}	
		for(ix = 0; ix < nx+3; ix++){
			dif_p_x[ip][ix] = new double*[ny+4];
			if(dif_p_x[ip][ix] == NULL){std::cerr<<"NULL\n";}
			for(iy = 0; iy < ny+4; iy++){
				dif_p_x[ip][ix][iy] = new double[nz+4];
				if(dif_p_x[ip][ix][iy] == NULL){std::cerr<<"NULL\n";}			
				for(iz = 0; iz < nz+4; iz++){
					dif_p_x[ip][ix][iy][iz] = 0.0;
				}
			}
		}
	}
	dif_p_y = new double***[5];
	for(ip = 0; ip < 5; ip++){
		dif_p_y[ip] = new double**[nx+4];
		for(ix = 0; ix < nx+4; ix++){
			dif_p_y[ip][ix] = new double*[ny+3];
			for(iy = 0; iy < ny+3; iy++){
				dif_p_y[ip][ix][iy] = new double[nz+4];
				for(iz = 0; iz < nz+4; iz++){
					dif_p_y[ip][ix][iy][iz] = 0.0;
				}
			}
		}
	}
	dif_p_z = new double***[5];
	for(ip = 0; ip < 5; ip++){
		dif_p_z[ip] = new double**[nx+4];
		for(ix = 0; ix < nx+4; ix++){
			dif_p_z[ip][ix] = new double*[ny+4];
			for(iy = 0; iy < ny+4; iy++){
				dif_p_z[ip][ix][iy] = new double[nz+3];
				for(iz = 0; iz < nz+3; iz++){
					dif_p_z[ip][ix][iy][iz] = 0.0;
				}
			}
		}
	}
	fp.close();
	return;
}
void free_all(){
	int ip, i_ax, ix, iy, nx, ny, nz;
	nx = n[0]; ny = n[1]; nz = n[2];
	for(i_ax = 0; i_ax < 3; i_ax++){
		delete[] coord[i_ax];
		delete[] dif[i_ax];	
	}	
	delete(coord);
	delete(dif);
	for(ip = 0; ip < 5; ip++){
		for(ix = 0; ix < nx+4; ix++){
			for(iy = 0; iy < ny+4; iy++){
				delete[] p[ip][ix][iy];
			}
			delete[] p[ip][ix];
		}
		delete[] p[ip];
	}
	delete[] p;
	
	for(ix = 0; ix < nx+4; ix++){
		for(iy = 0; iy < ny+4; iy++){
			delete[] lays[ix][iy];
		}
		delete[] lays[ix];
	}
	delete[] lays;	
	
	for(ip = 0; ip < 5; ip++){
		for(ix = 0; ix < nx+3; ix++){
			for(iy = 0; iy < ny+4; iy++){
				delete[] dif_p_x[ip][ix][iy];
			}
			delete[] dif_p_x[ip][ix];
		}
		delete[] dif_p_x[ip];	
	}
	delete[] dif_p_x;
	
	for(ip = 0; ip < 5; ip++){
		for(ix = 0; ix < nx+4; ix++){
			for(iy = 0; iy < ny+3; iy++){
				delete[] dif_p_y[ip][ix][iy];
			}
			delete[] dif_p_y[ip][ix];
		}
		delete[] dif_p_y[ip];	
	}
	delete[] dif_p_y;
	
	for(ip = 0; ip < 5; ip++){
		for(ix = 0; ix < nx+4; ix++){
			for(iy = 0; iy < ny+4; iy++){
				delete[] dif_p_z[ip][ix][iy];
			}
			delete[] dif_p_z[ip][ix];
		}
		delete[] dif_p_z[ip];	
	}
	delete[] dif_p_z;
	return;
}

void read_opt(const std::string opt_fn, int& stepmax, double& tmax, int& outp_size){
	std::string s;
	std::stringstream ss;
	std::ifstream opt(opt_fn.c_str());
	if(!opt){std::cerr<<"Cannot open file: "<<opt_fn<<"\n";	free_all();		exit(1);	}
	getline(opt, s);	ss.str(s);	ss>>stepmax;
	getline(opt, s);	ss.str(s);	ss>>tmax;
	getline(opt, s);	ss.str(s);	ss>>outp_size;
	opt.close();
}
bool recover(const std::string ini_f, int& istep, double& t){
	return false;
}
void read_val(const std::string& val_fn){
	std::ifstream fp;
	std::stringstream ss;
	std::string s;
	int nx, ny, nz;
	nx = n[0];	ny = n[1];	nz = n[2];
	fp.open(val_fn.c_str());
	if(!fp){std::cerr<<"Cannot open file with grid: "<<val_fn<<"\n";	free_all();		exit(1);	}
	for(int ix = 0; ix < nx+4; ix++){
		for(int iy = 0; iy < ny+4; iy++){
			getline(fp, s);	ss.clear();
			ss.str(s);
			int iz = 0;
			int tmp;
			while(ss>>tmp){lays[ix][iy][iz] = tmp;				iz++;			}
			if(iz != nz+4){std::cerr<<"Value file %s corrupted: layers"<<val_fn<<" "<<ix<<" "<<iy<<" "<<iz<<"\n";	
			free_all();		exit(1);	}
		}
	}
	for(int ip = 0; ip < 5; ip++){
		for(int ix = 0; ix < nx+4; ix++){
			for(int iy = 0; iy < ny+4; iy++){
				getline(fp, s);		ss.clear();	ss.str(s);
				int iz = 0;
				while(ss>>p[ip][ix][iy][iz]){				iz++;			}
				if(iz != nz+4){std::cerr<<"Value file "<<val_fn<<" "<<ip<<" "<<ix<<" "<<iy<<" "<<iz<<"\n";
				free_all();exit(1);}
			}
		}
	}
	fp.close();
	return;
}

void flux(double* par_ax, int ix, int iy, int iz, double* B, int ax, double& ddeltat){
	int alarm;
	double *par, rho1, u1, p1, ltan1, ltan2, rho2, u2, p2, rtan1, rtan2;
	par = new double[12];

	//forward
	riemparam2(par, ix, iy, iz, ax);	riemi(par, &alarm);
	if(alarm){
		std::cerr<<"function 'riemi' sends an alarm at ("<<ix<<" "<<iy<<" "<<iz<<'\n';}
	tanvel(par[10], rtan1, rtan2, ax, ix, iy, iz);
	updatetimestep(ax, ix, iy, iz, par[9], par[11], ddeltat);
	rho2 = par[6];		u2 = par[7];		p2 = par[8];

	//backward
	if(obstacle(ix, iy, iz, -ax) == 1){
		//tangential velosity components conserve near an obstacle
		riemparam2(par, ix, iy, iz, -ax);	riemi(par, &alarm);
		if(alarm){
			std::cerr<<"function 'riemi' sends an alarm at ("<<ix<<" "<<iy<<" "<<iz<<'\n';}
		updatetimestep(-ax, ix, iy, iz, par[9], par[11], ddeltat);
		par_ax[0] = par[6];		par_ax[1] = par[7];		par_ax[2] = par[8];		par_ax[3] = p[tan1(ax)][ix][iy][iz];		par_ax[4] = p[tan2(ax)][ix][iy][iz];
	}
	rho1 = par_ax[0];			u1 = par_ax[1];			p1 = par_ax[2];		ltan1 = par_ax[3];		ltan2 = par_ax[4];
	
	//fluxes
	B[0] = u2*rho2 - rho1*u1;
	switch(ax){
		case 'X':				
			B[1] = rho2*u2*u2 + p2 - (rho1*u1*u1 + p1);
			B[2] = rho2*u2*rtan1 -	(rho1*u1*ltan1);
			B[3] = rho2*u2*rtan2 -	(rho1*u1*ltan2);
			break;
		case 'Y':
			B[1] = rho2*u2*rtan1 - (rho1*u1*ltan1);
			B[2] = rho2*u2*u2 + p2 - (rho1*u1*u1 + p1);
			B[3] = rho2*u2*rtan2 - (rho1*u1*ltan2);
			break;
		case 'Z':
			B[1] = rho2*u2*rtan1 - rho1*u1*ltan1;
			B[2] = rho2*u2*rtan2 - rho1*u1*ltan2;
			B[3] = rho2*u2*u2 + p2 - (rho1*u1*u1 + p1);
			break;
	}
	B[4] = u2*(rho2*(inenergy(rho2,p2) + 0.5*(u2*u2+rtan1*rtan1+rtan2*rtan2)) + p2) - u1*(rho1*(inenergy(rho1,p1) + 0.5*(u1*u1+ltan1*ltan1+ltan2*ltan2)) + p1);
	
	delete[] par;
	
	//parameters' update
	par_ax[0] = rho2;		par_ax[1] = u2;		par_ax[2] = p2;		par_ax[3] = rtan1;		par_ax[4] = rtan2;
	return;
}
void step(double& t){
	int nx, ny, nz, alarm;
	double *A, 
		*Bx, *By, *Bz, 
		dt, 
		*parx, *pary, *parz;
	nx = n[0];	ny = n[1];	nz = n[2];
	parx = new double[(nz)*(ny)*5];
	pary = new double[(nz)*5];
	parz = new double[5];
	A = new double[5];
	Bx = new double[5];
	By = new double[5];
	Bz = new double[5];
	init_aux(parx, pary, parz);
	dt = 42;
	for(int ix = 2; ix < nx+2; ix++){
		for(int iy = 2; iy < ny+2; iy++){
			for(int iz = 2; iz < nz+2;iz++){
				if(obstacle(ix,iy,iz,0) > 0){
					flux(&parx[(iy-2)*nz + (iz-2)], ix, iy, iz, Bx, 'X', dt);
					flux(&pary[iz-2],			ix, iy, iz, By, 'Y', dt);
					flux(parz,				ix, iy, iz, Bz, 'Z', dt);
					f(ix, iy, iz, A);
					for(int i = 0; i < 5; i++){	A[i] = A[i] + dt*(Bx[i]/dif[0][ix-2] + By[i]/dif[1][iy-2] + Bz[i]/dif[2][iz-2]);	}
					f_inv(ix, iy, iz, A);
				}
			}
		}
	}
	recalc_dif_p();
	t = t + dt;
}
int main(){
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	int istep, nstep, outp_size;
	double t, tmax;

	std::string name("state");
	read_grid("KIRPICH.GRD");
	read_val("initval");
	read_opt("opt", nstep, tmax, outp_size);
	
	if(!recover("schG1.ini", istep, t)){istep = 0;	t = 0;}
	show(name);
	while(istep < nstep || t < tmax){
		step(t);
		istep++;
		std::stringstream ss; ss<<istep;
		std::string s(name);		s.append(ss.str());
		show(s);
	}

	free_all();
	return 0;
}