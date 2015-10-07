#pragma once
#include "global.h"

void recalc_dif_p();
void free_all();
void calc_dif();
int tan1(int ax);
int tan2(int ax);
void updatetimestep(int ax, int ix, int iy, int iz, double s1, double s3, double& ddeltat);
void init_aux(double* parx, double* pary, double* parz);
void tanvel(double s2, double& rtan1, double& rtan2, int ax, int ix, int iy, int iz);
int obstacle(int ix, int iy, int iz, int axis);
double inenergy(double rho, double pressure);
void riemparam1(double* param, int ix, int iy, int iz, int axis);
void riemparam2(double* param, int ix, int iy, int iz, int axis);
void f(int ix, int iy, int iz, double* A);
void f_inv(int ix, int iy, int iz, double* A);
void show(std::string);