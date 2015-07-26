#include <stdio.h>
#include <string.h>

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include "global.h"

void read_grid(char* grd_fn){
	FILE* fp;
	char buff[2048];
	char *buf;
	char tmp[32];
	char ax[2];
	int i_ax, i, nbytes;
	int ix, iy, iz, ip, nx, ny, nz;
	buf = (char*)malloc(sizeof(char)*2048);
	if((fp = fopen(grd_fn, "r")) == NULL){
		fprintf(stderr, "Cannot open file with grid: %s\n", grd_fn);
		exit(1);
	}
	PARMS_malloc(coord, 3, double*);
	PARMS_malloc(dif, 3, double*);
	for(i_ax = 0; i_ax < 3; i_ax++){
		ax[0] = 'X' + i_ax;
		ax[1] = '\0';
		while(fgets(buff, 256, fp) != NULL && sscanf(buff,"%s", tmp) && strcmp(tmp,ax)){}
			
		fgets(buff, 256, fp);
		sscanf(buff, "%d", &n[i_ax]);
		if(n[i_ax] < 0){
			fprintf(stderr, "Grid file %s: strange number of points in %c\n", grd_fn, 'X'+ i_ax);
			exit(1);
		}
		PARMS_malloc(coord[i_ax], n[i_ax]+1, double);
		PARMS_malloc(dif[i_ax], n[i_ax], double);
		i = 0;
		fgets(buff, 2048, fp);
		buf = buff;
		while(sscanf(buf, "%lf%n", &coord[i_ax][i], &nbytes) == 1){
			buf += nbytes;
			i++;
		}
		if(i != n[i_ax]+1){
			fprintf(stderr, "Grid file %s corrupted: points in %c less then expected\n", grd_fn, 'X'+i_ax);
			exit(1);
		}
	}
	nx = n[0];
	ny = n[1];
	nz = n[2];
	if(nx*ny*nz > 1e+7){
		fprintf(stderr, "Grid file %s too many points in area: \n", grd_fn, n[0]*n[1]*n[2]);
		exit(1);
	}
	
	PARMS_malloc(p, 5, double***);
	if(p == NULL){fprintf(stderr, "NULL\n");}
	for(ip = 0; ip < 5; ip++){
		PARMS_malloc(p[ip], nx+4, double**);
		if(p[ip] == NULL){fprintf(stderr, "NULL\n");}
		for(ix = 0; ix < nx+4; ix++){
			PARMS_malloc(p[ip][ix], ny+4, double*);
			if(p[ip][ix] == NULL){fprintf(stderr, "NULL\n");}
			for(iy = 0; iy < ny+4; iy++){
				PARMS_malloc(p[ip][ix][iy], nz+4, double);
				if(p[ip][ix][iy] == NULL){fprintf(stderr, "NULL\n");}
				for(iz = 0; iz < nz+4; iz++){
					p[ip][ix][iy][iz] = 0.0;
				}
			}
		}
	}
	
	PARMS_malloc(lays, nx+4, int**);	
	if(lays == NULL){fprintf(stderr, "NULL\n");}
	for(ix = 0; ix < nx+4; ix++){
		PARMS_malloc(lays[ix], ny+4, int*);
		if(lays[ix] == NULL){fprintf(stderr, "NULL\n");}
		for(iy = 0; iy < ny+4; iy++){
			PARMS_malloc(lays[ix][iy], nz+4, int);
			if(lays[ix][iy] == NULL){fprintf(stderr, "NULL\n");}
			for(iz = 0; iz < nz+4; iz++){
				lays[ix][iy][iz] = 0;
			}
		}
	}
	
	PARMS_malloc(dif_p_x, 5, double***);
	if(dif_p_x == NULL){fprintf(stderr, "NULL\n");}
	for(ip = 0; ip < 5; ip++){
		PARMS_malloc(dif_p_x[ip], nx+4, double**);	
		if(dif_p_x[ip] == NULL){fprintf(stderr, "NULL\n");}	
		for(ix = 0; ix < nx+4; ix++){
			PARMS_malloc(dif_p_x[ip][ix], ny+4, double*);
			if(dif_p_x[ip][ix] == NULL){fprintf(stderr, "NULL\n");}
			for(iy = 0; iy < ny+4; iy++){
				PARMS_malloc(dif_p_x[ip][ix][iy], nz+4, double);
				if(dif_p_x[ip][ix][iy] == NULL){fprintf(stderr, "NULL\n");}			
				for(iz = 0; iz < nz+4; iz++){
					dif_p_x[ip][ix][iy][iz] = 0.0;
				}
			}
		}
	}
	PARMS_malloc(dif_p_y, 5, double***);
	for(ip = 0; ip < 5; ip++){
		PARMS_malloc(dif_p_y[ip], nx+4, double**);	
		for(ix = 0; ix < nx+4; ix++){
			PARMS_malloc(dif_p_y[ip][ix], ny+4, double*);
			for(iy = 0; iy < ny+4; iy++){
				PARMS_malloc(dif_p_y[ip][ix][iy], nz+4, double);
				for(iz = 0; iz < nz+4; iz++){
					dif_p_y[ip][ix][iy][iz] = 0.0;
				}
			}
		}
	}
	PARMS_malloc(dif_p_z, 5, double***);
	for(ip = 0; ip < 5; ip++){
		PARMS_malloc(dif_p_z[ip], nx+4, double**);	
		for(ix = 0; ix < nx+4; ix++){
			PARMS_malloc(dif_p_z[ip][ix], ny+4, double*);
			for(iy = 0; iy < ny+4; iy++){
				PARMS_malloc(dif_p_z[ip][ix][iy], nz+4, double);
				for(iz = 0; iz < nz+4; iz++){
					dif_p_z[ip][ix][iy][iz] = 0.0;
				}
			}
		}
	}
	fclose(fp);
	//free(buf);
	return;
}
void free_all(){
	int ip, i_ax, ix, iy, iz, nx, ny, nz;
	nx = n[0]; ny = n[1]; nz = n[2];
	for(i_ax = 0; i_ax < 3; i_ax++){
		free(coord[i_ax]);
		free(dif[i_ax]);
	}	
	free(coord);
	free(dif);
	for(ip = 0; ip < 5; ip++){
		for(ix = 0; ix < nx+4; ix++){
			for(iy = 0; iy < ny+4; iy++){
				free(p[ip][ix][iy]);
			}
			free(p[ip][ix]);
		}
		free(p[ip]);
	}
	free(p);
	
	for(ix = 0; ix < nx+4; ix++){
		for(iy = 0; iy < ny+4; iy++){
			free(lays[ix][iy]);
		}
		free(lays[ix]);
	}
	free(lays);	
	
	for(ip = 0; ip < 5; ip++){
		for(ix = 0; ix < nx+4; ix++){
			for(iy = 0; iy < ny+4; iy++){
				free(dif_p_x[ip][ix][iy]);
			}
			free(dif_p_x[ip][ix]);
		}
		free(dif_p_x[ip]);	
	}
	free(dif_p_x);
	
	for(ip = 0; ip < 5; ip++){
		for(ix = 0; ix < nx+4; ix++){
			for(iy = 0; iy < ny+4; iy++){
				free(dif_p_y[ip][ix][iy]);
			}
			free(dif_p_y[ip][ix]);
		}
		free(dif_p_y[ip]);	
	}
	free(dif_p_y);
	
	for(ip = 0; ip < 5; ip++){
		for(ix = 0; ix < nx+4; ix++){
			for(iy = 0; iy < ny+4; iy++){
				free(dif_p_z[ip][ix][iy]);
			}
			free(dif_p_z[ip][ix]);
		}
		free(dif_p_z[ip]);	
	}
	free(dif_p_z);
	return;
}

void read_val(char* val_fn){
	FILE* fp;
	char buff[2048];
	char *buf;
	int nbytes;
	int ix, iy, iz,nx, ny, nz, nny, nnx;
	buf = (char*)malloc(sizeof(char)*2048);
	nx = n[0];	ny = n[1];	nz = n[2];
	nny = nz;	nnx = nz*ny;

	if((fp = fopen(val_fn, "r")) == NULL){
		fprintf(stderr, "Cannot open file with grid: %s\n", val_fn);
		free_all();
		exit(1);
	}
	for(ix = 0; ix < nx+2; ix++){
		for(iy = 0; iy < ny+2; iy++){
			fgets(buff, 2048, fp);
			iz = 0;
			buf = buff;
			while(sscanf(buf, "%d%n", &lays[ix][iy][iz], &nbytes) == 1){
				//fprintf(stderr, "%s\n", buf);
				buf += nbytes;
				iz++;
			}
			if(iz != nz+2){
				fprintf(stderr, "Value file %s corrupted: layers %d %d %d\n", val_fn, ix, iy, iz);
				free_all();
				exit(1);
			}
		}
	}
	for(ix = 0; ix < nx+2; ix++){
		for(iy = 0; iy < ny+2; iy++){
			fgets(buff, 2048, fp);
			buf = buff;
			iz = 0;
			while(sscanf(buf, "%lf%n", &p[0][ix][iy][iz], &nbytes) == 1){
				//fprintf(stderr, "%s %d\n", buf, iz);
				buf += nbytes;
				iz++;}
			if(iz != nz+2){
				fprintf(stderr, "Value file %s corrupted: density %d %d\n", val_fn, ix, iy);
				free_all();
				exit(1);
			}
		}
	}
	for(ix = 0; ix < nx+2; ix++){
		for(iy = 0; iy < ny+2; iy++){
			fgets(buff, 2048, fp);
			buf = buff;
			iz = 0;
			while(sscanf(buf, "%lf%n", &p[1][ix][iy][iz], &nbytes) == 1){
				//fprintf(stderr, "%s %s %d\n","p[1]", buf, iz);
				buf += nbytes;
				iz++;}
			if(iz != nz+2){
				fprintf(stderr, "Value file %s corrupted: x-velocity %d %d\n", val_fn, ix, iy);
				free_all();
				exit(1);
			}
		}
	}
	for(ix = 0; ix < nx+2; ix++){
		for(iy = 0; iy < ny+2; iy++){
			fgets(buff, 2048, fp);
			buf = buff;
			iz = 0;
			while(sscanf(buf, "%lf%n", &p[2][ix][iy][iz], &nbytes) == 1){
				//fprintf(stderr, "%s %s %d\n", "p[2]",buf, iz);
				buf += nbytes;
				iz++;}
			if(iz != nz+2){
				fprintf(stderr, "Value file %s corrupted: y-velocity %d %d\n", val_fn, ix, iy);
				free_all();
				exit(1);
			}
		}
	}
	for(ix = 0; ix < nx+2; ix++){
		for(iy = 0; iy < ny+2; iy++){
			fgets(buff, 2048, fp);
			buf = buff;
			iz = 0;
			while(sscanf(buf, "%lf%n", &p[3][ix][iy][iz], &nbytes) == 1){
				//fprintf(stderr, "%s %d\n", "p[3]",buf, iz);
				buf += nbytes;
				iz++;}
			if(iz != nz+2){
				fprintf(stderr, "Value file %s corrupted: z-velocity %d %d\n", val_fn, ix, iy);
				free_all();
				exit(1);
			}
		}
	}
	for(ix = 0; ix < nx+2; ix++){
		for(iy = 0; iy < ny+2; iy++){
			fgets(buff, 2048, fp);
			buf = buff;
			iz = 0;
			while(sscanf(buf, "%lf%n", &p[4][ix][iy][iz], &nbytes) == 1){
				//fprintf(stderr, "%s %s %d\n","p[4]", buf, iz);
				buf += nbytes;
				iz++;}
			if(iz != nz+2){
				fprintf(stderr, "Value file %s corrupted: pressure %d %d\n", val_fn, ix, iy);
				free_all();
				exit(1);
			}
		}
	}
	fclose(fp);
	//free(buf);
	return;
}

int main(){
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	read_grid("KIRPICH.GRD");
	read_val("initval");
	free_all();
	return 0;
}