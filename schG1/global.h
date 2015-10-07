#pragma once
#include <stdio.h>
#pragma once
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <memory.h>

#ifndef __GLOBAL__
#define __GLOBAL__
extern int n[3];
extern int ***lays;
extern double ****p;			/*(nx+4)*(ny+4)*(nz+4)*/
extern double ****dif_p_x;
extern double ****dif_p_y;
extern double ****dif_p_z;

extern double **coord;
extern double **dif;

#endif
#define CHECKERR(ierr) assert(!(ierr))
#define PARMS_malloc(base, nmem, type) {\
    (base) = (type *)malloc((nmem)*sizeof(type)); \
    CHECKERR((base) == NULL); \
}
