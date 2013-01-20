#include <stdio.h>

#ifndef __ASTRO_H__
#define __ASTRO_H__
void date_to_daynr(int y, int m, int d, int h, int i, int s, double *time);
void sunpos(double d, double *L, double *M, double *ra, double *decl, double *rad);
void moonpos(double d, double Ls, double Ms, double *ra, double *decl, double *rad);
void sidtime_and_ha(double L, double utoffset, double lon, double ra, double *sidtime, double *ha);
void sunaltazimuth(double L, double lat, double ha, double decl, double *azimuth, double *altitude);
#endif
