#include <math.h>
#include "astro.h"

#define PI          3.14159265358979323846
#define RADEG       (180.0/PI)
#define DEGRAD      (PI/180.0)
#define sind(x)     sin((x)*DEGRAD)
#define cosd(x)     cos((x)*DEGRAD)
#define tand(x)     tan((x)*DEGRAD)
#define asind(x)    (RADEG*asin(x))
#define acosd(x)    (RADEG*acos(x))
#define atand(x)    (RADEG*atan(x))
#define atan2d(y,x) (RADEG*atan2((y),(x)))

/* Limits the angle "x" to 360 degrees */
static double rev( double x )
{
	return  x - floor( x / 360.0 ) * 360.0;
}

/* Converts spherical coordinates ra, decl and r to rectangular coordinates x,
 * y, and z
 *
 * in:   ra and decl in degrees, r in "units"
 * out:  x, y, z in "units"
 */
void spherical_to_rectangular(double ra, double decl, double r, double *x, double *y, double *z)
{
	*x = r * cosd(ra) * cosd(decl);
	*y = r * sind(ra) * cosd(decl);
	*z = r * sind(decl);
}

/* Converts rectangular coordinates x, y and z to spherical coordinates ra,
 * decl and r
 *
 * in:  x, y, z in "units"
 * out: ra and decl in degrees, r in "units"
 */
void rectangular_to_spherical(double x, double y, double z, double *ra, double *decl, double *r)
{
	*r = sqrt( x * x + y * y + z * z );
	*ra = atan2d( y, x );
	*decl = atan2d( z, sqrt( x * x + y * y ) );
}

/* Rotates the rectangular coordinates x, y, and z around the X-axis by oblecl
 * degrees
 *
 * in:  x, y, z in "units", oblecl in degrees
 * out: xr, yr, xr in "units"
 */
static void rotatex(double x, double y, double z, double oblecl, double *xr, double *yr, double *zr)
{
	*xr = x;
	*yr = y * cosd(oblecl) - z * sind(oblecl);
	*zr = y * sind(oblecl) + z * cosd(oblecl);
}

/* Rotates the rectangular coordinates x, y, and z around the Y-axis by oblecl
 * degrees
 *
 * in:  x, y, z in "units", oblecl in degrees
 * out: xr, yr, xr in "units"
 */
static void rotatey(double x, double y, double z, double oblecl, double *xr, double *yr, double *zr)
{
	*xr = x * cosd(oblecl) - z * sind(oblecl);
	*yr = y;
	*zr = x * sind(oblecl) + z * cosd(oblecl);
}

/* Converts ecliptic to equatorial coordinates by rotating around the X-axis
 *
 * in:  x, y, z in "units", oblecl in degrees
 * out: xr, yr, xr in "units"
 */
void ecliptic_to_equatorial(double x, double y, double z, double oblecl, double *xr, double *yr, double *zr)
{
	return rotatex(x, y, z, oblecl, xr, yr, zr);
}

/* Converts equatorial to ecliptic coordinates by reverse-rotating around the X-axis
 *
 * in:  x, y, z in "units", oblecl in degrees
 * out: xr, yr, xr in "units"
 */
void equatorial_to_ecliptic(double x, double y, double z, double oblecl, double *xr, double *yr, double *zr)
{
	return rotatex(x, y, z, 0 - oblecl, xr, yr, zr);
}

/* Computes year, month, day, hour, minute and second to a "day number" from
 * 2000 Jan 0.0
 *
 * in:  y, m, d, h, i, s
 * out: time
 */
void date_to_daynr(int y, int m, int d, int h, int i, int s, double *time)
{
	int daynr;

	daynr = 367 * y - (7 * (y + ((m+9)/12))) / 4 + (275 * m)/9 + d - 730530;
	*time = daynr + ((h + ((i + (s / 60)) / 60)) / 24);
}

/* Calculates the Sun's Mean Longitude (L), Right Ascension (RA), Declination
 * (decl) and "rad"
 *
 * in: d (time)
 * out: L, RA, decl, rad */
void sunpos(double d, double *L, double *M, double *ra, double *decl, double *rad)
{
	double w, a, e, oblecl, E, x, y, r, v, z, ex, ey, ez, i_lon;

	w = 282.9404 + 4.70935e-5 * d;
	a = 1.000000;
	e = 0.016709 - 1.151e-9 * d;
	*M = 356.0470 + 0.9856002585 * d;
	*M = rev(*M);

	oblecl = 23.4393 - 3.563e-7 * d;

	*L = w + *M;
	*L = rev(*L); /* Sun's mean logitude */

	E = *M + (180/PI) * e * sind(*M) * (1 + e * cosd(*M));

	x = cosd(E) - e;
	y = sind(E) * sqrt(1 - e * e);

	r = sqrt(x*x + y*y);
	v = atan2d(y, x);

	i_lon = rev(v + w);

	/* reusing x, y, and z here */
	x = r * cosd(i_lon);
	y = r * sind(i_lon);
	z = 0.0;

	ecliptic_to_equatorial(x, y, z, oblecl, &ex, &ey, &ez);
	rectangular_to_spherical(ex, ey, ez, ra, decl, rad);
}

void moonpos(double d, double Ls, double Ms, double *ra, double *decl, double *rad)
{
	double N, i, w, a, e, M;
	double E, E0, E1;
	double x, y, z, r, v, ex, ey, ez, lon, lat;
	double Lm, Mm, D, F;

	N = rev(125.1228 - 0.0529538083 * d);
	i = rev(  5.1454);
	w = rev(318.0634 + 0.1643573223 * d);
	a = rev( 60.2666);
	e = rev(  0.054900);
	M = rev(115.3654 + 13.0649929509 * d);

	E0 = M + (180/M_PI) * e * sind(M) * (1 + e * cosd(M));
	do {
		E1 = E0 - (E0 - (180/M_PI) * e * sind(E0) - M) / (1 - e * cosd(E0));
		if (fabs(E1 - E0) < 0.0005) {
			break;
		}
		E0 = E1;
	} while (1);
	E = E1;

	x = a * (cosd(E) - e);
	y = a * sqrt(1 - e*e) * sind(E);

	r = sqrt(x*x + y*y);
	v = rev(atan2d(y, x));

	ex = r * (cosd(N) * cosd(v + w) - sind(N) * sind(v + w) * cosd(i));
	ey = r * (sind(N) * cosd(v + w) + cosd(N) * sind(v + w) * cosd(i));
	ez = r * sind(v + w) * sind(i);

	rectangular_to_spherical(ex, ey, ez, &lon, &lat, rad);

	Lm = rev(N + w + M);
	Mm = M;
	D = rev(Lm - Ls);
	F = rev(Lm - N);

	lon = lon
		- 1.274 * sind(Mm - 2*D)
		+ 0.658 * sind(2 * D)
		- 0.186 * sind(Ms)
		- 0.059 * sind(2*Mm - 2*D)
		- 0.057 * sind(Mm - 2*D + Ms)
		+ 0.053 * sind(Mm + 2*D)
		+ 0.046 * sind(2*D - Ms)
		+ 0.041 * sind(Mm - Ms)
		- 0.035 * sind(D)
		- 0.031 * sind(Mm + Ms)
		- 0.015 * sind(2*F - 2*D)
		+ 0.011 * sind(Mm - 4*D);

	lat = lat
		- 0.173 * sind(F - 2*D)
		- 0.055 * sind(Mm - F - 2*D)
		- 0.046 * sind(Mm + F - 2*D)
		+ 0.033 * sind(F + 2*D)
		+ 0.017 * sind(2*Mm + F);

	*rad = *rad
		- 0.58 * cosd(Mm - 2*D)
		- 0.46 * cosd(2*D);
	
	spherical_to_rectangular(lon, lat, 1, &x, &y, &z);
	ecliptic_to_equatorial(x, y, z, (23.4393 - 3.563e-7 * d), &ex, &ey, &ez);
	rectangular_to_spherical(ex, ey, ez, ra, decl, rad);
	*ra = rev(*ra);
}

/* Calculates the Siderial Time (sidtime) and Hour Angle (HA)
 * in: L, utoffset, lon, ra
 * out: sidtime, ha
 */
void sidtime_and_ha(double L, double utoffset, double lon, double ra, double *sidtime, double *ha)
{
	double GMST0;

	/* Sidereal time and hour angle, Altitude and azimuth */
	GMST0 = rev( L + 180 ); /* in degrees */
	*sidtime = GMST0 + ( utoffset * 15 ) /* UT */ + lon /* 15°E */;

	*ha = *sidtime - ra;
}

/* Calculates the Sun's Altitude and Azimuth
 * in: L, lat, ha, decl
 * out: azimuth, altitude */
void sunaltazimuth(double L, double lat, double ha, double decl, double *azimuth, double *altitude)
{
	double x, y, z;
	double ex, ey, ez; /* temp vars */

	spherical_to_rectangular(ha, decl, 1, &x, &y, &z);
	rotatey(x, y, z, 90 - lat, &ex, &ey, &ez);

	*azimuth = atan2d( ey, ex ) + 180;
	*altitude = atan2d( ez, sqrt( ex*ex + ey*ey ) );
}
