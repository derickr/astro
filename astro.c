#include <math.h>

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


double rev( double x )
{
	return  x - floor( x / 360.0 ) * 360.0;
}


double cbrt( double x )
{
	if ( x > 0.0 )
	{
		return exp( log(x) / 3.0 );
	}
	else if ( x < 0.0 )
	{
		return -cbrt(-x);
	}
	else /* x == 0.0 */
	{
		return 0.0;
	}
}

/*
 * in:   ra and decl in degrees, r in "units"
 * out:  x, y, z in "units"
 */
void spherical_to_rectangular(double ra, double decl, double r, double *x, double *y, double *z)
{
	*x = r * cosd(ra) * cosd(decl);
	*y = r * sind(ra) * cosd(decl);
	*z = r * sind(decl);
}

/*
 * in:  x, y, z in "units"
 * out: ra and decl in degrees, r in "units"
 */
void rectangular_to_spherical(double x, double y, double z, double *ra, double *decl, double *r)
{
	*r = sqrt( x * x + y * y + z * z );
	*ra = atan2d( y, x );
	*decl = atan2d( z, sqrt( x * x + y * y ) );
}


/*
 * in:  x, y, z in "units", oblecl in degrees
 * out: xr, yr, xr in "units"
 */
static void rotate(double x, double y, double z, double oblecl, double *xr, double *yr, double *zr)
{
	*xr = x;
	*yr = y * cosd(oblecl) - z * sind(oblecl);
	*zr = y * sind(oblecl) + z * cosd(oblecl);
}

/*
 * in:  x, y, z in "units", oblecl in degrees
 * out: xr, yr, xr in "units"
 */
static void rotatey(double x, double y, double z, double oblecl, double *xr, double *yr, double *zr)
{
	*xr = x * cosd(oblecl) - z * sind(oblecl);
	*yr = y;
	*zr = x * sind(oblecl) + z * cosd(oblecl);
}

/*
 * in:  x, y, z in "units", oblecl in degrees
 * out: xr, yr, xr in "units"
 */
void ecliptic_to_equatorial(double x, double y, double z, double oblecl, double *xr, double *yr, double *zr)
{
	return rotate(x, y, z, oblecl, xr, yr, zr);
}

/*
 * in:  x, y, z in "units", oblecl in degrees
 * out: xr, yr, xr in "units"
 */
void equatorial_to_ecliptic(double x, double y, double z, double oblecl, double *xr, double *yr, double *zr)
{
	return rotate(x, y, z, 0 - oblecl, xr, yr, zr);
}

/*
 * in:  y, m, d
 * out: daynr
 */
void date_to_daynr(int y, int m, int d, int *daynr)
{
	*daynr = 367 * y - (7 * (y + ((m+9)/12))) / 4 + (275 * m)/9 + d - 730530;
}


void sunpos(int d, double *lon, double *ra, double *decl, double *rad)
{
	double w, a, e, M, oblecl, L, E, x, y, r, v, z, ex, ey, ez;
	double GMST0, SIDTIME, HA;
	double LONG, LAT, azimuth, altitude;

	w = 282.9404 + 4.70935e-5 * d;
	a = 1.000000;
	e = 0.016709 - 1.151e-9 * d;
	M = 356.0470 + 0.9856002585 * d;
	M = rev(M);

	oblecl = 23.4393 - 3.563e-7 * d;

	L = w + M;
	L = rev(L);

	E = M + (180/PI) * e * sind(M) * (1 + e * cosd(M));

	x = cosd(E) - e;
	y = sind(E) * sqrt(1 - e * e);

	r = sqrt(x*x + y*y);
	v = atan2d(y, x);

	*lon = rev(v + w);

	/* reusing x, y, and z here */
	x = r * cosd(*lon);
	y = r * sind(*lon);
	z = 0.0;

	ecliptic_to_equatorial(x, y, z, oblecl, &ex, &ey, &ez);
	rectangular_to_spherical(ex, ey, ez, ra, decl, rad);

	LONG = -0.127;
	LAT  = 51.50;

	GMST0 = rev( L + 180 ); /* in degrees */
	SIDTIME = GMST0 + ( 8.2 * 15 ) /* UT */ + LONG /* 15°E */;

	HA = SIDTIME - *ra;

	/* reusing x, y, and z here again */
	spherical_to_rectangular(HA, *decl, 1, &x, &y, &z);
	rotatey(x, y, z, 90 - LAT, &ex, &ey, &ez);

	azimuth = atan2d( ey, ex ) + 180;
	altitude = atan2d( ez, sqrt( ex*ex + ey*ey ) );

	printf("%f, %f\n", azimuth, altitude);
}



int main(void)
{
/*
	double x, y, z, ex, ey, ez, ra, decl, r;

	spherical_to_rectangular(90, 0, 1, &x, &y, &z);
	ecliptic_to_equatorial(x, y, z, 23.4, &ex, &ey, &ez);
	rectangular_to_spherical(ex, ey, ez, &ra, &decl, &r);

	printf("%f, %f, %f\n", ra, decl, r);
*/
	int daynr;
	double lon, ra, decl, r;

	date_to_daynr(2013, 1, 11, &daynr);

	sunpos(daynr, &lon, &ra, &decl, &r);
}
