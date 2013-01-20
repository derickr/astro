#include "astro.h"

int main(void)
{
	double daynr;
	double lat, lon;
	int    year, month, day, h, m;
	double L, M, ra, decl, rad, SIDTIME, HA;
	double azimuth, altitude;

	year = 2013;
	month = 1;
	day = 12;
	h = 8;
	m = 0;
	lat = 51.50;
	lon = -0.127;

	year = 1990;
	month = 4;
	day = 19;
	h = 0;
	m = 0;
	lat = 60;
	lon = 15;

	date_to_daynr(year, month, day, h, m, 0, &daynr);

	sunpos(daynr, &L, &M, &ra, &decl, &rad);
	printf("%f, %f, %f, %f\n", L, ra, decl, rad);

	moonpos(daynr, L, M, &ra, &decl, &rad);
	printf("MOON: %f, %f, %f\n", ra, decl, rad);

	sidtime_and_ha(L, h + (m/60), lon, ra, &SIDTIME, &HA);

	sunaltazimuth(L, lat, HA, decl, &azimuth, &altitude);
	printf("%f, %f\n", azimuth, altitude);
}

