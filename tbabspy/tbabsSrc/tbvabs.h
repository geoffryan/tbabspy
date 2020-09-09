#ifndef TBVABS_H

#define TBVABS_H

void tbnew(const double *ear, int ne, const double *param, int ifl, double *photar,
		 double *photer, const char *init);
void tbvabsfeo(const double *ear, int ne, const double *param, int ifl, double *photar,
	       double *photer, const char *init);
void tbgas(const double *ear, int ne, const double *param, int ifl, double *photar,
		 double *photer, const char *init);
void tbvabspcf(const double *ear, int ne, const double *param, int ifl, double *photar,
	       double *photer, const char *init);
void tbvabsrel(const double *ear, int ne, const double *param, int ifl, double *photar,
	       double *photer, const char *init);

#endif
