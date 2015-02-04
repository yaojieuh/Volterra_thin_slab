void greenfunction(double *r1, double *r2, double k,  double* G0);
void greenfunction2(double x, double z, double k,  double* G00);
void prodcomp( double *z1, double *z2, double *zres );
void greenfunction3(double x, double z, double k,  double* G0);
void Pincforward( int iw, int nx, double dx, int nk, double dk, double k, int iz,double dz, double *P1r,double *P1i, double *P0r,double *P0i);
void Greentable( int iw, int nx2, double dx,double z, double k,    double *G0r,double *G0i);
void P1num(FILE *filefprintf, int nw, int nx, int nz,double dx,double dz,double c0, int *ps, double *sourcefren,double* fren, double* vpe, double *P1r,double *P1i);
void Pwtot(int nw, int nx,  int nz,   double dw,double dx, double dz, double dt, int nt,  double* fren, double *P1r, double *P1i, double *Pret);
void greenfunctionV(double x, double z, double k,  double* G00);
