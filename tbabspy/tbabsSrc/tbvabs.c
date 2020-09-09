// define LINEARCACHE to have the cache in a linear array
#define LINEARCACHE

#ifdef LINEARCACHE
#define EXACT 0
#define INTERVAL 1
#endif

// define TREECACHE to store the cache in a (balanced, on Linux) tree
// #define TREECACHE

#include "tbvabs.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef TREECACHE
#include <search.h>
#endif

#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


////////////////////////////////////////////////////////////////
typedef struct {
  double e;
  double *sigma;
} cachesig;


double horner(double e,double *coeff,int n);

double vernabs(double e,int z);

double habs(int z,double en);
double heabs(double en);
double h2abs(double en);

double oabs(double en);
double neabs(double en);
double feabs(double en);
double siginter(int nstep,double *energ,double *sig, double en);
double photabs(int zind,double en);

double jwfano(double q,double nu,double gamma,double ener);


#ifdef TREECACHE
int cachecompare(const void *left, const void *right);
#endif

cachesig *sc_new(double en);
void sc_free(cachesig *sig);
int sc_insert(cachesig *sig);
int sc_binsrch(double e,int precision);
void sc_dump(int zind);

double dgami(double x,double y);
double fgabnd(const char *ele);

int binsrch(double e,double *ener,int n);

double photavg(int zind,double elo,double ehi);



//////////////////////////////////////////////////////////////////////////////
// Fortran routines from XSPEC
double dgami_(double *x,double *y);

// 
float FGABND(char *ele);

void phfit2_(int *z,int *ion,int *shell,float *energ, float *sigma);


//////////////////////////////////////////////////////////////////////////////

// Mbarn -> cm^2
#define MBARN2CM2 1e-18

// Atomic mass unit in grams
#define AMU 1.6605655E-24

// number of elements
#define NELE 18

// number of species cached (=NELE plus whatever molecules etc. we use)
#define NION 19

// elements taken into account
// Missing compared to the paper: P, Ti, Mn
// 1000: H2
const int atomic[]={1,2,6,7,8,10,11,12,13,14,16,17,18,20,24,26,27,28,1000};

const char *celts[]={"H", "He", "C", "N", "O", "Ne", "Na", "Mg",
		     "Al", "Si", "S", "Cl", "Ar", "Ca", "Cr", "Fe",
		     "Co", "Ni", "H2"};
const double weight[]={  1.,  4.,  12., 14.,  16.,  20.,  23.,  24.,
			 27., 28.,  32., 35.,  40.,  40.,  52.,  56.,
			 59., 59., 2.};


// Initial and maximum cache size
#define MAXCACHE 0xFFFFF

#ifdef LINEARCACHE
#define STARTSIZE 1024
#endif

// cache for the cross section values
#ifdef TREECACHE
void *sigcache=NULL;
size_t cachenum=0; // number of elements in the cache

#endif
#ifdef LINEARCACHE
cachesig *(*sigcache)=NULL;
size_t cachesize=0; // current size of linear cache. This always points at the
                    // next free element
size_t cachelen=0;

#endif

////////////////////////////////////////////////////////////////

void tbnew(const double *ear, int ne, const double *param, int ifl, double *photar,
	    double *photer, const char *init) {
  //
  //     Compute X-ray absorptivity of cold gas using the formalism given
  //     by Wilms, Allen, and McCray, 2000, ApJ 542, 914-924
  //
  //     The parameters are
  //
  //     param( 0)   hydrogen column 
  //                 positive: column in 1e22/cm^2
  //                 negative: column in 1/cm^2
  //
  //     param(1..29): 
  //          POSITIVE: relative abundance of the element wrt the currently
  //                    set abundance table
  //          NEGATIVE: abs(param(i)) is the column of the element (1/cm^2)
  //                    (makes only sense for param(1)>0)
  //     param( 1)   helium      (2)
  //     param( 2)   carbon      (6)
  //     param( 3)   nitrogen    (7)
  //     param( 4)   oxygen      (8)
  //     param( 5)   neon        (10)
  //     param( 6)   sodium      (11)
  //     param( 7)   magnesium   (12)
  //     param( 8)   aluminum    (13)
  //     param( 9)   silicon     (14)
  //     param(10)   sulphur     (16)
  //     param(11)   chlorine    (17)
  //     param(12)   argon       (18)
  //     param(13)   calcium     (20)
  //     param(14)   chromium    (24)
  //     param(15)   iron        (26)
  //     param(16)   cobalt      (27)
  //     param(17)   nickel      (28)
  //     param(18): fraction of hydrogen existing in molecular form
  //     param(19): density of dust in g/cm^3
  //        if param(10) <= 0 then no dust is included in the modeling
  //        (i.e., the depletion factor is also not taken into account)
  //     param(20): minimum thickness of dust (mum)
  //     param(21): maximum thickness of dust (mum)
  //     param(22): PL index for dust size distribution (implicit minus-sign!)
  //        if param(20) = param(21): use single size grain
  //
  //     param(23..41): depletion factor (ratio of GAS to total ISM abundance)
  //     param(23)   hydrogen    (1)
  //     param(24)   helium      (2)
  //     param(25)   carbon      (6)
  //     param(26)   nitrogen    (7)
  //     param(27)   oxygen      (8)
  //     param(28)   neon        (10)
  //     param(29)   sodium      (11)
  //     param(30)   magnesium   (12)
  //     param(31)   aluminium   (13)
  //     param(32)   silicon     (14)
  //     param(33)   sulphur     (16)
  //     param(34)   chlorine    (17)
  //     param(35)   argon       (18)
  //     param(36)   calcium     (20)
  //     param(37)   chromium    (24)
  //     param(38)   iron        (26)
  //     param(39)   cobalt      (27)
  //     param(40)   nickel      (28)
  //
  //     param(41):  Redshift
  

  int i;
  
  double abun[NELE];
  double beta[NELE];

  static int startup=1;

  if (startup) {
    //
    // print out the copyright (once)
    //

    fprintf(stderr,"tbvabs Version 2.3\n");
    fprintf(stderr,"Cosmic absorption with grains and H2, modified from\n");
    fprintf(stderr,"Wilms, Allen, & McCray, 2000, ApJ 542, 914-924\n");
    fprintf(stderr,"Questions: Joern Wilms\n");
    fprintf(stderr,"joern.wilms@sternwarte.uni-erlangen.de\n");
    fprintf(stderr,"joern.wilms@fau.de\n");
    fprintf(stderr,"http://pulsar.sternwarte.uni-erlangen.de/wilms/research/tbabs/\n\n");
    fprintf(stderr,"PLEASE NOTICE:\n");
    fprintf(stderr,"To get the model described by the above paper\n");
    fprintf(stderr,"you will also have to set the abundances:\n");
    fprintf(stderr,"   abund wilm\n\n");
    fprintf(stderr,"Note that this routine ignores the current cross section setting\n");
    fprintf(stderr,"as it always HAS to use the Verner cross sections as a baseline.\n");
    startup=0;
  }
  
  // Put parameters into more meaningful variables
  // column
  double nh=param[0];

  // this variable is always a column
  if (nh<0) {
    nh=fabs(nh)/1e22;
  }

  // because of the way how we normalize, we need at
  // least one H-atom in our way.
  if (nh<1e-22) {
    nh=1e-22;
  }

  // H2
  double h2=param[18];
  if (h2<0.) {
    fprintf(stderr,"The H2 fraction is negative. Setting to 0.\n");
    h2=0.;
  }
  if (h2>1.) {
    fprintf(stderr,"The H2 fraction is greater than 1. Setting to 1.\n");
    h2=1.;
  }

  // grains
  // amin,amax: min and max thickness in cm!
  // pl: power law index
  // rho: density
  double rho=param[19];
  double amin=param[20]*1E-4;
  double amax=param[21]*1E-4;
  double pl=param[22];

  // set the elemental abundances
  abun[0]=1.;
  double nh22=nh*1e22;
  for (i=1; i<NELE; i++) {
    if (param[i]>0) {
      // positive abundance: measure relative to H abundance
      // of current data set. So things are easy
      abun[i]=param[i]*fgabnd(celts[i]);
    } else {
      // negative abundance: parameter contains the COLUMN of this
      // element, so its relative abundance with respect to H
      // is just the ratio of param(i) to nh:
      abun[i]=fabs(param[i])/nh22;
    }

  }

  // Fraction of atoms in the grain phase
  if (rho>0) {
    for (i=0; i<NELE; i++) {
      beta[i]=1.0-param[23+i];
      if (beta[i]<0.) {
	fprintf(stderr,"Limiting gas/dust to 0 for Z=%i (%s)\n",atomic[i],celts[i]);
	beta[i]=0.;
      }
      if (beta[i]>1.) {
	fprintf(stderr,"Limiting gas/dust to 1 for Z=%i (%s)\n",atomic[i],celts[i]);
	beta[i]=1.;
      }
    }
  } else {
    // no grains
    for (i=0; i<NELE; i++) {
      beta[i]=0.;
    }
  }

  double *siggas=calloc(ne,sizeof(double));
  if (siggas==NULL) {
    fprintf(stderr,"Out of memory: siggas\n");
    exit(102);
  }
  double *siggrains=calloc(ne,sizeof(double));
  if (siggrains==NULL) {
    fprintf(stderr,"Out of memory: siggrains\n");
    exit(102);
  }
  double *sigmol=calloc(ne,sizeof(double));
  if (sigmol==NULL) {
    fprintf(stderr,"Out of memory: sigmol\n");
    exit(102);
  }
  double *sigavg=calloc(ne,sizeof(double));
  if (sigavg==NULL) {
    fprintf(stderr,"Out of memory: sigavg\n");
    exit(102);
  }

  // Redshift: move energies into source frame
  double *earz=malloc((ne+1)*sizeof(double));
  if (earz==NULL) {
    fprintf(stderr,"Out of memory: earz\n");
    exit(102);
  }
  double zfac1000=1000.*(1.+param[41]);
  for (i=0;i<=ne;i++) {
    earz[i]=ear[i]*zfac1000; // convert to rest frame and eV
  }

  
  /////////////////////////////////////////////////////////////

  // mean weights for dust
  double normal=0.; // normalization
  double massgrain=0.;  // mass of grains

  int zind;
  for (zind=0; zind<NELE; zind++) {

    // abundance of this element in dust and gas
    double abdust=abun[zind]*beta[zind];
    double abgas=abun[zind]*(1.-beta[zind]);
    // take care of molecular H
    if (atomic[zind]==1) {
      abgas=abgas*(1.-h2);
    }

    // dust propeties
    massgrain+=abdust*weight[zind];
    normal+=abdust;
    
    //  sum total cross section for the dust and gas phase
    for (i=0; i<ne; i++) {
      double phsig=photavg(zind,earz[i],earz[i+1]);
      siggas[i]+=abgas*phsig;
      sigavg[i]+=abdust*phsig;
    }
  }

  //
  // absorption in the grain phase
  //
  if (rho>0. && normal>0. ) {
    // convert mass to g
    massgrain*=AMU;

    // Mean weight per atom in g
    double meanmol=massgrain/normal;
        
    // normalize cross section
    for (i=0;i<ne; i++) {
      sigavg[i]=sigavg[i]/normal;
    }
    
    // Shielding calculation
    if (amin==amax) {
      //
      // single sized grains
      //
      double a=amin;
      
      // Number of grains per H atom
      double ngrains=massgrain/(rho*(4./3.)*M_PI*pow(a,3.));
      
      // Column of each grain (atoms/cm^2)
      double nbar=(4./3.)*rho*a/meanmol;
      
      // grain cross-section per H atom
      double fac=ngrains*M_PI*a*a;
      for (i=0; i<ne; i++) {
	siggrains[i]=fac*(1.-exp(-sigavg[i]*nbar));
      }
    } else {
      //
      // Case of grains with a power-law distribution
      //
      
      // Normalization of PL Distribution
      double k;
      if (pl==1.) {
	k=1./log(amax/amin);
      } else {
	k=(1.-pl) / ( pow(amax, 1.-pl) - pow(amin, 1.-pl) );
      }
      
      // Number of grains per H atom
      double ngrains=massgrain/(rho*(4./3.)*M_PI*k);
      
      if (pl==4.) {
	ngrains=ngrains/log(amax/amin);
      } else {
	ngrains=ngrains*(4.-pl)/( pow(amax, 4.-pl) - pow(amin, 4.-pl ));
      }
      
      //
      // Formulae from the appendix of the Wilms et al. paper
      // incl. conversion to Mbarns
      //
      double pmin=pow(amin,3.-pl);
      double pmax=pow(amax,3.-pl);
      double signorm=ngrains* M_PI *k/(pl-3.);
      double sigd=(4./3.)*rho/meanmol;
      for (i=0; i<ne; i++) {
	double sigdummy=sigd*sigavg[i];
	double sigmax=sigdummy*amax;
	double sigmin=sigdummy*amin;
	siggrains[i]=signorm* (pmin*(1.-exp(-sigmin)) - pmax*(1.-exp(-sigmax)) + 
			       pow(sigdummy,pl-3.) * ( 
						      dgami(4.-pl,sigmax) - dgami(4.-pl,sigmin)
						       ) 
			       );
      }
    }
  }
  
  //
  // Molecular phase
  // (right now we have H2 only, but other molecules can be added)
  //
  if (h2>0.) {
    // ... not that H will ever be depleted
    double ab=abun[0]*(1.-beta[0]);
    // abundance of H2 by number (two H go into one H2)
    ab=ab*h2*0.5;
    for (i=0;i<ne;i++){
      sigmol[i]=photavg(18,earz[i],earz[i+1])*ab;
    }
  }
  
  // add everything up 
  for (i=0; i<ne; i++) {
    photar[i]=exp(-nh22*(siggas[i]+siggrains[i]+sigmol[i]));
  }

  free(siggas);
  free(siggrains);
  free(sigmol);
  free(sigavg);
  free(earz);
  
}

double oabs(double en) {
  //
  // cross section of oxygen
  //

#include "ophoto.h"

  int nstep=sizeof(energ)/sizeof(energ[0]);

  //  outside of the edge or very coarse response?
  //  if yes: use verner cross section
  if (en<energ[0] || en>energ[nstep-1]) {
    // yes: use Verner cross section
    return vernabs(en,8);
  } 
  // no: interpolate in cross section table
  return siginter(nstep,energ,sig,en);
}

double neabs(double en) {
  //
  // cross section of Neon
  //

#include "nephoto.h"

  int nstep=sizeof(energ)/sizeof(energ[0]);

    //  outside of the edge or very coarse response?
  //  if yes: use verner cross section
  if (en<energ[0] || en>energ[nstep-1]) {
    // yes: use Verner cross section
    return vernabs(en,10);
  } 
  // no: interpolate in cross section table
  return siginter(nstep,energ,sig,en);
}

double feabs(double en) {
  //
  // cross section for iron
  //

  #include "fephoto.h"

  int nstep=sizeof(energ)/sizeof(energ[0]);

  //  outside of the edge?
  if (en<energ[0] || en>energ[nstep-1]) {
    // yes: use Verner cross section / power law towards verner
    if ( en>=745. &&  en<=795. ) {
      return 1.645E-18*pow(en/720.74, -3.2);
    }
    return vernabs(en,26);
  }

  // no: interpolate in cross section table
  return siginter(nstep,energ,sig,en);
}

double siginter(int nstep,double *energ,double *sig, double en) {
  //
  //  Linear interpolation in a general cross section table
  //

  assert(en>energ[0] && en<energ[nstep-1]);

  // search for the proper energy bin
  int ii=binsrch(en,energ,nstep);
  assert(energ[ii]<=en && en<energ[ii+1]);
  double sum=sig[ii]+
    (en-energ[ii])/(energ[ii+1]-energ[ii])*(sig[ii+1]-sig[ii]);
  return sum*MBARN2CM2;
}


double vernabs(double e,int z) {
  //
  // calculate cross section per Verner et al.
  // returns sigma in cm^2
  //
  int shell;
  
  float ee=e;

  // always loop over all shells, since phfit2_ returns 0 very quickly
  // when called with an invalid shell
  double sigma=0.;
  for (shell=1;shell<=7;shell++) {
    float sig;
    phfit2_(&z,&z,&shell,&ee,&sig);
    sigma+=sig;
  }
  return(sigma*MBARN2CM2);
}


int cachecompare(const void *left, const void *right) {
  //
  // Comparison function for tsearch and tfind
  //
  double eleft = ((cachesig *) left)->e;
  double eright= ((cachesig *) right)->e;
  
  if (eleft<eright) return -1;
  if (eleft>eright) return +1;
  return 0;
}

double photabs(int zind,double en) {
  //
  // return the photoabsorption cross section of element atomic[zind] at energy en
  // energies are measured in eV
  //
  // z<1000: nuclear charge of the element in the gas phase
  // z>=1000: molecules.
  //      1000: H2
  //

  int i;
  
#ifdef TREECACHE
  static cachesig *oldcall=NULL;
#endif
  
#ifdef LINEARCACHE
  static int lastndx=-1;
#endif
  
  if (en>=100E3) {
    return 0.;
  }
  
#ifdef TREECACHE
  if (cachenum>0) {
    if (oldcall!=NULL && (oldcall->e==en)) {
      return oldcall->sigma[zind];
    }
    cachesig findme;
    findme.e=en;
    
    cachesig **cachel= tfind(&findme,&sigcache,cachecompare);
    
    // element is already in the cache
    if (cachel!=NULL) {
      oldcall=*cachel;
      return (*cachel)->sigma[zind];
    }
  }
#endif
    
#ifdef LINEARCACHE
  if (cachesize>0) {
    if (sigcache[lastndx]->e==en) {

      //  printf("ONE %f %e %i\n",en,sigcache[lastndx]->sigma[zind],lastndx);
      return sigcache[lastndx]->sigma[zind];
    }
    // because this routine is generally called from photavg, we will typically
    // call it in a sequence of en(i),en(i+1), etc., 
    lastndx++;
    if (lastndx<cachesize) {
      if (sigcache[lastndx]->e==en) {
	//  printf("TWO %f %e %i\n",en,sigcache[lastndx]->sigma[zind],lastndx);
	return sigcache[lastndx]->sigma[zind];
      }
    }
    
    // not found in previous call, perform a binary search
    int ndx=sc_binsrch(en,EXACT);
    if(ndx!=-1) {
      lastndx=ndx;
      //  printf("FOUND %f %e %i\n",en,sigcache[lastndx]->sigma[zind],lastndx);

      return sigcache[lastndx]->sigma[zind];
    }
  }
#endif

  // cache empty or energy bin not found in cache:
  // calculate cross sections and put them into the cache

  // create the new cache element
  cachesig *sig=sc_new(en);

  // default: calculate all cross sections for this energy
  int istart=0;
  int istop=NION;

  // Cache is full:  just calculate the desired cross section
#ifdef LINEARCACHE
  if (cachesize==MAXCACHE) {
    istart=zind;
    istop=istart+1;
  }
#endif
#ifdef TREECACHE
  if (cachenum==MAXCACHE) {
    istart=zind;
    istop=istart+1;
  }
#endif


  // loop over the ions for which want to calculate the cross section
  for (i=istart; i<istop; i++) {
    switch (atomic[i]) {
    case 1:
      sig->sigma[i]=habs(1,en);
      break;
    case 2:
      sig->sigma[i]=heabs(en);
      break;
    case 8:
      sig->sigma[i]=oabs(en);
      break;
    case 10:
      sig->sigma[i]=neabs(en);
      break;
    case 26:
      sig->sigma[i]=feabs(en);
      break;
    case 1000:
      sig->sigma[i]=h2abs(en);
      break;
    default:
      sig->sigma[i]=vernabs(en,atomic[i]);
      break;
    }
  }

  // the wanted cross section 
  double sigma=sig->sigma[zind];

  
#ifdef TREECACHE
  cachesig **cmp=NULL;
  if (cachenum<MAXCACHE) {
    // insert new element into cache
    // returns NULL if the insertion failed
    cmp=tsearch(sig,&sigcache,cachecompare);
    cachenum++;
  }
  // no insertion (either because we are full or memory is exhausted)
  // -> free the element
  if (cmp==NULL) {
    sc_free(sig);
  }
#endif

#ifdef LINEARCACHE
  // insert element into cache, free it in case
  // insertion was not successful
  lastndx=sc_insert(sig);
  if (lastndx==-1) {
    sc_free(sig);
  }
#endif

  return(sigma);
}

cachesig *sc_new(double en) {
  // initialize a new cache element
  cachesig *sig=malloc(sizeof(cachesig));
  if (sig==NULL) {
    fprintf(stderr,"Out of memory in caching. Sorry.\n");
    exit(100);
  }
  sig->sigma=calloc(NION,sizeof(double));
  if (sig->sigma==0) {
    fprintf(stderr,"Out of memory in caching. Sorry.\n");
    exit(100);
  }
  sig->e=en;
  return(sig);
}

void sc_free(cachesig *sig) {
  // empty a cache element
  free(sig->sigma);
  free(sig);
}

#ifdef LINEARCACHE
int sc_binsrch(double e,int precision) {
  //
  // search for energy e in the cache array
  // If precision==EXACT: and return bin k for which e[k]==e or -1
  // If precision==INTERVAL: return bin k for which e[k]<=e<e[k+1]
  //   (i.e., -1 if e<ener[0])
  int klo=0;
  int khi=cachesize-1;
  while (klo<=khi) {
    int m=((unsigned int) klo+(unsigned int) khi)/2;
    double miden=sigcache[m]->e;
    if (sigcache[m]->e==e) {
      return(m);
    }
    if (miden<e) {
      klo=m+1;
    } else {
      if (miden>e) {
	khi=m-1;
      } else {
	return m;
      }
    }
  }

  if(precision==EXACT) {
    return(-1);
  }
  return khi;
}

int sc_insert(cachesig *sig) {
  // insert element into the linear cache
  // returns the index where the element was inserted or -1 in case we
  // were not successful
  int i;
  if (cachelen>=MAXCACHE) return(-1);
  
  // grow cache?
  // NB this also serves as the initialization of the cache!
  if (cachelen==cachesize) {
    cachelen=(cachelen==0) ? STARTSIZE : cachelen*1.5;
    if (cachelen>MAXCACHE) {
      cachelen=MAXCACHE;
    }
    
    sigcache=realloc(sigcache, (size_t) cachelen*sizeof(sigcache[0]) );
    if (sigcache==NULL) {
      fprintf(stderr,"Cannot grow or initialize linear cache\n");
      exit(102);
    }
  }

  // find insertion point
  int ndx=sc_binsrch(sig->e,INTERVAL);

  // this should never happen when called from inside tbvabs
  if (ndx>-1 && sigcache[ndx]->e == sig->e) {
    printf("This should never happen!\n");
    return (ndx);
  }

  // move elements out of the way
  for (i=cachesize-1;i>ndx;i--) {
    sigcache[i+1]=sigcache[i];
  }
  sigcache[ndx+1]=sig;

  cachesize++;

  return(ndx+1);
}

void sc_dump(int zind) {
  //
  // debugging function: dump the cache for index zind
  //
  int i;
  
  if (sigcache==NULL) {
    return;
  }
  printf("CACHE VALUES FOR INDEX %i (z=%02i [%s]):\n",zind,atomic[zind],celts[zind]);
  for ( i=0;i<cachesize;i++) {
    printf("e: %10.3f s: %13.5e\n",sigcache[i]->e,sigcache[i]->sigma[zind]);
  }
}

#endif

double photavg(int zind,double elo,double ehi) {
  //
  // calculate average photoabsorption cross section between elo and ehi for
  // element atomic[zind]
  //

  // at the moment we naively average at the energy boundaries
  // we will have to do this in a more elegant way in the area of resonances...

  // note that even though this routine will be called for [e(i),e(i+1)] and
  // [e(i+1),e(i+2)], we're still fast because of caching

  double sig1=photabs(zind,elo);
  double sig2=photabs(zind,ehi);

  return((sig1+sig2)/2.);
}


double habs(int z,double en) {
  //
  // Photoionization cross-section for hydrogenic ions
  // after Band et al., 1990, A&A 237, 267
  //
  // e: energy in eV
  // sig: cross-section in cm^2
  //
  // Version 1.0:
  //   J. Wilms, wilms@astro.uni-tuebingen.de, 1999/04/21
  // Version 2.0:
  //   J. Wilms, joern.wilms@sternwarte.uni0erlangen.de, 2016/03/23
  //   conversion to C, now returns sigma in cm^2
  //

  // ... threshold
  double e1s=13.599*z*z;
  if (en<e1s) {
    return 0.;
  }

  double y=2.*en/e1s;
  double denom=z*z*pow(y,1.5)*pow(1.+sqrt(y),4.);
  double sigma=605./denom;

  return(sigma*MBARN2CM2);
}

double heabs(double ener) {
  //
  // Photoionization cross-section for helium after
  // Yan et al., 1998, ApJ 496, 1044-1050
  //
  // ener: energy in eV
  // sig: cross-section in cm^2
  //
  // Version 1.0
  // J. Wilms, wilms@astro.uni-tuebingen.de, 1999/04/22
  //
  // Version 2.0
  // J. Wilms, joern.wilms@sternwarte.uni-erlangen.de, 2016/03/23
  // Conversion to C
  int i;
  
  // Fitting coefficients (Yan et al., Tab. 4)
  static double a[]={1.0,-4.7416,14.8200,-30.8678,37.3584, -23.4585,5.9133};
  // ... parameters Q for resonances (Fernley et al. 1987)
  static double q[]={2.81, 2.51, 2.45, 2.44};

  // ... parameters NU for resonances (Oza 1986)
  static double nu[]={1.610, 2.795, 3.817, 4.824};

  // ... parameters GAMMA for resonances (Oza 1986)
  static double gamma[]={2.64061E-3, 6.20116E-4, 2.56061E-4,1.320159E-4};

  // Threshold
  if (ener<24.58) {
    return(0.);
  }

  // 1./x**(1/2) in Yan notation
  double sq=sqrt(24.58/ener);

  // equation 14
  double sig=horner(sq,a,7);
  sig=sig*733e-6*pow(ener/1000.,-3.5);

  // Add in autoionization resonances
  for (i=0; i<4; i++) {
    sig=sig*jwfano(q[i],nu[i],gamma[i],ener);
  }

  return(sig*MBARN2CM2);
}

double jwfano(double q,double nu,double gamma,double ener) {
  //
  //
  // Source :
  //   Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
  //      J. Phys. B., 20, 6457.
  //   Oza, D. H., (1986), Phys. Rev. A, 33,  824.
  //
  // Description :
  // returns a Fano line profile for use by HELIUM. The form of
  // this line shape is taken from Fernley, Taylor and Seaton (above).
  //
  // Deficiencies :
  //   works in the energy range from 30 eV to 10,000 eV
  // Bugs :
  //  if any are found please report to the authors
  //
  // History :
  //   04 jan 93 : original, based on Jelinsky's program written in
  //               C (EUVE Archive) - (19775::MBC)
  //
  //   23 feb 93 : comments added and modified to remove VAX
  //               fortran 77 extensions (19775::MC)
  //
  //   21 sep 93 : further remnants of VAX fortran 77 extensions
  //               removed (19775::MBC)
  //
  //   22 sep 93 : bug fixed in calculations of EPSI - divided by 1.807317
  //               not added (19775::MBC)
  //
  //   16 jul 2000: modified by Joern Wilms (wilms@astro.uni-tuebingen.de)
  //                fourth argument is now energy in eV
  //                renamed to jwfano
  //   23 Mar 2016: Conversion to C by 
  //                joern.wilms@sternwarte.uni-erlangen.de


  // Energy in rydbergs
  double eps=ener/13.598;
  double epsi=3.0 - 1./(nu*nu) + 1.807317;
  double x=2.0*(eps - epsi)/gamma;
  return((x-q)*(x-q)/(1.0+x*x));
}

double horner(double x,double *coeff,int n) {
  //
  // evaluate the polynomial
  // horner(x)=sum_{i=0}^(n-1) coeff[i] x^i
  // using a Horner scheme
  //
  // for speed reasons this helper function does NOT do
  // any error checking!
  //
  int i;
  double horner=coeff[n-1];
  for (i=n-2;i>=0;i--) {
    horner=horner*x+coeff[i];
  }
  return(horner);
}

double h2abs(double e) {
  //   
  // Cross section for the H2 molecule according to
  // Yan et al., 1998, ApJ 496, 1044-1050
  //
  // e: energy in eV
  // sig: cross-section in cm^2
  //
  // Version 1.0, 1999/05/05 Joern Wilms,
  // wilms@astro.uni-tuebingen.de
  //
  // Version 1.1, 1999/07/06 Joern Wilms
  // Updated to cross section of Wilms, Allen, McCray,1999, ApJ 542, 914-924
  // for energies between 30 and 85eV
  //
  // Version 1.2, 2016/03/23 Joern Wilms
  // Conversion to C, optimization of speed.
  // This version also fixes an error in the Fortran implementation of this
  // routine, where the cross section was off by 8 orders of magnitude between
  // 30 and 80eV. This was not noticed for 17 years...
  //
  // An ifdef allows to use the polynomial coefficients of
  // Yan et al., ApJ 559, 1194, i.e., the erratum to the original article.
  // Unfortunately, these partly result in negative cross sections
  // (e.g., Liu & Shemansky, 2004, ApJ 614:1132), since for the 15.4-18eV interval
  // the 2nd coefficient should not be -197.448 as listed in the erratum,
  // but -197.25 (Aleman & Gruenwald, 2011, A&A 528 A74)
  //
  //

#ifdef YANERRATUM  
  // yan 15.4-18.0 eV
  static double coeff0[]={1,-197.25,+438.823,-260.481,+17.915};
  // yan 18.0-30.0 eV
  static double coeff1[]={-145.528,+351.394,-274.294,+74.320};
  // yan 30.0-85.0 eV
  static double coeff2[]={65.304,-91.762,+51.778,-9.364};
#endif
  // yanold 15.4-18.0eV
  static double oldcoeff0[]={-37.895,+99.723,-87.227,+25.4};
  // yan old 18.0-30.0eV
  static double oldcoeff1[]={0.071,-0.673,+1.977,-0.692};
  // wilms 30.-85eV
  static double wilmscoeff[]={0.664,-11.768,+78.118,-231.339,+368.053,-189.953};
  // yan >85eV
  static double coeff3[]={1.,-2.003,-4.806,+50.577,-171.044,+231.608,-81.885};

  if (e<=15.4) {
    return 0.;
  }
  double x=e/15.4;
  double xsq=1./sqrt(x);
  double sig;

#ifdef YANERRATUM
  if (e<18.) {
    sig=1E7*horner(xsq,coeff0,5);
  } else {
    double e35=pow(e/1000.,-3.5);
    if (e<30.) {
      sig=e35*horner(xsq,coeff1,4);
    } else {
      if (e<85.) {
	sig=e35*horner(xsq,coeff2,4);
      } else {
	sig=45.57*e35*horner(xsq,coeff3,7);
      }
    }
  }
#else
    if (e<18.) {
      sig=1E8*horner(x,oldcoeff0,4);
  } else {
    if (e<30.) {
      sig=2e7*pow(x,-0.252)*horner(1./x,oldcoeff1,4);
    } else {
      if (e<85.) {
	sig=1e6*horner(1./x,wilmscoeff,6);
      } else {
	double e35=pow(e/1000.,-3.5);
	sig=45.57*e35*horner(xsq,coeff3,7);
      }
    }
  }
#endif
  
  
  // barn -> MBarn -> cm2
  sig*=1E-6*MBARN2CM2;

  return(sig);
  
}

double dgami(double x,double y) {
  return dgami_(&x,&y);
}

double fgabnd(const char *ele) {
  return (double) (FGABND((char *) ele));
}


int binsrch(double e,double *ener,int n) {
  //
  // search for energy e in array ener(0:n-1) and return bin k for
  //  which ener[k]<=e<ener[k+1]
  int klo=0;
  int khi=n-1;
  while ((khi-klo)>1) {
    int k=((unsigned int) klo+(unsigned int) khi)/2;
    if(ener[k]>e) {
      khi=k;
    } else {
      klo=k;
    }
  }

  assert(klo>=0 && klo<n);
  assert(ener[klo]<=e && e<ener[klo+1]);

  return(klo);

}


void tbvabsfeo(const double *ear, int ne, const double *param, int ifl, double *photar,
	    double *photer, const char *init) {

  //
  // Wrapper around the ISM absorption code, allows to fit NH and the O and Fe abundances
  // All other parameters are left at their default values
  //
  int i;
  double newparam[42];

  // abundances
  for (i=1; i<=17;i++) {
    newparam[i]=1.;
  }
  // molecular H2
  newparam[18]=0.2;
  // dust
  newparam[19]=1.0; // g/cm^2
  newparam[20]=0.025; // amin, mum
  newparam[21]=0.25; // amax, mum
  newparam[22]=3.5; // PL index

  // depletions
  newparam[23]=1.0; // H
  newparam[24]=1.  ;// He
  newparam[25]=0.5 ;// C
  newparam[26]=1.  ;// N
  newparam[27]=0.6 ;// O
  newparam[28]=1.  ;// Ne
  newparam[29]=0.25;// Na
  newparam[30]=0.2 ;// Mg
  newparam[31]=0.02;// Al
  newparam[32]=0.1 ;// Si
  newparam[33]=0.6 ;// S
  newparam[34]=0.5 ;// Cl
  newparam[35]=1.  ;// Ar
  newparam[36]=0.003;// Ca
  newparam[37]=0.03;// Cr
  newparam[38]=0.3 ;// Fe
  newparam[39]=0.05;// Co
  newparam[40]=0.04;// Ni

  newparam[41]=param[3]; // redshift
  

  newparam[0]=param[0];  // NH
  newparam[4]=param[1];  // O
  newparam[15]=param[2]; // Fe
  
  tbnew(ear,ne,newparam,ifl,photar,photer,init);
}


void tbgas(const double *ear, int ne, const double *param, int ifl, double *photar,
	   double *photer, const char *init) {
  //
  // Wrapper around the ISM absorption code, assuming pure gas absorption
  // (no H2, no dust)
  //
  int i;
  double newparam[42];

  // abundances
  for (i=1; i<=17;i++) {
    newparam[i]=1.;
  }
  // molecular H2
  newparam[18]=0.;
  // dust parameters: rho=0 switches the dust phase off,
  // but just in case we initialize all depletions to 1 as well
  newparam[19]=0.; // g/cm^2 -- 
  newparam[20]=0.025; // amin, mum
  newparam[21]=0.25; // amax, mum
  newparam[22]=3.5; // PL index

  // depletions
  for (i=23; i<=40;i++) {
    newparam[i]=1.0;
  }

  newparam[0]=param[0];  // NH
  newparam[41]=param[1]; // redshift
  
  tbnew(ear,ne,newparam,ifl,photar,photer,init);
}


void tbvabspcf(const double *ear, int ne, const double *param, int ifl, double *photar,
	    double *photer, const char *init) {

  //
  // Wrapper around the ISM absorption code, allows to fit NH and the O and Fe abundances
  // All other parameters are left at their default values
  //
  int i;
  double newparam[42];

  // abundances
  newparam[0]=param[0]; // NH
  for (i=1; i<=17;i++) {
    newparam[i]=1.;
  }
  // molecular H2
  newparam[18]=0.2;
  // dust
  newparam[19]=1.0; // g/cm^2
  newparam[20]=0.025; // amin, mum
  newparam[21]=0.25; // amax, mum
  newparam[22]=3.5; // PL index

  // depletions
  newparam[23]=1.0; // H
  newparam[24]=1.  ;// He
  newparam[25]=0.5 ;// C
  newparam[26]=1.  ;// N
  newparam[27]=0.6 ;// O
  newparam[28]=1.  ;// Ne
  newparam[29]=0.25;// Na
  newparam[30]=0.2 ;// Mg
  newparam[31]=0.02;// Al
  newparam[32]=0.1 ;// Si
  newparam[33]=0.6 ;// S
  newparam[34]=0.5 ;// Cl
  newparam[35]=1.  ;// Ar
  newparam[36]=0.003;// Ca
  newparam[37]=0.03;// Cr
  newparam[38]=0.3 ;// Fe
  newparam[39]=0.05;// Co
  newparam[40]=0.04;// Ni

  newparam[41]=param[2] ; // redshift
  
  tbnew(ear,ne,newparam,ifl,photar,photer,init);

  for (i=0;i<ne;i++) {
    photar[i]=1.+param[1]*(photar[i]-1.0);
  }
  
}

void tbvabsrel(const double *ear, int ne, const double *param, int ifl, double *photar,
	    double *photer, const char *init) {

  //
  // Wrapper around the ISM absorption code, allowing a relative column
  //
  int i;
  double newparam[42];

  // abundances
  int signh=(param[0]<0)? -1: +1;
  
  newparam[0]=fabs(param[0]); // NH
  for (i=1; i<=41;i++) {
    newparam[i]=param[i];
  }
  tbnew(ear,ne,newparam,ifl,photar,photer,init);

  if (signh<0) {
    for (i=0;i<ne;i++) {
      photar[i]=1./photar[i];
    }
  }
  
}

