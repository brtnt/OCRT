/*
 * snf_mie_gen.c -- standalone SnF/OPAC .inp -> OCRT high-resolution .mie generator
 *
 * Reconstructed from the OCRT V1 --mie-gen implementation (Bohren-Huffman/6SV
 * EXSCPHASE lineage). Generates the OCRT canonical 20-wavelength P11/P12/P33
 * cache on a user-selected uniform angular grid (default 0.5 degree, 361 points).
 *
 * Build: gcc -std=c11 -O3 -fopenmp snf_mie_gen.c -lm -o snf_mie_gen
 * Usage: snf_mie_gen INPUT.inp OUTPUT.mie [dtheta_deg] [rlogpas]
 */
#define _USE_MATH_DEFINES
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#define MAX_COMPONENTS 10
#define NWL 20
#define LINE_SIZE 16384
#define MIE_CONVERGENCE_TEST 1.0e-14

static const double WL[NWL] = {
    0.350, 0.400, 0.412, 0.443, 0.470, 0.488, 0.515, 0.550,
    0.590, 0.633, 0.670, 0.694, 0.760, 0.860, 1.240, 1.536,
    1.650, 1.950, 2.250, 3.750
};
static const char *COLNAMES[6] = {
    "Nor_Ext_Co", "Nor_Sca_Co", "Sg_Sca_Alb", "Asymm_Para",
    "Extinct_Co", "Scatter_Co"
};

typedef struct {
    double rmean, sigma, frac;
    double nr[NWL], ni[NWL];
} component_t;

typedef struct {
    double rmin, rmax;
    int ncomp;
    component_t c[MAX_COMPONENTS];
} input_t;

typedef struct {
    int nang;
    double *angles;
    double spectral[NWL][6];
    double *P11, *P12, *P33; /* [nang][NWL] */
} mie_out_t;

static int starts_numeric(const char *s) {
    while (*s && isspace((unsigned char)*s)) ++s;
    return *s=='+' || *s=='-' || *s=='.' || isdigit((unsigned char)*s);
}

static int parse_values(const char *s, double *out, int maxn) {
    int n=0;
    while (*s && n<maxn) {
        while (*s && isspace((unsigned char)*s)) ++s;
        if (!*s) break;
        errno=0;
        char *end=NULL;
        double v=strtod(s,&end);
        if (end==s || errno==ERANGE) break;
        out[n++]=v;
        s=end;
    }
    return n;
}

static int next_numeric_line(FILE *f, char *buf, size_t n) {
    while (fgets(buf,(int)n,f)) {
        if (starts_numeric(buf)) return 1;
    }
    return 0;
}

static int read_input(const char *path, input_t *in) {
    memset(in,0,sizeof(*in));
    FILE *f=fopen(path,"r");
    if (!f) { perror(path); return -1; }
    char buf[LINE_SIZE]; double v[NWL];
    if (!next_numeric_line(f,buf,sizeof(buf)) || parse_values(buf,v,3)<3) {
        fprintf(stderr,"ERROR: invalid header in %s\n",path); fclose(f); return -1;
    }
    in->rmin=v[0]; in->rmax=v[1]; in->ncomp=(int)llround(v[2]);
    if (!(in->rmin>0.0) || !(in->rmax>in->rmin) || in->ncomp<1 || in->ncomp>MAX_COMPONENTS) {
        fprintf(stderr,"ERROR: invalid r range/component count in %s\n",path); fclose(f); return -1;
    }
    for (int i=0;i<in->ncomp;i++) {
        if (!next_numeric_line(f,buf,sizeof(buf)) || parse_values(buf,v,3)<3) {
            fprintf(stderr,"ERROR: component %d header in %s\n",i+1,path); fclose(f); return -1;
        }
        in->c[i].rmean=v[0]; in->c[i].sigma=v[1]; in->c[i].frac=v[2];
        if (!next_numeric_line(f,buf,sizeof(buf)) || parse_values(buf,in->c[i].nr,NWL)<NWL) {
            fprintf(stderr,"ERROR: component %d nr in %s\n",i+1,path); fclose(f); return -1;
        }
        if (!next_numeric_line(f,buf,sizeof(buf)) || parse_values(buf,in->c[i].ni,NWL)<NWL) {
            fprintf(stderr,"ERROR: component %d ni in %s\n",i+1,path); fclose(f); return -1;
        }
    }
    fclose(f); return 0;
}

static double lognormal(double r,double rm,double sigma) {
    double ls=log10(sigma), a=log10(r/rm);
    double den=sqrt(2.0*M_PI)*log(10.0)*r*ls;
    return exp(-a*a/(2.0*ls*ls))/den;
}

static int mie_series_order(double Y) {
    int N=(int)(0.5*(-1.0+sqrt(1.0+4.0*Y*Y)))+1;
    if (N==1) N=2;
    int Np=N;
    double Up=2.0*Y/(2.0*Np+1.0);
    double den=1.0-Up;
    if (fabs(den)<1e-14) den=(den<0?-1:1)*1e-14;
    int mu1=(int)(Np+30.0*(0.10+0.35*Up*(2.0-Up*Up)/2.0/den));
    int mu2=1000000;
    int Np_alt=(int)(Y-0.5+sqrt(30.0*0.35*Y));
    if (Np_alt>N) {
        double Up2=2.0*Y/(2.0*Np_alt+1.0);
        double den2=1.0-Up2;
        if (fabs(den2)<1e-14) den2=(den2<0?-1:1)*1e-14;
        mu2=(int)(Np_alt+30.0*(0.10+0.35*Up2*(2.0-Up2*Up2)/2.0/den2));
    }
    int mu=mu1<mu2?mu1:mu2;
    if (mu<2) mu=2;
    return mu;
}

/* Exact OCRT V1/6SV-lineage single-sphere computation. */
static int exscphase(double x,double nr,double ni,const double *cosmu,int nang,
                     double *Qext,double *Qsca,double *p11,double *q11,double *u11) {
    double Ren=nr/(nr*nr+ni*ni), Imn=ni/(nr*nr+ni*ni);
    double Y=x*sqrt(nr*nr+ni*ni);
    int mu=mie_series_order(Y);
    size_t N=(size_t)mu+2;
    double *mem=(double*)calloc(N*14,sizeof(double));
    if (!mem) return -1;
    double *xj=mem+N*0, *xy=mem+N*1, *Rn=mem+N*2;
    double *RDnY=mem+N*3, *IDnY=mem+N*4, *RDnX=mem+N*5;
    double *RGnX=mem+N*6, *IGnX=mem+N*7;
    double *RAn=mem+N*8, *IAn=mem+N*9, *RBn=mem+N*10, *IBn=mem+N*11;
    double *PIn=mem+N*12, *TAUn=mem+N*13;

    Rn[mu]=0.0;
    int k=mu+1,mub=0;
    xj[mu+1]=0.0;
    while (1) {
        --k; xj[k]=0.0;
        xj[k]=0.0;
        Rn[k-1]=x/(2.0*k+1.0-x*Rn[k]);
        if (k==2) { mub=mu; xj[mub+1]=0.0; xj[mub]=1.0; break; }
        if (Rn[k-1]>1.0) { mub=k-1; xj[mub+1]=Rn[mub]; xj[mub]=1.0; break; }
    }
    for (int kk=mub;kk>=1;--kk) xj[kk-1]=(2.0*kk+1.0)*xj[kk]/x-xj[kk+1];
    double coxj=(xj[0]-x*xj[1])*cos(x)+x*xj[0]*sin(x);

    RDnY[mu]=IDnY[mu]=RDnX[mu]=0.0;
    for (int kk=mu;kk>=1;--kk) {
        RDnX[kk-1]=kk/x-1.0/(RDnX[kk]+kk/x);
        double xr=RDnY[kk]+Ren*kk/x, xi=IDnY[kk]+Imn*kk/x;
        double xd=xr*xr+xi*xi;
        RDnY[kk-1]=kk*Ren/x-xr/xd;
        IDnY[kk-1]=kk*Imn/x+xi/xd;
    }

    double xy_m1=sin(x)/x;
    xy[0]=-cos(x)/x;
    RGnX[0]=0.0; IGnX[0]=-1.0;
    *Qext=0.0; *Qsca=0.0;
    int muf=mu;
    for (int kk=1;kk<=mu;++kk) {
        if (kk<=mub) xj[kk]/=coxj; else xj[kk]=Rn[kk-1]*xj[kk-1];
        if (kk==1) xy[kk]=(2.0*kk-1.0)*xy[kk-1]/x-xy_m1;
        else       xy[kk]=(2.0*kk-1.0)*xy[kk-1]/x-xy[kk-2];
        double xJonH=xj[kk]/(xj[kk]*xj[kk]+xy[kk]*xy[kk]);
        double gd=(RGnX[kk-1]-kk/x)*(RGnX[kk-1]-kk/x)+IGnX[kk-1]*IGnX[kk-1];
        RGnX[kk]=(kk/x-RGnX[kk-1])/gd-kk/x;
        IGnX[kk]=IGnX[kk-1]/gd;

        double a1=RDnY[kk]-nr*RDnX[kk];
        double a2=IDnY[kk]+ni*RDnX[kk];
        double ad1=RDnY[kk]-nr*RGnX[kk]-ni*IGnX[kk];
        double ad2=IDnY[kk]+ni*RGnX[kk]-nr*IGnX[kk];
        double ad=ad1*ad1+ad2*ad2;
        double arb=(a1*ad1+a2*ad2)/ad;
        double aib=(-a1*ad2+a2*ad1)/ad;
        RAn[kk]=xJonH*(xj[kk]*arb-xy[kk]*aib);
        IAn[kk]=xJonH*(xy[kk]*arb+xj[kk]*aib);

        double b1=nr*RDnY[kk]+ni*IDnY[kk]-RDnX[kk];
        double b2=nr*IDnY[kk]-ni*RDnY[kk];
        double bd1=nr*RDnY[kk]+ni*IDnY[kk]-RGnX[kk];
        double bd2=nr*IDnY[kk]-ni*RDnY[kk]-IGnX[kk];
        double bd=bd1*bd1+bd2*bd2;
        double brb=(b1*bd1+b2*bd2)/bd;
        double bib=(-b1*bd2+b2*bd1)/bd;
        RBn[kk]=xJonH*(xj[kk]*brb-xy[kk]*bib);
        IBn[kk]=xJonH*(xy[kk]*brb+xj[kk]*bib);

        double temp=RAn[kk]*RAn[kk]+IAn[kk]*IAn[kk]+RBn[kk]*RBn[kk]+IBn[kk]*IBn[kk];
        if (temp/kk<MIE_CONVERGENCE_TEST) { muf=kk; break; }
        double xp=2.0/(x*x)*(2.0*kk+1.0);
        *Qsca += xp*temp;
        *Qext += xp*(RAn[kk]+RBn[kk]);
    }

    for (int j=0;j<nang;++j) {
        double muang=cosmu[j];
        double RS1=0,RS2=0,IS1=0,IS2=0;
        PIn[0]=0.0; PIn[1]=1.0; TAUn[1]=muang;
        for (int kk=1;kk<=muf;++kk) {
            double cn=(2.0*kk+1.0)/(kk*(kk+1.0));
            RS1 += cn*(RAn[kk]*PIn[kk]+RBn[kk]*TAUn[kk]);
            RS2 += cn*(RAn[kk]*TAUn[kk]+RBn[kk]*PIn[kk]);
            IS1 += cn*(IAn[kk]*PIn[kk]+IBn[kk]*TAUn[kk]);
            IS2 += cn*(IAn[kk]*TAUn[kk]+IBn[kk]*PIn[kk]);
            PIn[kk+1]=((2.0*kk+1.0)*muang*PIn[kk]-(kk+1.0)*PIn[kk-1])/kk;
            TAUn[kk+1]=(kk+1.0)*muang*PIn[kk+1]-(kk+2.0)*PIn[kk];
        }
        p11[j]=(RS1*RS1+IS1*IS1+RS2*RS2+IS2*IS2)/(x*x)/2.0;
        q11[j]=(RS2*RS2+IS2*IS2-RS1*RS1-IS1*IS1)/(x*x)/2.0;
        u11[j]=(2.0*RS2*RS1+2.0*IS2*IS1)/(x*x)/2.0;
    }
    free(mem); return 0;
}

static int alloc_output(mie_out_t *o,double dtheta) {
    memset(o,0,sizeof(*o));
    double q=180.0/dtheta;
    long nstep=lround(q);
    if (fabs(q-nstep)>1e-10 || nstep<1) return -1;
    o->nang=(int)nstep+1;
    o->angles=(double*)malloc((size_t)o->nang*sizeof(double));
    size_t n=(size_t)o->nang*NWL;
    o->P11=(double*)calloc(n,sizeof(double));
    o->P12=(double*)calloc(n,sizeof(double));
    o->P33=(double*)calloc(n,sizeof(double));
    if (!o->angles||!o->P11||!o->P12||!o->P33) return -1;
    for (int k=0;k<o->nang;k++) o->angles[k]=180.0-k*dtheta;
    o->angles[o->nang-1]=0.0;
    return 0;
}
static void free_output(mie_out_t *o) { free(o->angles);free(o->P11);free(o->P12);free(o->P33); }
#define IDX(k,l) ((size_t)(k)*NWL+(l))

static double trap(const double *x,const double *y,int n) {
    double s=0; for(int i=1;i<n;i++) s+=0.5*(x[i]-x[i-1])*(y[i]+y[i-1]); return s;
}

static int compute(const input_t *in,mie_out_t *o,double rlogpas) {
    int nc=in->ncomp, na=o->nang;
    double *cosmu=(double*)malloc((size_t)na*sizeof(double));
    double *ext=(double*)calloc((size_t)nc*NWL,sizeof(double));
    double *sca=(double*)calloc((size_t)nc*NWL,sizeof(double));
    size_t psz=(size_t)nc*na*NWL;
    double *pc=(double*)calloc(psz,sizeof(double)),*qc=(double*)calloc(psz,sizeof(double)),*uc=(double*)calloc(psz,sizeof(double));
    double *vi=(double*)calloc((size_t)nc,sizeof(double));
    if(!cosmu||!ext||!sca||!pc||!qc||!uc||!vi) return -1;
    for(int k=0;k<na;k++) cosmu[k]=cos(o->angles[k]*M_PI/180.0);
    double sm=pow(10.0,rlogpas)-1.0;
    /* Volume integral exactly on the same left-point geometric radius grid. */
    for(int i=0;i<nc;i++) {
        for(double r=in->rmin;r<in->rmax;) {
            double dr=r*sm, dn=lognormal(r,in->c[i].rmean,in->c[i].sigma);
            vi[i]+=r*r*r*dn*dr; r+=dr;
        }
    }

    int failed=0;
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<nc;i++) for(int l=0;l<NWL;l++) {
        double *pt=(double*)malloc((size_t)na*sizeof(double));
        double *qt=(double*)malloc((size_t)na*sizeof(double));
        double *ut=(double*)malloc((size_t)na*sizeof(double));
        if(!pt||!qt||!ut) {
            #pragma omp atomic write
            failed=1;
        } else {
            for(double r=in->rmin;r<in->rmax;) {
                double dr=r*sm, dn=lognormal(r,in->c[i].rmean,in->c[i].sigma);
                double xnd=dn*dr*M_PI*r*r;
                double x=2.0*M_PI*r/WL[l],qe=0,qs=0;
                if(exscphase(x,in->c[i].nr[l],in->c[i].ni[l],cosmu,na,&qe,&qs,pt,qt,ut)!=0) {
                    #pragma omp atomic write
                    failed=1; break;
                }
                ext[i*NWL+l]+=xnd*qe; sca[i*NWL+l]+=xnd*qs;
                for(int k=0;k<na;k++) {
                    size_t z=((size_t)i*na+k)*NWL+l;
                    pc[z]+=4.0*pt[k]*xnd; qc[z]+=4.0*qt[k]*xnd; uc[z]+=4.0*ut[k]*xnd;
                }
                r+=dr;
            }
        }
        free(pt);free(qt);free(ut);
    }
    if(failed) { free(cosmu);free(ext);free(sca);free(pc);free(qc);free(uc);free(vi);return -1; }

    double denom=0,cij[MAX_COMPONENTS];
    for(int i=0;i<nc;i++) { vi[i]*=4.0*M_PI/3.0; denom+=in->c[i].frac/vi[i]; }
    for(int i=0;i<nc;i++) cij[i]=(in->c[i].frac/vi[i])/denom;
    for(int l=0;l<NWL;l++) {
        double ee=0,ss=0;
        for(int i=0;i<nc;i++) { ee+=cij[i]*ext[i*NWL+l]; ss+=cij[i]*sca[i*NWL+l]; }
        o->spectral[l][4]=ee; o->spectral[l][5]=ss;
    }
    for(int k=0;k<na;k++) for(int l=0;l<NWL;l++) {
        double p=0,q=0,u=0;
        for(int i=0;i<nc;i++) { size_t z=((size_t)i*na+k)*NWL+l; p+=cij[i]*pc[z];q+=cij[i]*qc[z];u+=cij[i]*uc[z]; }
        double ss=o->spectral[l][5]; o->P11[IDX(k,l)]=p/ss;o->P12[IDX(k,l)]=q/ss;o->P33[IDX(k,l)]=u/ss;
    }
    double ext550=o->spectral[7][4],sca550=o->spectral[7][5];
    double *mu=(double*)malloc((size_t)na*sizeof(double)),*pv=(double*)malloc((size_t)na*sizeof(double)),*pmu=(double*)malloc((size_t)na*sizeof(double));
    if(!mu||!pv||!pmu) return -1;
    /* Input angles descend, so cos(theta) already ascends from -1 to +1. */
    for(int k=0;k<na;k++) mu[k]=cosmu[k];
    for(int l=0;l<NWL;l++) {
        o->spectral[l][0]=o->spectral[l][4]/ext550;
        o->spectral[l][1]=o->spectral[l][5]/sca550;
        o->spectral[l][2]=o->spectral[l][5]/o->spectral[l][4];
        for(int k=0;k<na;k++) { pv[k]=o->P11[IDX(k,l)];pmu[k]=pv[k]*mu[k]; }
        o->spectral[l][3]=trap(mu,pmu,na)/trap(mu,pv,na);
    }
    free(mu);free(pv);free(pmu);free(cosmu);free(ext);free(sca);free(pc);free(qc);free(uc);free(vi);return 0;
}

static int write_mie(const mie_out_t *o,const char *path) {
    FILE *f=fopen(path,"w"); if(!f){perror(path);return -1;}
    fprintf(f," %d\n",o->nang);
    fprintf(f,"   Wlgth"); for(int c=0;c<6;c++) fprintf(f,"  %10s",COLNAMES[c]); fprintf(f,"\n");
    for(int l=0;l<NWL;l++) fprintf(f,"    %.4f     %.4f        %.4f        %.4f        %.4f     %.4E  %.4E\n",
        WL[l],o->spectral[l][0],o->spectral[l][1],o->spectral[l][2],o->spectral[l][3],o->spectral[l][4],o->spectral[l][5]);
    fprintf(f,"\n\n");
    const char *labels[3]={"Phase Function (P11)","Phase Function (P12 / Q-polarization)","Phase Function (P33 / U-polarization)"};
    const double *b[3]={o->P11,o->P12,o->P33};
    for(int ib=0;ib<3;ib++) {
        fprintf(f,"                    %s \n",labels[ib]);
        fprintf(f,"   TETA");for(int l=0;l<NWL;l++)fprintf(f,"  %8.4f",WL[l]);fprintf(f,"\n");
        for(int k=0;k<o->nang;k++) { fprintf(f," %7.2f",o->angles[k]);for(int l=0;l<NWL;l++)fprintf(f,"  %+.4E",b[ib][IDX(k,l)]);fprintf(f,"\n"); }
        fprintf(f,"\n");
    }
    fclose(f);return 0;
}

int main(int argc,char **argv) {
    if(argc<3||argc>5){fprintf(stderr,"Usage: %s INPUT.inp OUTPUT.mie [dtheta_deg=0.5] [rlogpas=0.011]\n",argv[0]);return 2;}
    double dtheta=argc>=4?strtod(argv[3],NULL):0.5;
    double rlogpas=argc>=5?strtod(argv[4],NULL):0.011;
    if(!(dtheta>0&&dtheta<=180&&rlogpas>0&&rlogpas<1)){fprintf(stderr,"ERROR: invalid dtheta/rlogpas\n");return 2;}
    input_t in; mie_out_t out;
    if(read_input(argv[1],&in)!=0||alloc_output(&out,dtheta)!=0)return 1;
    fprintf(stderr,"[mie-gen] %s: components=%d r=[%.6g,%.6g] dtheta=%.6g nang=%d rlogpas=%.6g threads=%d\n",
        argv[1],in.ncomp,in.rmin,in.rmax,dtheta,out.nang,rlogpas,
#ifdef _OPENMP
        omp_get_max_threads()
#else
        1
#endif
    );
    if(compute(&in,&out,rlogpas)!=0){fprintf(stderr,"ERROR: Mie computation failed\n");free_output(&out);return 1;}
    if(write_mie(&out,argv[2])!=0){free_output(&out);return 1;}
    fprintf(stderr,"[mie-gen] wrote %s\n",argv[2]);
    free_output(&out);return 0;
}
