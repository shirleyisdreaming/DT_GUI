#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>
#include <nlopt.h>
#include <stdarg.h>
#include <fcntl.h>
#include <ctype.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <stdint.h>
#include <tiffio.h>
#include <libgen.h>


#define nFitVals 12
#define NrValsFitOutput 12
struct my_profile_func_data{
    int NrPtsForFit;
    double *Rs;
    double *PeakShape;
};

static
double problem_function_profile(unsigned n,const double *x,double *grad,void* f_data_trial,double *Rstarts, int nRadFits)
{
    
    struct my_profile_func_data *f_data = (struct my_profile_func_data *) f_data_trial;
    int NrPtsForFit = f_data->NrPtsForFit;
    double *Rs, *PeakShape;
    Rs = &(f_data->Rs[0]);
    PeakShape = &(f_data->PeakShape[0]);
    double Rcen, Mu, SigmaG, SigmaL, Imax, BG;
    Rcen = x[0];
    Mu = x[1];
    SigmaG = x[2];
    SigmaL = x[2];
    Imax = x[3];
    BG = x[4];
    double TotalDifferenceIntensity=0,CalcIntensity;
    int i,j,k;
    double L, G;
    for (i=0;i<NrPtsForFit;i++){
        L = (1/(((Rs[i]-Rcen)*(Rs[i]-Rcen)/(SigmaL*SigmaL))+(1)));
        G = (exp((-0.5)*(Rs[i]-Rcen)*(Rs[i]-Rcen)/(SigmaG*SigmaG)));
        CalcIntensity = BG + Imax*((Mu*L)+((1-Mu)*G));
        TotalDifferenceIntensity += (CalcIntensity - PeakShape[i])*(CalcIntensity - PeakShape[i]);
    }
#ifdef PRINTOPT
    printf("Peak profiler intensity difference: %f\n",TotalDifferenceIntensity);
#endif
    return TotalDifferenceIntensity;
}

static
double CalcIntegratedIntensity(
    const double *x,
    void* f_data_trial)
{
    struct my_profile_func_data *f_data = (struct my_profile_func_data *) f_data_trial;
    int NrPtsForFit = f_data->NrPtsForFit;
    double *Rs;
    Rs = &(f_data->Rs[0]);
    double Rcen, Mu, SigmaG, SigmaL, Imax, BG;
    Rcen = x[0];
    Mu = x[1];
    SigmaG = x[2];
    SigmaL = x[2];
    Imax = x[3];
    BG = x[4];
    double TotalIntensity=0;
    int i,j,k;
    double L, G;
    for (i=0;i<NrPtsForFit;i++){
        L = (1/(((Rs[i]-Rcen)*(Rs[i]-Rcen)/(SigmaL*SigmaL))+(1)));
        G = (exp((-0.5)*(Rs[i]-Rcen)*(Rs[i]-Rcen)/(SigmaG*SigmaG)));
        TotalIntensity += Imax*((Mu*L)+((1-Mu)*G));
    }
#ifdef PRINTOPT2
    printf("Peak fit intensity value: %lf\n",TotalIntensity);
#endif
    return TotalIntensity;
}

void FitPeakShape(int NrPtsForFit, double Rs[NrPtsForFit], double PeakShape[NrPtsForFit],
                double *Rfit)
{
   
    unsigned n = 5;
    double x[n],xl[n],xu[n];
    struct my_profile_func_data f_data;
    f_data.NrPtsForFit = NrPtsForFit;
    f_data.Rs = &Rs[0];
    f_data.PeakShape = &PeakShape[0];
    double BG0 = (PeakShape[0]+PeakShape[NrPtsForFit-1])/2;
    if (BG0 < 0) BG0=0;
    double MaxI=-100000, TotInt = 0;
    double Rmean = (Rs[0] + Rs[NrPtsForFit-1])/2;
    double RTemp;
    int i;
    for (i=0;i<NrPtsForFit;i++){
        TotInt += PeakShape[i];
        if (PeakShape[i] > MaxI){
            MaxI=PeakShape[i];
            RTemp = Rs[i];
        }
    }
    if (fabs(RTemp-Rmean)<2.0) Rmean = RTemp;
    MaxI -= BG0;
    x[0] = Rmean; xl[0] = Rs[0];    xu[0] = Rs[NrPtsForFit-1];
    x[1] = 0.5;   xl[1] = 0;        xu[1] = 1;
    x[2] = 1;     xl[2] = 0.05;  xu[2] = 100;
    //~ x[3] = 1;     xl[3] = 0.05;  xu[3] = 100;
    x[3] = MaxI;  xl[3] = MaxI/100; xu[3] = MaxI*3;
    x[4] = BG0;   xl[4] = 0;        xu[4] = BG0*3;
    struct my_profile_func_data *f_datat;
    f_datat = &f_data;
    void* trp = (struct my_profile_func_data *) f_datat;
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, n);
    nlopt_set_lower_bounds(opt, xl);
    nlopt_set_upper_bounds(opt, xu);
    nlopt_set_min_objective(opt, problem_function_profile, trp);
    double minf,MeanDiff;
    nlopt_optimize(opt, x, &minf);
    nlopt_destroy(opt);
    MeanDiff = sqrt(minf)/(NrPtsForFit);
    Rfit[0] = x[0];
    Rfit[1] = x[1];
    Rfit[2] = x[2];
    Rfit[3] = x[2];
    Rfit[4] = x[3];
    Rfit[5] = MaxI; // Input Max Intensity
    Rfit[6] = x[4];
    Rfit[7] = BG0;
    Rfit[8] = MeanDiff;
    Rfit[9] = CalcIntegratedIntensity(x,trp); // Calculate integrated intensity
    Rfit[10] = TotInt; // Total intensity
    Rfit[11] = TotInt - (BG0*NrPtsForFit); // Total intensity after removing background
    //~ for (i=0;i<NrValsFitOutput;i++) printf("%lf ",Rfit[i]); printf("\n");
}

struct my_profile_func_data_multipeak{
    int NrPtsForFit;
    double *Rs;
    double *PeakShape;
    int NrCenter;
    double *Rcenter;
};

static
double problem_function_profile_multipeak(
    unsigned n,
    const double *x,
    double *grad,
    void* f_data_trial)
{
    int i,j,k;
    double L, G;
    struct my_profile_func_data_multipeak *f_data = (struct my_profile_func_data_multipeak *) f_data_trial;
    int NrPtsForFit = f_data->NrPtsForFit;
    int NrCenter=f_data->NrCenter;
    double *Rs, *PeakShape,*Rcenter,*CalcIntensity;
    Rs = &(f_data->Rs[0]);
    PeakShape = &(f_data->PeakShape[0]);
    Rcenter=&(f_data->Rcenter[0]);
    CalcIntensity = calloc(NrPtsForFit,sizeof(*CalcIntensity));
    for (i=0;i<NrPtsForFit;i++){CalcIntensity[i]=0;}
    double Rcen, Mu, SigmaG, SigmaL, Imax, BG;
    double TotalDifferenceIntensity=0;
   
    for(i=0;i<NrCenter;i++){
        Rcen = x[0+5*i];
        Mu = x[1+5*i];
        SigmaG = x[2+5*i];
        SigmaL = x[2+5*i];
        Imax = x[3+5*i];
        BG = x[4+5*i];
        for (j=0;j<NrPtsForFit;j++){
            L = (1/(((Rs[j]-Rcen)*(Rs[j]-Rcen)/(SigmaL*SigmaL))+(1)));
            G = (exp((-0.5)*(Rs[j]-Rcen)*(Rs[j]-Rcen)/(SigmaG*SigmaG)));
            CalcIntensity[j] =CalcIntensity[j]+BG + Imax*((Mu*L)+((1-Mu)*G));
        }
    }
    for (i=0;i<NrPtsForFit;i++){
        TotalDifferenceIntensity += (CalcIntensity[i] - PeakShape[i])*(CalcIntensity[i] - PeakShape[i]);
    }
	free(CalcIntensity);
#ifdef PRINTOPT
    printf("Peak profiler intensity difference: %f\n",TotalDifferenceIntensity);
#endif
    return TotalDifferenceIntensity;
}


static
double CalcIntegratedIntensity_multipeak(
    const double *x,
    void* f_data_trial)
{
    struct my_profile_func_data_multipeak *f_data = (struct my_profile_func_data_multipeak *) f_data_trial;
    int NrPtsForFit = f_data->NrPtsForFit;
    int NrCenter=f_data->NrCenter;
    double *Rs, *PeakShape,*Rcenter,*CalcIntensity;
    Rs = &(f_data->Rs[0]);
    PeakShape = &(f_data->PeakShape[0]);
    Rcenter=&(f_data->Rcenter[0]);
    double Rcen, Mu, SigmaG, SigmaL, Imax, BG;
    double TotalIntensity=0;
    CalcIntensity = calloc(NrPtsForFit,sizeof(*CalcIntensity));
    int i,j,k;
    double L, G;
    for (i=0;i<NrPtsForFit;i++){CalcIntensity[i]=0;}
    
   
    for(i=0;i<NrCenter;i++){
        Rcen = x[0+5*i];
        Mu = x[1+5*i];
        SigmaG = x[2+5*i];
        SigmaL = x[2+5*i];
        Imax = x[3+5*i];
        BG = x[4+5*i];
        for (j=0;j<NrPtsForFit;j++){
            L = (1/(((Rs[j]-Rcen)*(Rs[j]-Rcen)/(SigmaL*SigmaL))+(1)));
            G = (exp((-0.5)*(Rs[j]-Rcen)*(Rs[j]-Rcen)/(SigmaG*SigmaG)));
            CalcIntensity[j] =CalcIntensity[j]+BG + Imax*((Mu*L)+((1-Mu)*G));
        }
    }
    for (i=0;i<NrPtsForFit;i++){
        TotalIntensity += CalcIntensity[i];
    }
    free(CalcIntensity);
#ifdef PRINTOPT2
    printf("Peak fit intensity value: %lf\n",TotalIntensity);
#endif
    return TotalIntensity;
}


void FitMultiPeakShape(int NrPtsForFit, double Rs[NrPtsForFit], double PeakShape[NrPtsForFit],
                double *Rfit,int NrCenter,double Rcenter[NrCenter], double Rwidth)
{
    struct my_profile_func_data_multipeak f_data;
    f_data.NrPtsForFit = NrPtsForFit;
    f_data.Rs = &Rs[0];
    f_data.PeakShape = &PeakShape[0];
    f_data.NrCenter=NrCenter;
    f_data.Rcenter=&Rcenter[0];
    double Rmean[NrCenter];
    int PosTemp[NrCenter];
    double RTemp[NrCenter];
    unsigned n = 5*NrCenter;
    double x[n],xl[n],xu[n];
    double MaxI[NrCenter], TotInt[NrCenter];
    memset(MaxI,-10000,NrCenter*sizeof(double));
    memset(TotInt,0,NrCenter*sizeof(double));
    double BG0 = (PeakShape[0]+PeakShape[NrPtsForFit-1])/2;
    if (BG0 < 0) BG0=0;
    int i,j,k;
    for (i=0;i<NrCenter;i++){
        for (j=0;j<(NrPtsForFit-1);j++){
            if((Rcenter[i]>=Rs[j])&&(Rcenter[i]<Rs[j+1])){
                MaxI[i] = PeakShape[j];}}}
    for (i=0;i<NrCenter;i++){
        for (j=0;j<(NrPtsForFit);j++){
            TotInt[i]+=PeakShape[j];}}
    for(i=0;i<NrCenter;i++){
        //~ x[0+5*i] = Rmean[i]; xl[0+5*i] = Rcenter[i];    xu[0+5*i] = Rcenter[i+1];
        x[0+5*i] = Rcenter[i]; xl[0+5*i] = Rcenter[i]-Rwidth;    xu[0+5*i] = Rcenter[i]+Rwidth;
        x[1+5*i] = 0.5;   xl[1+5*i] = 0;        xu[1+5*i] = 1;
        x[2+5*i] = 1;     xl[2+5*i] = 0.05;  xu[2+5*i] = 100;
        //~ x[3] = 1;     xl[3] = 0.05;  xu[3] = 100;
        x[3+5*i] = MaxI[i];
        xl[3+5*i] = MaxI[i]/100;
        xu[3+5*i] = MaxI[i]*10;
        x[4+5*i] = BG0;   xl[4+5*i] = 0;        xu[4+5*i] = BG0*3;
    }
    struct my_profile_func_data_multipeak *f_datat;
    f_datat = &f_data;
    void* trp = (struct my_profile_func_data_multipeak *) f_datat;
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, n);
    nlopt_set_lower_bounds(opt, xl);
    nlopt_set_upper_bounds(opt, xu);
    nlopt_set_min_objective(opt, problem_function_profile_multipeak, trp);
    double minf,MeanDiff;
    nlopt_optimize(opt, x, &minf);
    nlopt_destroy(opt);
    MeanDiff = sqrt(minf)/(NrPtsForFit);
    for(i=0;i<(12*NrCenter);i++){Rfit[i]=0;}
    for(i=0;i<NrCenter;i++)
     {
        Rfit[0+12*i] = x[0+5*i];
         Rfit[1+12*i] = x[1+5*i];
         Rfit[2+12*i] = x[2+5*i];
        Rfit[3+12*i] = x[2+5*i];
         Rfit[4+12*i] = x[3+5*i];
        Rfit[5+12*i] = MaxI[i]; // Input Max Intensity
         Rfit[6+12*i] = x[4+5*i];
        Rfit[7+12*i] = BG0;
         Rfit[8+12*i] = MeanDiff;
        Rfit[9+12*i] = CalcIntegratedIntensity_multipeak(x,trp); // Calculate integrated intensity
         Rfit[10+12*i] = TotInt[i]; // Total intensity
         Rfit[11+12*i] = TotInt[i] - (BG0*NrPtsForFit); // Total intensity after removing background
     }
    
}
    
static inline
int StartsWith(const char *a, const char *b)
{
    if (strncmp(a,b,strlen(b)) == 0) return 1;
    return 0;
}

int main(int argc, char **argv)
{
    if (argc != 2){
        printf("Usage: ./peakfit paramFN\n");
        return 1;
    }
   
    
    char *paramFN = argv[1],FileStem[4096],darkFN[4096],ext[4096];
    int startNr,endNr;
    int pad;
    int numProcs;
    int nFrames;
    int nFiles;
    int nFits;
    int multipeak;
    double Rstarts[200];
    int nRadFits = 0, ReconSize,nEtaFits=0,nElsPerRad=0,nRcenter=0;
    double Rmin;
    int NrCenter=0;
    double Rcenter[300];
    char  peakfittingFN[4096],PeakfitresultFN[4096],MultiPeakfitresultFN[4096];
    FILE *paramFile;
    char aline[4096], dummy[4096], *str,outputFolder[4096],outputFolder_peakfit[4096];
    int separateFolder = 0;
    paramFile = fopen(paramFN,"r");
    double Rstep;
    double Rwidth;
    double radiiToFit[200][6], etasToFit[200][4];
    double *ResultArr,*ResultArrAll;
    
   
    while (fgets(aline,4096,paramFile) != NULL){
        str = "fStem ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, FileStem);
        }
        str = "ext ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, ext);
        }
        str = "startNr ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &startNr);
        }
        str = "endNr ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &endNr);
        }
        str = "pad ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &pad);
        }
        str = "numProcs ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &numProcs);
        }
        str = "nFrames ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &nFrames);
        }
        str = "nFiles ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &nFiles);
        }
        str = "RawDataPeakFN ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, peakfittingFN);
        }
        str = "multipeak ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &multipeak);
        }
        str = "PeakFitResultFN ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, PeakfitresultFN);
        }
        str = "EtaToFit ";
        if (StartsWith(aline,str) == 1){
            nEtaFits++;
        }
        
        str = "RadiusToFit ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &Rstarts[nRadFits]);
            nRadFits++;
        }
        str = "RMin ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &Rmin);
        }
        str = "RBinSize ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &Rstep);
        }
        str = "ReconSize ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &ReconSize);
        }
        str = "nElsPerRad ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &nElsPerRad);
        }
     
        str = "Rcenter ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &Rcenter[NrCenter]);
            NrCenter++;
        }
        str = "fitWidth ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &Rwidth);
        }
        str = "OutFolder_Peakfit ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, outputFolder_peakfit);
            separateFolder = 1;
        }
        str = "OutFolder_Integration ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, outputFolder);
            separateFolder = 1;
        }
        
    }

    if(multipeak==1){Rstarts[0]=Rmin;}
//method 1 or method 2
//method 1 from rawdataforpeakfit calculated from tomo
//method 2 nEta*nRad*12*(nframe*nfile) from lineoutferfile
    
//do peakfit for method 1//
    nFits = (ReconSize*ReconSize)*nEtaFits*nRadFits;
    if(multipeak==1)
    {ResultArrAll = calloc(nFitVals*nFits*NrCenter,sizeof(*ResultArrAll));}
    if(multipeak==0){ResultArrAll = calloc(nFitVals*nFits,sizeof(*ResultArrAll));}
    double *rawData;
    rawData = calloc(nFits*nElsPerRad,sizeof(double));
//printf("%d %d %d %d %d\n",ReconSize,nEtas,nRadFits,nFits,nElsPerRad);
    FILE *peakFile,*fitresultFile;
    peakFile=fopen(peakfittingFN,"rb");
    fread(rawData,nFits*nElsPerRad*sizeof(double),1,peakFile);
    fitresultFile=fopen(PeakfitresultFN,"w");
    int fitNr;
    double *RsAll;
    RsAll = calloc(nElsPerRad*nFits,sizeof(*RsAll));
    
	#pragma omp parallel for num_threads(numProcs) private(fitNr) schedule(dynamic)
    for (fitNr=0;fitNr<nFits;fitNr++){
     
        double *Peakshape;
        Peakshape = &rawData[fitNr*nElsPerRad];
        int i, iRad;
        double *Rs;
        Rs = &RsAll[fitNr*nElsPerRad];
        memset(Rs,0,nElsPerRad*sizeof(double));
		//~ iRad = fitNr % nRadFits;
        for (i=0;i<nElsPerRad;i++){Rs[i]=Rstarts[0]+i*Rstep;}
        if(multipeak==1){ResultArr = &ResultArrAll[fitNr*nFitVals*NrCenter];}
        if(multipeak==0){ResultArr = &ResultArrAll[fitNr*nFitVals];}
        if(multipeak==0){FitPeakShape(nElsPerRad,Rs,Peakshape,ResultArr);}
        if(multipeak==1){FitMultiPeakShape(nElsPerRad,Rs,Peakshape,ResultArr,NrCenter,Rcenter,Rwidth);}
        printf("%d %d\n",fitNr,nFits);
        //fwrite(ResultArr,nFitVals*NrCenter*sizeof(double),1,fitresultFile);
    }
    
  
if(multipeak==1)
{fwrite(ResultArrAll,NrCenter*nFitVals*nFits*sizeof(double),1,fitresultFile);}
if(multipeak==0)
{fwrite(ResultArrAll,nFitVals*nFits*sizeof(double),1,fitresultFile);}
fclose(fitresultFile);
    
   
 
}
