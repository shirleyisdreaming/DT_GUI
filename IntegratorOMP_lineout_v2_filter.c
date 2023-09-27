//
// Copyright (c) 2014, UChicago Argonne, LLC
// See LICENSE file.
//

// Integrator.c
//
// Hemant Sharma
// Dt: 2017/07/26
//
// TODO: Add option to give QbinSize instead of RbinSize

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


typedef double pixelvalue;

//~ #define PRINTOPT2
//~ #define PRINTOPT
#define SetBit(A,k)   (A[(k/32)] |=  (1 << (k%32)))
#define TestBit(A,k)  (A[(k/32)] &   (1 << (k%32)))
#define rad2deg 57.2957795130823
#define maxNFits 300
#define NrValsFitOutput 12

static inline double atand(double x){return rad2deg*(atan(x));}

static void
check (int test, const char * message, ...)
{
    if (test) {
        va_list args;
        va_start (args, message);
        vfprintf (stderr, message, args);
        va_end (args);
        fprintf (stderr, "\n");
        exit (EXIT_FAILURE);
    }
}

static inline
double**
allocMatrix(int nrows, int ncols)
{
    double** arr;
    int i;
    arr = malloc(nrows * sizeof(*arr));
    if (arr == NULL ) {
        return NULL;
    }
    for ( i = 0 ; i < nrows ; i++) {
        arr[i] = malloc(ncols * sizeof(*arr[i]));
        if (arr[i] == NULL ) {
            return NULL;
        }
    }
    return arr;
}

struct data {
    int y;
    int z;
    double frac;
};

struct data *pxList;
int *nPxList;

int ReadBins(){
    int fd;
    struct stat s;
    int status;
    size_t size;
    const char * file_name = "Map.bin";
    int rc;
    fd = open (file_name, O_RDONLY);
    check (fd < 0, "open %s failed: %s", file_name, strerror (errno));
    status = fstat (fd, & s);
    check (status < 0, "stat %s failed: %s", file_name, strerror (errno));
    size = s.st_size;
    int sizelen = 2*(int)sizeof(int) + (int)sizeof(double);
    //~ printf("Map size in bytes: %lld, each element size: %d, total elements: %lld. \n",(long long int)size,sizelen,(long long int)(size/sizelen));
    pxList = mmap (0, size, PROT_READ, MAP_SHARED, fd, 0);
    check (pxList == MAP_FAILED, "mmap %s failed: %s",file_name, strerror (errno));

    int fd2;
    struct stat s2;
    int status2;
    const char* file_name2 = "nMap.bin";
    fd2 = open (file_name2, O_RDONLY);
    check (fd2 < 0, "open %s failed: %s", file_name2, strerror (errno));
    status2 = fstat (fd2, & s2);
    check (status2 < 0, "stat %s failed: %s", file_name2, strerror (errno));
    size_t size2 = s2.st_size;
    nPxList = mmap (0, size2, PROT_READ, MAP_SHARED, fd2, 0);
    //printf("nMap size in bytes: %lld, each element size: %d, total elements: %lld. \n",(long long int)size2,2*(int)sizeof(int),2*(long long int)(size2/sizeof(int)));
    fflush(stdout);
    check (nPxList == MAP_FAILED, "mmap %s failed: %s",file_name, strerror (errno));
    return 1;
}

static inline
int StartsWith(const char *a, const char *b)
{
    if (strncmp(a,b,strlen(b)) == 0) return 1;
    return 0;
}

static inline void Transposer (double *x, int n1, int n2, double *y)
{
    int i,j;
    for (i=0;i<n1;i++){
        for (j=0;j<n2;j++){
            y[(i*n2)+j] = x[(j*n1)+i];
        }
    }
}

static inline
void
REtaMapper(
    double Rmin,
    double EtaMin,
    int nEtaBins,
    int nRBins,
    double EtaBinSize,
    double RBinSize,
    double *EtaBinsLow,
    double *EtaBinsHigh,
    double *RBinsLow,
    double *RBinsHigh)
{
    int i, j, k, l;
    for (i=0;i<nEtaBins;i++){
        EtaBinsLow[i] = EtaBinSize*i      + EtaMin;
        EtaBinsHigh[i] = EtaBinSize*(i+1) + EtaMin;
    }
    for (i=0;i<nRBins;i++){
        RBinsLow[i] =  RBinSize * i     + Rmin;
        RBinsHigh[i] = RBinSize * (i+1) + Rmin;
    }
}

static inline void DoImageTransformations (int NrTransOpt, int TransOpt[10], pixelvalue *ImageIn, pixelvalue *ImageOut, int NrPixelsY, int NrPixelsZ)
{
    int i,j,k,l,m;
    if (NrTransOpt == 0 || (NrTransOpt==1 && TransOpt[0]==0)){
        memcpy(ImageOut,ImageIn,NrPixelsY*NrPixelsZ*sizeof(*ImageIn)); // Nothing to do
        return;
    }
    for (i=0;i<NrTransOpt;i++){
        if (TransOpt[i] == 1){
            for (k=0;k<NrPixelsY;k++){
                for (l=0;l<NrPixelsZ;l++){
                    ImageOut[l*NrPixelsY+k] = ImageIn[l*NrPixelsY+(NrPixelsY-k-1)]; // Invert Y
                }
            }
        }else if (TransOpt[i] == 2){
            for (k=0;k<NrPixelsY;k++){
                for (l=0;l<NrPixelsZ;l++){
                    ImageOut[l*NrPixelsY+k] = ImageIn[(NrPixelsZ-l-1)*NrPixelsY+k]; // Invert Z
                }
            }
        }
    }
}

int fileReader (FILE *f,char fn[], int dType, int NrPixels, double *returnArr)
{
    int i;
    if (dType == 1){ // Binary with uint16
        uint16_t *readData;
        readData = calloc(NrPixels,sizeof(*readData));
        fread(readData,NrPixels*sizeof(*readData),1,f);
        for (i=0;i<NrPixels;i++){
            returnArr[i] = (double) readData[i];
        }
        free(readData);
        return 0;
    } else if (dType == 2){ // Binary with double
        double *readData;
        readData = calloc(NrPixels,sizeof(*readData));
        fread(readData,NrPixels*sizeof(*readData),1,f);
        for (i=0;i<NrPixels;i++){
            returnArr[i] = (double) readData[i];
        }
        free(readData);
        return 0;
    } else if (dType == 3){ // Binary with float
        float *readData;
        readData = calloc(NrPixels,sizeof(*readData));
        fread(readData,NrPixels*sizeof(*readData),1,f);
        for (i=0;i<NrPixels;i++){
            returnArr[i] = (double) readData[i];
        }
        free(readData);
        return 0;
    } else if (dType == 4){ // Binary with uint32
        uint32_t *readData;
        readData = calloc(NrPixels,sizeof(*readData));
        fread(readData,NrPixels*sizeof(*readData),1,f);
        for (i=0;i<NrPixels;i++){
            returnArr[i] = (double) readData[i];
        }
        free(readData);
        return 0;
    } else if (dType == 5){ // Binary with int32
        int32_t *readData;
        readData = calloc(NrPixels,sizeof(*readData));
        fread(readData,NrPixels*sizeof(*readData),1,f);
        for (i=0;i<NrPixels;i++){
            returnArr[i] = (double) readData[i];
        }
        free(readData);
        return 0;
    } else if (dType == 6){ // TIFF with uint32 format
        TIFFErrorHandler oldhandler;
        oldhandler = TIFFSetWarningHandler(NULL);
        printf("%s\n",fn);
        TIFF* tif = TIFFOpen(fn, "r");
        TIFFSetWarningHandler(oldhandler);
        if (tif){
            uint32 imagelength;
            tsize_t scanline;
            TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&imagelength);
            scanline = TIFFScanlineSize(tif);
            tdata_t buf;
            buf = _TIFFmalloc(scanline);
            uint32_t *datar;
            int rnr;
            for (rnr=0;rnr<imagelength;rnr++){
                TIFFReadScanline(tif,buf,rnr,1);
                datar = (uint32_t*)buf;
                for (i=0;i<scanline/sizeof(uint32_t);i++){
                    returnArr[rnr*(scanline/sizeof(uint32_t)) + i] = (double) datar[i];
                }
            }
            _TIFFfree(buf);
        }
        return 0;
    } else if (dType == 7){ // TIFF with uint8 format
        TIFFErrorHandler oldhandler;
        oldhandler = TIFFSetWarningHandler(NULL);
        printf("%s\n",fn);
        TIFF* tif = TIFFOpen(fn, "r");
        TIFFSetWarningHandler(oldhandler);
        if (tif){
            uint32 imagelength;
            tsize_t scanline;
            TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&imagelength);
            scanline = TIFFScanlineSize(tif);
            tdata_t buf;
            buf = _TIFFmalloc(scanline);
            uint8_t *datar;
            int rnr;
            for (rnr=0;rnr<imagelength;rnr++){
                TIFFReadScanline(tif,buf,rnr,1);
                datar = (uint8_t*)buf;
                for (i=0;i<scanline/sizeof(uint8_t);i++){
                    if (datar[i] == 1){
                        returnArr[rnr*(scanline/sizeof(uint8_t)) + i] = 1;
                    }
                }
            }
            _TIFFfree(buf);
        }
        return 0;
    } else {
        return 127;
    }
}


void BackgroundRemoval(int nRadFits,int nEtaFits,int nElsTot,double *peakIntensities,int width,int maxIteration,double *peakIntensity_without_bg)
{
    int inElsTot,inEtaFits,inRadFits,inElsPerRad,i,j;
    int nElsPerRad;
    nElsPerRad=(int)(nElsTot/nRadFits);
    //printf("%d\n",nElsPerRad);
    double peakIntensity_now[nElsPerRad], background_temp[nElsPerRad],background[nElsPerRad],peakIntensities_nElsTot[nElsTot];
    double peakIntensity_average,peakIntensity_min;
    
    for(inEtaFits=0;inEtaFits<nEtaFits;inEtaFits++){
        for(inElsTot=0;inElsTot<nElsTot;inElsTot++){peakIntensities_nElsTot[inElsTot]=peakIntensities[inEtaFits*nElsTot+inElsTot];}
        
        for(inRadFits=0;inRadFits<nRadFits;inRadFits++){
            for(inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){
                peakIntensity_now[inElsPerRad]=peakIntensities_nElsTot[inRadFits*nElsPerRad+inElsPerRad];
                background_temp[inElsPerRad]=0;
            }
            
            for (inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){
                if(peakIntensity_now[inElsPerRad]<peakIntensity_min){peakIntensity_min=peakIntensity_now[inElsPerRad];}
                peakIntensity_average+=peakIntensity_now[inElsPerRad];
            }
            peakIntensity_average=peakIntensity_average/nElsPerRad;
            
            for (inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){
                background[inElsPerRad]=peakIntensity_now[inElsPerRad];
                
                if (background[inElsPerRad]>peakIntensity_average+2*(peakIntensity_average-peakIntensity_min)){
                    background[inElsPerRad]=peakIntensity_average+2*(peakIntensity_average-peakIntensity_min);
                }
            }
            
            //for(inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){if(isnan(background_temp[inElsPerRad])==1){                    {printf("step1%d %lf %lf\n",inElsPerRad,background_temp[inElsPerRad],background[inElsPerRad]);}}}
            
            
            for(i=0;i<maxIteration;i++){
                for(inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){
                    if((inElsPerRad>=width)&&(inElsPerRad<=(nElsPerRad-width-1))){
                        for(j=0;j<width;j++){
                            
                               // if((i==0)){                    {printf("step1 %d %d %lf %lf\n",i,inElsPerRad,background_temp[inElsPerRad],background[inElsPerRad]);}}
                                
                            
                            
                            background_temp[inElsPerRad]=background_temp[inElsPerRad]+background[inElsPerRad-j]+background[inElsPerRad+j];
                            
                            {if((isnan(background_temp[inElsPerRad])==1)&&(i==0)){                    {printf("step2 %d %d %lf %lf %lf %lf\n",i,inElsPerRad,background_temp[inElsPerRad],background[inElsPerRad],background[inElsPerRad-j],background[inElsPerRad+j]);}}}
                        }
                        background_temp[inElsPerRad]=(1.0/(2*width))*background_temp[inElsPerRad];
                        //{if((isnan(background_temp[inElsPerRad])==1)&&(i==0)){                    {printf("step3 %d %d %lf %lf\n",i,inElsPerRad,background_temp[inElsPerRad],background[inElsPerRad]);}}}
                    }
                    if(inElsPerRad<width){
                        for(j=0;j<width;j++){background_temp[inElsPerRad]=background_temp[inElsPerRad]+background[inElsPerRad]+background[inElsPerRad+j];}
                        background_temp[inElsPerRad]=(1.0/(2*width))*background_temp[inElsPerRad];
                      
                    }
                    if(inElsPerRad>(nElsPerRad-width-1)){
                        for(j=0;j<width;j++){background_temp[inElsPerRad]=background_temp[inElsPerRad]+background[inElsPerRad-j]+background[inElsPerRad];}
                        background_temp[inElsPerRad]=(1.0/(2*width))*background_temp[inElsPerRad];
                    }
                 
                   
                    
                        
                    if(background_temp[inElsPerRad]>background[inElsPerRad]){background_temp[inElsPerRad]=background[inElsPerRad];}
                    
                    if(background_temp[inElsPerRad]<0.1){background_temp[inElsPerRad]=0.1;}
                   // printf("%d %lf %lf %lf\n",inElsPerRad,peakIntensity_now[inElsPerRad],background[inElsPerRad],background_temp[inElsPerRad]);
                    
                }
                //for(inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){{printf("step2 %d %lf %lf %lf\n",inElsPerRad,peakIntensity_now[inElsPerRad],background[inElsPerRad],background_temp[inElsPerRad]);}}
                
                for(inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){background[inElsPerRad]=background_temp[inElsPerRad];background_temp[inElsPerRad]=0;}
                
                //for(inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){{printf("step2 %d %lf %lf %lf\n",inElsPerRad,peakIntensity_now[inElsPerRad],background[inElsPerRad],background_temp[inElsPerRad]);}}
                
                
            }
            for(inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){peakIntensity_without_bg[inEtaFits*nElsTot+(inRadFits*nElsPerRad+inElsPerRad)]=peakIntensity_now[inElsPerRad]-background[inElsPerRad];

            }
            //for(inElsPerRad=0;inElsPerRad<nElsPerRad;inElsPerRad++){if(isnan(peakIntensity_without_bg[inElsPerRad])==1){printf("step3%d %lf %lf\n",inElsPerRad,peakIntensity_now[inElsPerRad],background[inElsPerRad]);}}
        }
    }
}


int integrator(char *ParamFN, char *darkFN, char *imageFN, int nFrames)
{
    
    clock_t start, end, start0, end0;
    start0 = clock();
    double diftotal;
    double RMax, RMin, RBinSize, EtaMax, EtaMin, EtaBinSize, Lsd, px;
    int NrPixelsY = 2048, NrPixelsZ = 2048, Normalize = 1;
    int nEtaBins, nRBins;
    FILE *paramFile;
    char aline[4096], dummy[4096], *str;
    paramFile = fopen(ParamFN,"r");
    int HeadSize = 8192;
    int NrTransOpt=0;
    long long int GapIntensity=0, BadPxIntensity=0;
    int TransOpt[10];
    int makeMap = 0;
    size_t mapMaskSize = 0;
    int *mapMask;
    int dType = 1;
    int filter=0;
    char GapFN[4096], BadPxFN[4096], outputFolder[4096], timestr[4096];
    int sumImages=0, separateFolder=0,newOutput=0, binOutput = 0;
    double radiiToFit[maxNFits][6], etasToFit[maxNFits][4];
    size_t nRadFits = 0, nEtaFits = 0;
    int multipeak,dointegration,dolineout;
    while (fgets(aline,4096,paramFile) != NULL){
        str = "GapFile ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, GapFN);
            makeMap = 2;
        }
        str = "OutFolder_Integration ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, outputFolder);
            separateFolder = 1;
        }
        str = "timestr ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, timestr);
        }
        str = "BadPxFile ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, BadPxFN);
            makeMap = 2;
        }
        str = "RadiusToFit ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf %lf", dummy, &radiiToFit[nRadFits][0], &radiiToFit[nRadFits][1]);
            radiiToFit[nRadFits][2] = radiiToFit[nRadFits][0] - radiiToFit[nRadFits][1];
            radiiToFit[nRadFits][3] = radiiToFit[nRadFits][0] + radiiToFit[nRadFits][1];
            nRadFits++;
        }
        str = "EtaToFit ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf %lf", dummy, &etasToFit[nEtaFits][0], &etasToFit[nEtaFits][1]);
            etasToFit[nEtaFits][2] = etasToFit[nEtaFits][0] - etasToFit[nEtaFits][1];
            etasToFit[nEtaFits][3] = etasToFit[nEtaFits][0] + etasToFit[nEtaFits][1];
            nEtaFits++;
        }
        str = "EtaBinSize ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &EtaBinSize);
        }
        str = "RBinSize ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &RBinSize);
        }
        str = "DataType ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &dType);
        }
        str = "HeadSize ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &HeadSize);
        }
        str = "RMax ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &RMax);
        }
        str = "RMin ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &RMin);
        }
        str = "EtaMax ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &EtaMax);
        }
        str = "EtaMin ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &EtaMin);
        }
        str = "Lsd ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &Lsd);
        }
        str = "px ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lf", dummy, &px);
        }
        str = "NrPixelsY ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &NrPixelsY);
        }
        str = "NrPixelsZ ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &NrPixelsZ);
        }
        str = "Normalize ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &Normalize);
        }
        str = "NrPixels ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &NrPixelsY);
            sscanf(aline,"%s %d", dummy, &NrPixelsZ);
        }
        str = "GapIntensity ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lld", dummy, &GapIntensity);
            makeMap = 1;
            continue;
        }
        str = "BadPxIntensity ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %lld", dummy, &BadPxIntensity);
            makeMap = 1;
            continue;
        }
        str = "ImTransOpt ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &TransOpt[NrTransOpt]);
            NrTransOpt++;
            continue;
        }
        str = "SumImages ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &sumImages);
            continue;
        }
        str = "multipeak ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &multipeak);
        }
        str = "doIntegration ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &dointegration);
        }
        str = "doLineout ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &dolineout);
        }
        str = "Filter ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &filter);
        }
    }
   

    if (separateFolder == 0){
        sprintf(outputFolder,".");
        separateFolder = 1;
    }

    double *AreaMapPixels;
    size_t nFits = nEtaFits * nRadFits;
    size_t iRadFit, iEtaFit, nEls, nElsTot=0;
    for (iRadFit=0;iRadFit<nRadFits;iRadFit++){
        nEls = (int)(ceil(radiiToFit[iRadFit][1]*2 / RBinSize)) + 2;
        radiiToFit[iRadFit][4] = nElsTot;
        radiiToFit[iRadFit][5] = nEls;
        nElsTot += nEls;
    }
    
    AreaMapPixels = calloc(NrPixelsY*NrPixelsZ,sizeof(*AreaMapPixels));
   
    nRBins = (int) ceil((RMax-RMin)/RBinSize);
    nEtaBins = (int)ceil((EtaMax - EtaMin)/EtaBinSize);
    double *EtaBinsLow, *EtaBinsHigh;
    double *RBinsLow, *RBinsHigh;
    EtaBinsLow = malloc(nEtaBins*sizeof(*EtaBinsLow));
    EtaBinsHigh = malloc(nEtaBins*sizeof(*EtaBinsHigh));
    RBinsLow = malloc(nRBins*sizeof(*RBinsLow));
    RBinsHigh = malloc(nRBins*sizeof(*RBinsHigh));
    REtaMapper(RMin, EtaMin, nEtaBins, nRBins, EtaBinSize, RBinSize, EtaBinsLow, EtaBinsHigh, RBinsLow, RBinsHigh);

    int i,j,k,l;
    //~ printf("NrTransOpt: %d\n",NrTransOpt);
    for (i=0;i<NrTransOpt;i++){
        if (TransOpt[i] < 0 || TransOpt[i] > 2){printf("TransformationOptions can only be 0, 1, 2.\nExiting.\n");return 0;}
        //~ printf("TransformationOptions: %d ",TransOpt[i]);
        //~ if (TransOpt[i] == 0) printf("No change.\n");
        //~ else if (TransOpt[i] == 1) printf("Flip Left Right.\n");
        //~ else if (TransOpt[i] == 2) printf("Flip Top Bottom.\n");
    }
    double *Image;
    pixelvalue *ImageIn;
    pixelvalue *DarkIn;
    pixelvalue *ImageInT;
    pixelvalue *DarkInT;
    double *AverageDark;
    DarkIn = malloc(NrPixelsY*NrPixelsZ*sizeof(*DarkIn));
    DarkInT = malloc(NrPixelsY*NrPixelsZ*sizeof(*DarkInT));
    AverageDark = calloc(NrPixelsY*NrPixelsZ,sizeof(*AverageDark));
    ImageIn = malloc(NrPixelsY*NrPixelsZ*sizeof(*ImageIn));
    ImageInT = malloc(NrPixelsY*NrPixelsZ*sizeof(*ImageInT));
    Image = malloc(NrPixelsY*NrPixelsZ*sizeof(*Image));
    size_t pxSize;
    if (dType == 1){ // Uint16
        pxSize = sizeof(uint16_t);
    } else if (dType == 2){ // Double
        pxSize = sizeof(double);
    } else if (dType == 3){ // Float
        pxSize = sizeof(float);
    } else if (dType == 4){ // Uint32
        pxSize = sizeof(uint32_t);
    } else if (dType == 5){ // Int32
        pxSize = sizeof(int32_t);
    } else if (dType == 6){ // Tiff Uint32
        pxSize = sizeof(uint32_t);
        HeadSize = 0;
    } else if (dType == 7){ // Tiff Uint8
        pxSize = sizeof(uint8_t);
        HeadSize = 0;
    }
    size_t SizeFile = pxSize * NrPixelsY * NrPixelsZ;
    //~ int nFrames;
    size_t sz;
    int Skip = HeadSize;
    FILE *fp, *fd;
    int nrdone = 0;
    fd = fopen(darkFN,"rb");
    fseek(fd,0L,SEEK_END);
    sz = ftell(fd);
    rewind(fd);
    int nFramesDark = sz / (SizeFile);
    printf("Reading dark file:      %s, nFrames: %d, skipping first %d bytes.\n",darkFN,nFrames,Skip);
    fseek(fd,Skip,SEEK_SET);
    for (i=0;i<nFramesDark;i++){
        int retCode = fileReader(fd,darkFN,dType,NrPixelsY*NrPixelsZ,DarkInT);
        DoImageTransformations(NrTransOpt,TransOpt,DarkInT,DarkIn,NrPixelsY,NrPixelsZ);
        if (makeMap == 1){
            mapMaskSize = NrPixelsY;
            mapMaskSize *= NrPixelsZ;
            mapMaskSize /= 32;
            mapMaskSize ++;
            mapMask = calloc(mapMaskSize,sizeof(*mapMask));
            for (j=0;j<NrPixelsY*NrPixelsZ;j++){
                if (DarkIn[j] == (pixelvalue) GapIntensity || DarkIn[j] == (pixelvalue) BadPxIntensity){
                    SetBit(mapMask,j);
                    nrdone++;
                }
            }
            //~ printf("Nr mask pixels: %d\n",nrdone);
            makeMap = 0;
        }
        for(j=0;j<NrPixelsY*NrPixelsZ;j++) AverageDark[j] += (double)DarkIn[j]/nFramesDark;
    }
    printf("Dark file read\n");
    if (makeMap == 2){
        mapMaskSize = NrPixelsY;
        mapMaskSize *= NrPixelsZ;
        mapMaskSize /= 32;
        mapMaskSize ++;
        mapMask = calloc(mapMaskSize,sizeof(*mapMask));
        double *mapper;
        mapper = calloc(NrPixelsY*NrPixelsZ,sizeof(*mapper));
        double *mapperOut;
        mapperOut = calloc(NrPixelsY*NrPixelsZ,sizeof(*mapperOut));
        fileReader(fd,GapFN,7,NrPixelsY*NrPixelsZ,mapper);
        DoImageTransformations(NrTransOpt,TransOpt,mapper,mapperOut,NrPixelsY,NrPixelsZ);
        for (i=0;i<NrPixelsY*NrPixelsZ;i++){
            if (mapperOut[i] != 0){
                SetBit(mapMask,i);
                mapperOut[i] = 0;
                nrdone++;
            }
        }
        fileReader(fd,BadPxFN,7,NrPixelsY*NrPixelsZ,mapper);
        DoImageTransformations(NrTransOpt,TransOpt,mapper,mapperOut,NrPixelsY,NrPixelsZ);
        for (i=0;i<NrPixelsY*NrPixelsZ;i++){
            if (mapperOut[i] != 0){
                SetBit(mapMask,i);
                mapperOut[i] = 0;
                nrdone++;
            }
        }
        free(mapper);
        free(mapperOut);
        printf("Nr mask pixels: %d\n",nrdone);
    }
    fp = fopen(imageFN,"rb");
    fseek(fp,0L,SEEK_END);
    sz = ftell(fp);
    rewind(fp);
    fseek(fp,Skip,SEEK_SET);
    //~ nFrames = (sz-Skip) / SizeFile;
    //~ printf("Number of eta bins: %d, number of R bins: %d. Number of frames in the file: %d\n",nEtaBins,nRBins,(int)nFrames);
    long long int Pos;
    int nPixels, dataPos;
    struct data ThisVal;
    char outfn[4096];
    char outfn2[4096];
    FILE *out,*out2;
    char outFN1d[4096];
    char dmyt[10000];
    FILE *out1d;
    double Intensity, totArea, ThisInt;
    size_t testPos;
    double RMean, EtaMean;
    double Int1d;
    int n1ds;
    double *sumMatrix;
   
    double *outArr, *outThisArr, *out1dArr;
    char *outext;
    outext = ".csv";
    size_t bigArrSize = nEtaBins*nRBins;
    double *IntArrPerFrame;
    IntArrPerFrame = calloc(bigArrSize,sizeof(*IntArrPerFrame));
    FILE *out3, *outPeakFit,*out4;
    int etaFitNr, radFitNr, radBinFitNr,nEtaF;
    size_t PeakIntPos, peakPos, nElsThis, peakValArrPos;
    double *RPosArr, *RFitThis, *IntensityThis, *ResultArr;
    RPosArr = calloc(nElsTot,sizeof(*RPosArr));
    ResultArr = malloc(NrValsFitOutput*sizeof(*ResultArr));
    FILE *outPeak,*outPeakwithoutbg,*outPeakmedian;
    FILE *outPeak2theta,*outPeakwithoutbg2theta,*outPeakmedian2theta;
    char outfnAll[4096];
    char outfnAll_1d[4096];
    char outfnFit[4096];
    char fn3[4096];
    sprintf(fn3,"%s",imageFN);
    char *bname3;
    bname3 = basename(fn3);
    
    if (dointegration == 1){
        sprintf(outfnAll,"%s/Integrated/%s_integrated_%s.bin",outputFolder,bname3,timestr);
        out3 = fopen(outfnAll,"wb");
        sprintf(outfnAll_1d,"%s/Integrated_1d/%s_integrated_1d_%s.bin",outputFolder,bname3,timestr);
        out4 = fopen(outfnAll_1d,"wb");}
    
    if(dolineout==1){
        char outpeakfn[4096];
        sprintf(outpeakfn,"%s/Lineout_mean/%s_LineOuts_mean_%s.bin",outputFolder,bname3,timestr);
        outPeak = fopen(outpeakfn,"wb");
        char outpeakwithoutbgfn[4096];
        sprintf(outpeakwithoutbgfn,"%s/Lineout_mean_withoutbg/%s_LineOuts_mean_withoutbg_%s.bin",outputFolder,bname3,timestr);
        outPeakwithoutbg = fopen(outpeakwithoutbgfn,"wb");
        char outpeakfn_median[4096];
        sprintf(outpeakfn_median,"%s/Lineout_median/%s.LineOuts_median_%s.bin",outputFolder,bname3,timestr);
        outPeakmedian = fopen(outpeakfn_median,"wb");
        
        //char outpeakfn2theta[4096];
        //sprintf(outpeakfn2theta,"%s/Lineout_mean_2theta/%s_LineOuts_mean_2theta_%s.txt",outputFolder,bname3,timestr);
        //outPeak2theta = fopen(outpeakfn2theta,"wb");
        //char outpeakwithoutbgfn2theta[4096];
        //sprintf(outpeakwithoutbgfn2theta,"%s/Lineout_mean_withoutbg_2theta/%s_LineOuts_mean_2theta_withoutbg_%s.txt",outputFolder,bname3,timestr);
        //outPeakwithoutbg2theta = fopen(outpeakwithoutbgfn2theta,"wb");
        //char outpeakfn_median2theta[4096];
       // sprintf(outpeakfn_median2theta,"%s/Lineout_median_2theta/%s.LineOuts_median_2theta_%s.txt",outputFolder,bname3,timestr);
        //outPeakmedian2theta = fopen(outpeakfn_median2theta,"wb");
        
        
        
    }
        
    printf("Processing file %s\n",imageFN);
    for (i=0;i<nFrames;i++){
        int rcode = fileReader(fp,imageFN,dType,NrPixelsY*NrPixelsZ,ImageInT);
        DoImageTransformations(NrTransOpt,TransOpt,ImageInT,ImageIn,NrPixelsY,NrPixelsZ);
        for (j=0;j<NrPixelsY*NrPixelsZ;j++){
            Image[j] = (double)ImageIn[j] - AverageDark[j];
            if(Image[j]<0)
            {Image[j]=0;}
        }
        
        if (i==0 && dointegration==1){
            if (separateFolder==0) sprintf(outfn2,"%s.REtaAreaMap_%s.csv",imageFN,timestr);
            else{
                char fn2[4096];
                sprintf(fn2,"%s",imageFN);
                char *bnname;
                bnname = basename(fn2);
                sprintf(outfn2,"%s/REtaMap/%s.REtaAreaMap_%s.csv",outputFolder,bnname,timestr);
            }
            out2 = fopen(outfn2,"w");
            fprintf(out2,"%%nEtaBins:\t%d\tnRBins:\t%d\n%%Radius(px)\t2Theta(degrees)\tEta(degrees)\tBinArea\n",nEtaBins,nRBins);
        }
        
        double *peakIntensities,*peakIntensity_without_bg;
        
        double **peakIntensitiesarray;
        
        int *NrofIntg;
        int CurrentNr=0;
        int counthere;
        int countherej,countherek;
        double *peakIntensities_median;
        int defaultcount=2000;
        peakIntensitiesarray= (double**)malloc(sizeof(double*)*nElsTot*nEtaFits);
        for (counthere=0;counthere<nElsTot*nEtaFits;counthere++)
        {
            peakIntensitiesarray[counthere]=(double*)malloc(sizeof(double)*defaultcount);
           
        }
        
        
        peakIntensities = calloc(nElsTot*nEtaFits,sizeof(*peakIntensities));
        peakIntensity_without_bg=calloc(nElsTot*nEtaFits,sizeof(*peakIntensity_without_bg));
        peakIntensities_median= calloc(nElsTot*nEtaFits,sizeof(*peakIntensities_median));
        NrofIntg=calloc(nElsTot*nEtaFits,sizeof(*NrofIntg));
        
        
        
        memset(IntArrPerFrame,0,bigArrSize*sizeof(double));
        memset(peakIntensities,0,nElsTot*nEtaFits*sizeof(double));
        memset(peakIntensity_without_bg,0,nElsTot*nEtaFits*sizeof(double));
        memset(NrofIntg,0,nElsTot*nEtaFits*sizeof(int));
        memset(peakIntensities_median,0,nElsTot*nEtaFits*sizeof(double));
        for(counthere=0;counthere<nElsTot*nEtaFits;counthere++)
        {
            for(countherej=0;countherej<(defaultcount);countherej++)
            {
                peakIntensitiesarray[counthere][countherej]=0;
            }
        }
        
        
        for (j=0;j<nRBins;j++){
            RMean = (RBinsLow[j]+RBinsHigh[j])/2;
            Int1d = 0;
            n1ds = 0;
            radFitNr = -1;
            for (k=0;k<nRadFits;k++){
                if (RBinsHigh[j] >= radiiToFit[k][2] && RBinsLow[j] <= radiiToFit[k][3]){
                    radFitNr = k;
                    radBinFitNr = (int) ((RBinsHigh[j]-radiiToFit[k][2])/RBinSize);
                    RPosArr[(int)radiiToFit[k][4]+ radBinFitNr] = RMean;
                    if (radBinFitNr >= (int)radiiToFit[k][5]){
                        printf("Something went wrong in fitting calculation, exiting %d %lf.\n",radBinFitNr,radiiToFit[k][5]);
                        return 1;
                    }
                }
            }
            for (k=0;k<nEtaBins;k++){
                etaFitNr = -1;
                if (radFitNr > -1){
                    for (nEtaF=0;nEtaF<nEtaFits;nEtaF++){
                        if (EtaBinsHigh[k] >= etasToFit[nEtaF][2] && EtaBinsLow[k] <= etasToFit[nEtaF][3]){
                            etaFitNr = nEtaF;
                        }
                    }
                }
                Pos = j*nEtaBins + k;
                nPixels = nPxList[2*Pos + 0];
                dataPos = nPxList[2*Pos + 1];
                Intensity = 0;
                totArea = 0;
                for (l=0;l<nPixels;l++){
                    ThisVal = pxList[dataPos + l];
                    testPos = ThisVal.z;
                    testPos *= NrPixelsY;
                    testPos += ThisVal.y;
                    if (mapMaskSize!=0){
                        if (TestBit(mapMask,testPos)){
                            continue;
                        }
                    }
                    if (i==0){
                        AreaMapPixels[testPos] += ThisVal.frac;
                    }
                    ThisInt = Image[testPos]; // The data is arranged as y(fast) and then z(slow)
                    Intensity += ThisInt*ThisVal.frac;
                    totArea += ThisVal.frac;
                }
                if (Intensity != 0){
                    if (Normalize == 1){
                        Intensity /= totArea;
                    }
                }
                if (etaFitNr > -1){
                    PeakIntPos = nElsTot;
                    PeakIntPos *= etaFitNr;
                    PeakIntPos += (int)radiiToFit[radFitNr][4];
                    PeakIntPos += radBinFitNr;
                    peakIntensities[PeakIntPos] += Intensity;
                    NrofIntg[PeakIntPos]+=1;
                    CurrentNr=NrofIntg[PeakIntPos];
                    peakIntensitiesarray[PeakIntPos][CurrentNr]=Intensity;
                    
                }
                EtaMean = (EtaBinsLow[k]+EtaBinsHigh[k])/2;
                Int1d += Intensity;
                n1ds ++;
              if(dointegration==1){
                    if (i==0){
                        fprintf(out2,"%lf\t%lf\t%lf\t%lf\n",RMean,atand(RMean*px/Lsd),EtaMean,totArea);
                    }
                    IntArrPerFrame[j*nEtaBins+k] = Intensity;
                }
              
            }
            Int1d /= n1ds;
            if (dointegration == 1) fprintf(out4,"%lf %lf %lf\n",RMean,atand(RMean*px/Lsd),Int1d);
        }
        if(filter==1){
            
            int m;
            int middlenumber;
            int realcount;
            realcount=CurrentNr;
            for(counthere=0;counthere<nElsTot*nEtaFits;counthere++)
            {
                for(countherek=0;countherek<(realcount-1);countherek++)
                {
                    for(countherej=0;countherej<(realcount-1-countherek);countherej++)
                    {
                        
                        if(peakIntensitiesarray[counthere][countherej]>peakIntensitiesarray[counthere][countherej+1])
                        {
                            m=peakIntensitiesarray[counthere][countherej];
                            peakIntensitiesarray[counthere][countherej]=peakIntensitiesarray[counthere][countherej+1];
                            peakIntensitiesarray[counthere][countherej+1]=m;
                        }
                        
                    }
                }
            }
            
            
            
            if(realcount%2==0)
            {
                middlenumber=realcount/2;
                for(counthere=0;counthere<nElsTot*nEtaFits;counthere++)
                {
                    peakIntensities_median[counthere]=(peakIntensitiesarray[counthere][middlenumber-1]+peakIntensitiesarray[counthere][middlenumber])/2.0;
                }
            }
            if(realcount%2==1)
            {
                middlenumber=(defaultcount-1)/2;
                for(counthere=0;counthere<nElsTot*nEtaFits;counthere++)
                {
                    peakIntensities_median[counthere]=peakIntensitiesarray[counthere][middlenumber];
                }
            }
        }
        
        if (dointegration == 1){
            fwrite(IntArrPerFrame,bigArrSize*sizeof(*IntArrPerFrame),1,out3);}
        
        if (i==0){
            if(dointegration==1)
            {
                FILE *areamap;
                char areamapfn[4096];
                sprintf(areamapfn,"%s/AreaFractionPixels_%s.bin",outputFolder,timestr);
                areamap = fopen(areamapfn,"wb");
                fwrite(AreaMapPixels,NrPixelsY*NrPixelsZ*sizeof(double),1,areamap);
                fclose(areamap);
                free(AreaMapPixels);
                fclose(out2);
            }
            if(dolineout==1)
            {
                //~ printf("%zu\n",sizeof(size_t));
                //fwrite(&nElsTot,sizeof(size_t),1,outPeakwithoutbg);
               // fwrite(&nEtaFits,sizeof(size_t),1,outPeakwithoutbg);
               // fwrite(&nRadFits,sizeof(size_t),1,outPeakwithoutbg);
            }
        }
        if(dolineout==1){
            //printf("%d,%d,%d\n",i,nElsTot,nEtaFits);
            fwrite(peakIntensities,nElsTot*nEtaFits*sizeof(double),1,outPeak);
            if(filter==1){
                fwrite(peakIntensities_median,nElsTot*nEtaFits*sizeof(double),1,outPeakmedian);
            }
            BackgroundRemoval(nRadFits,nEtaFits,nElsTot,peakIntensities,20,100,peakIntensity_without_bg);
            fwrite(peakIntensity_without_bg,nElsTot*nEtaFits*sizeof(double),1,outPeakwithoutbg);
                        
        }
        
        int countF=0;
      //  if(dolineout==1){
            //printf("%d,%d,%d\n",i,nElsTot,nEtaFits);
            
           // for(countF=0;countF<nElsTot*nEtaFits;countF++)
           // {
              //  fprintf(outPeak2theta,"%lf %lf\n",atand((RMin+(countF-1)*RBinSize)*px/Lsd),peakIntensities[countF]);
              //  fprintf(outPeakmedian2theta,"%lf %lf\n",atand((RMin+(countF-1)*RBinSize)*px/Lsd),peakIntensities_median[countF]);
               // fprintf(outPeakwithoutbg2theta,"%lf %lf\n",atand((RMin+(countF-1)*RBinSize)*px/Lsd),peakIntensity_without_bg[countF]);
    
           // }
       // }
        
       
        
        
        free(peakIntensities);
        free(peakIntensity_without_bg);
        free(peakIntensities_median);
        for (counthere=0;counthere<nElsTot*nEtaFits;counthere++)
        free(peakIntensitiesarray[counthere]);
        free(peakIntensitiesarray);
        
    }
    if (dointegration == 1){
        
        fclose(out3);
        fclose(out4);
    }
    if(dolineout==1){
        fclose(outPeak);
        fclose(outPeakwithoutbg);
        fclose(outPeakmedian);
      //  fclose(outPeak2theta);
       // fclose(outPeakmedian2theta);
      //  fclose(outPeakwithoutbg2theta);
        
    }
   
    end0 = clock();
    diftotal = ((double)(end0-start0))/CLOCKS_PER_SEC;
    //~ printf("Total time elapsed:\t%f s.\n",diftotal);
    
    free(EtaBinsLow);
    free(EtaBinsHigh);
    free(RBinsLow);
    free(RBinsHigh);
    free(DarkIn);
    free(DarkInT);
    free(AverageDark);
    free(ImageIn);
    free(ImageInT);
    free(Image);
    free(IntArrPerFrame);
    free(RPosArr);
    free(ResultArr);
   
    if (mapMaskSize != 0) free(mapMask);
   
    return 0;
}


int main(int argc, char **argv)
{
    
    if(argc!=2)
    {
        printf("Usage: ./IntegratorOMP_lineout_v2 ParamFN\n");
        return 1;
    }
    
    system("cp Map.bin nMap.bin /dev/shm");
    char *pfn = argv[1];
    char FileStem[4096],darkFN[4096],ext[4096],RawFolder[4096];
    int startNr,endNr;
    int pad;
    int nFrames;
    int numProcs;
    char aline[4096], dummy[4096], *str, outputFolder[4096];
    FILE *paramFile;
    paramFile = fopen(pfn,"r");
  
    
    while (fgets(aline,4096,paramFile) != NULL){
        
        str = "RawFolder ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, FileStem);
                    }
        str = "FileStem ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, RawFolder);
                    }
       
        str = "Dark ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %s", dummy, darkFN);
                   }
        str = "Ext ";
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
        str = "Padding ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &pad);
                    }
        str = "nFrames ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &nFrames);
        }
        str = "numProcs ";
        if (StartsWith(aline,str) == 1){
            sscanf(aline,"%s %d", dummy, &numProcs);
        }
    }
 
    int nFiles = (endNr - startNr + 1);
    strcat(FileStem,RawFolder);
    int rc = ReadBins();
    int fileNr;
    
#pragma omp parallel for num_threads(numProcs) private(fileNr) schedule(dynamic)
    for (fileNr=startNr;fileNr<=endNr;fileNr++)
    {
       
        char FN[4096];
        sprintf(FN,"%s_%0*d%s",FileStem,pad,fileNr,ext);
       
        int rt = integrator(pfn,darkFN,FN,nFrames);
    }
    return 0;
}

