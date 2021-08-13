/*-----------------Fast double-difference cross-correlation (FDTCC)-----------
                           Min Liu & Miao Zhang
                       m.liu@dal.ca & miao.zhang@dal.ca
                             Dalhouise University         
                                 Nov. 15 2019		

# create the dt.cc file from raw continuous SAC file or cut SAC file
#
#  -------|-----------|------------------|------
#             (wb)   pick      (wa)
#
# For P phase CC, if your "wa" is larger than 0.9*(ts-tp), it will be replaced 
# by 0.9*(ts-tp) to make sure you don't include S phase.
# For S phase CC, if your "wb" is larger than 0.5*(ts-tp), it will be replaced
# by 0.5*(ts-tp) to make sure you don't include P phase.
---------------------------------------------------------------------------*/
#include "sac.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// define name of ouput files
#define INPUT1 "Input.p"
#define INPUT2 "Input.s1"
#define INPUT3 "Input.s2"
#define OUTPUT "dt.cc"
#define SAC_BANDPASS "BP"
#define SAC_BUTTERWORTH "BU"

typedef struct picks {
    char sta[8];
    char phase[5];
    float arr1;
    float arr2;
    double ccv;
    float shift;
    int quality;
    float diff;
} PICK;
typedef struct event_pair {
    int event1;
    int event2;
    PICK* pk;

} PAIR;

typedef struct stationinfo {
    float stlo;
    float stla;
    char net[5];
    char sta[8];
    char comp[4];
    float elev;
} STATION;

typedef struct eventinfo {
    char date[10];
    char time[10];
    double sec;
    float evlo;
    float evla;
    float evdp;
    int event;
    float* pSNR;
    float* sSNR;
} EVENT;

typedef struct phases {
    char sta[10];
    char phase[5];
} PHASE;

typedef struct phaseinfo {
    float* time;
    PHASE* pa;
} PHASEINF;

typedef struct ttable {
    double gdist;
    double dep;
    double ptime;
    double stime;
    double prayp;
    double srayp;
    double phslow;
    double shslow;
    char pphase[10];
    char sphase[10];
} TTT;

// declaration of function
void SubccS(PAIR*, float**, float**, int*, int, int);
void SubccP(PAIR*, float**, int*, int, int);
void Correct_Sshift(PAIR*, float*, char**, int*);
void Correct_Pshift(PAIR*, float*, char**, int*);
void Transfer_sec(EVENT*, int);
void Cal_tt(PAIR*, PAIR*, EVENT*, STATION*);
void Search_event(PAIR*, EVENT*, int*, int);
void ddistaz(double, double, double, double, double*);
void Cal_sSNR(float**, float**, int*);
void Cal_pSNR(float**, int*);
void Replace(PAIR*, PHASEINF*, int, int);
void taper(float*, int, float, float);
void rtrend(float*, int);
void bpcc(float*, SACHEAD, float, float);
void xapiir(float*, int, char*, double, double, int, char*, double, double, double, int);

TTT* TB;
STATION* ST;
PAIR* PO;
PAIR* PT;
EVENT* EVE;
PHASEINF* PIN;

// globe parameters
int NP = 1700000;
int NS = 120;
int NE = 30000;
int np; // number of event pairs
int ns; // number of stations
int ne; // number of events
int ntb = 100000; // number of lines in ttt
float wb;
float wa;
float wf;
float wbs;
float was;
float wfs;
float delta = 0.01;
float threshold = 0.5;
float trx;
float tdx;
float tdh;
float trh;
float thre_SNR = 1.0;
float thre_shift = 1.5;
float timezone = 0 * 3600; //not tested

int main(int argc, char** argv)
{
    FILE *fp1, *fp2, *fp3;
    int i, j, error = 0;
    float low, high;
    char lines[500], flag[10];
    char staDir[100], tttDir[100], wavDir[100], eveDir[100], dctDir[100],
        paDir[100];
    float jk, jk1, jk2, jk3, jk4, jk5, jk6, jk7, jk8, jk9, jk10, jk11, jk12, jk13;
    int ID_event;
    int k = 0, kk = 0;
    int f = -3;
    char **staP, **staS1, **staS2;
    float *ptriger, *s1triger, *s2triger;
    float **waveP, **waveS1, **waveS2;
    int *labelP, *labelS1, *labelS2, *markP, *markS1, *markS2;
    double memory_require;
    char **la_staP, **la_staS1, **la_staS2;
    SACHEAD hd1, hd2, hd3;
    int size, npp, nss;
    int ife, ifd, ifp;
    extern int np, ne, ns, ntb;
    extern float threshold, delta, trx, tdx, trh, tdh, wa, wb, wf, was, wbs, wfs;

    // read parametets
    for (i = 1; !error && i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
            case 'C':
                sscanf(&argv[i][2], "%d/%d/%d", &ife, &ifd, &ifp);
            case 'W':
                sscanf(&argv[i][2], "%f/%f/%f/%f/%f/%f", &wb, &wa, &wf, &wbs, &was,
                    &wfs);
                break;
            case 'D':
                sscanf(&argv[i][2], "%f/%f/%f/%f", &delta, &threshold, &thre_SNR,
                    &thre_shift);
                break;
            case 'G':
                sscanf(&argv[i][2], "%f/%f/%f/%f", &trx, &trh, &tdx, &tdh);
                break;
            case 'B':
                sscanf(&argv[i][2], "%f/%f", &low, &high);
                break;
            case 'F':
                sscanf(&argv[i][2], "%d", &f);
                break;
            default:
                error = 1;
                break;
            }
        }
    }

    if (argc < 9 || error == 1) {
        fprintf(stderr, "Usage:  Fast double-difference cross-correlation (FDTCC) \n");
        fprintf(stderr, "        Create the dt.cc file for hypoDD, tomoDD, Growclust, etc.");
        fprintf(stderr, "(Authors: Min Liu & Miao Zhang)\n");
        fprintf(stderr,
            "	FDTCC -C(ife/ifd/ifp) -W(wb/wa/wf/wbs/was/wfs) "
            "-D(delta/threshold/thre_SNR/thre_shift)\n"
            "        -G(trx/trh/tdx/tdh) -B(low/high) -F(f)\n");
        fprintf(stderr, "                   "
                        "-------------------explanation------------------\n");
        fprintf(stderr, "	-C: specify the path of event.sel, dt.ct and "
                        "phase.dat (1: yes, 0: default names)\n");
        fprintf(stderr, "	-W: waveform window length before and after "
                        "picks and their maximum shift length\n");
        fprintf(stderr, "	-D: sampling interval, CC threshold, SNR threshold, max arrival time diff of the pick pair\n");
        fprintf(stderr, "	-G: ranges and grids in horizontal direction and depth (in traveltime table)\n");
        fprintf(stderr, "	-F: input data format (0: continuous data; 1: event segments)\n");
        fprintf(stderr, "	-B: waveform bandpass filtering (e.g., 2/8; "
                        "-1/-1: no filter applied).\n");
        fprintf(stderr, "         SAC name format: date/net.sta.comp, e.g., 20210101/AA.BBBB.HHZ\n"
                        "                           or eventID/net.sta.comop, e.g., 8/AA.BBBB.HHZ).\n");
        exit(-1);
    }
    strcpy(staDir, argv[7]);
    strcpy(tttDir, argv[8]);
    strcpy(wavDir, argv[9]);
    if (ife == 1)
        strcpy(eveDir, argv[10]);
    else
        strcpy(eveDir, "./event.sel");
    if (ifd == 1)
        strcpy(dctDir, argv[11]);
    else
        strcpy(dctDir, "./dt.ct");
    if (ifp == 1)
        strcpy(paDir, argv[12]);
    else
        strcpy(paDir, "./phase.dat");
    // read phase inforamtion phase.dat
    PIN = (PHASEINF*)malloc(sizeof(PHASEINF) * (NE));
    for (i = 0; i < NE; i++) {
        PIN[i].pa = (PHASE*)malloc(sizeof(PHASE) * (NS * 2));
        PIN[i].time = (float*)malloc(sizeof(float) * (NS * 2));
    }
    if ((fp1 = fopen(paDir, "r")) == NULL) {
        fprintf(stderr, "Unable to open phasefile\n");
        exit(-1);
    }
    while (fgets(lines, 150, fp1) != NULL) {
        sscanf(lines, "%s", flag);
        if (strcmp(flag, "#") == 0) {
            sscanf(lines, "%s %f %f %f %f %f %f %f %f %f %f %f %f %f %d", flag, &jk2,
                &jk3, &jk4, &jk5, &jk6, &jk7, &jk8, &jk9, &jk10, &jk11, &jk12,
                &jk13, &jk13, &ID_event);
            j = 0;
        } else {
            sscanf(lines, "%s %f %f %s", PIN[ID_event].pa[j].sta,
                &PIN[ID_event].time[j], &jk, PIN[ID_event].pa[j].phase);
            j++;
        }
    }
    fclose(fp1);
    // read observed event pairs (dt.ct)
    PO = (PAIR*)malloc(sizeof(PAIR) * (NP));
    for (i = 0; i < NP; i++) {
        PO[i].pk = (PICK*)malloc(sizeof(PICK) * (NS * 2));
    }
    if ((fp1 = fopen(dctDir, "r")) == NULL) {
        fprintf(stderr, "Unable to open dt.ct\n");
        exit(-1);
    }
    j = 0;
    k = 0;
    while (fgets(lines, 100, fp1) != NULL) {
        sscanf(lines, "%s", flag);
        if (strcmp(flag, "#") == 0) {
            sscanf(lines, "%s %d %d", flag, &PO[k].event1, &PO[k].event2);
            k++;
            if (k > NP) {
                fprintf(stderr, "Number of event-pairs exceeds the preset NP, please "
                                "recompile with a larger NP!\n");
                exit(-1);
            }
            kk = 0;
        } else {
            sscanf(lines, "%s %.2f %.2f %f %s", PO[k - 1].pk[kk].sta,
                &PO[k - 1].pk[kk].arr1, &PO[k - 1].pk[kk].arr2, &jk,
                PO[k - 1].pk[kk].phase);
            kk++;
            j++;
            if (kk > NS) {
                fprintf(stderr, "Number of stations exceeds the preset NS, please "
                                "recompile with a larger NS!\n");
                exit(-1);
            }
        }
    }
    fclose(fp1);
    np = k;
    printf("	FDTCC reads %d event-pairs\n", np);
    // read stations (REAL format)
    if ((fp1 = fopen(staDir, "r")) == NULL) {
        fprintf(stderr, "Unable to open stations\n");
        exit(-1);
    }
    ST = (STATION*)malloc(sizeof(STATION) * NS);
    k = 0;
    while (fgets(lines, 100, fp1) != NULL) {
        sscanf(lines, "%f %f %s %s %s %f", &ST[k].stlo, &ST[k].stla, ST[k].net,
            ST[k].sta, ST[k].comp, &ST[k].elev);
        k++;
        if (k > NS) {
            fprintf(stderr, "Number of stations exceeds the preset NS, please "
                            "recompile with a larger NS!\n");
            exit(-1);
        }
    }
    ns = k;
    fclose(fp1);
    printf("	FDTCC reads %d stations\n", k);
    // read event information (event.sel)
    if ((fp1 = fopen(eveDir, "r")) == NULL) {
        fprintf(stderr, "Unable to open event.sel\n");
        exit(-1);
    }
    EVE = (EVENT*)malloc(sizeof(EVENT) * NE);
    for (i = 0; i < NE; i++) {
        EVE[i].pSNR = (float*)malloc(sizeof(float) * ns);
        EVE[i].sSNR = (float*)malloc(sizeof(float) * ns);
    }
    k = 0;
    while (fgets(lines, 100, fp1) != NULL) {
        sscanf(lines, "%s %s %f %f %f %f %f %f %f %d", EVE[k].date, EVE[k].time,
            &EVE[k].evla, &EVE[k].evlo, &EVE[k].evdp, &jk, &jk1, &jk2, &jk3,
            &EVE[k].event);
        if (EVE[k].evdp > trh) {
            fprintf(stderr, "event out of the travel-time table, please update it\n");
            fprintf(stderr,
                "	maximum distance and depth are %f and %f, respectively\n", trx,
                trh);
            exit(-1);
        }
        k++;
        if (k > NE) {
            fprintf(stderr, "Number of events exceeds the preset NE, please "
                            "recompile with a larger NE!\n");
            exit(-1);
        }
    }
    ne = k;
    printf("	FDTCC reads %d events\n", k);
    Transfer_sec(EVE, ne);
    fclose(fp1);
    // read tt table (REAL format)
    if ((fp1 = fopen(tttDir, "r")) == NULL) {
        fprintf(stderr, "Unable to open tttDir\n");
        exit(-1);
    }
    TB = (TTT*)malloc(sizeof(TTT) * ntb);
    k = 0;
    while (fgets(lines, 300, fp1) != NULL) {
        sscanf(lines, "%lf %lf %lf %lf %lf %lf %lf %lf %s %s", &TB[k].gdist,
            &TB[k].dep, &TB[k].ptime, &TB[k].stime, &TB[k].prayp, &TB[k].srayp,
            &TB[k].phslow, &TB[k].shslow, TB[k].pphase, TB[k].sphase);
        k++;
        if (k > ntb) {
            fprintf(stderr, "Line of travel-time exceeds the preset ntb, please "
                            "recompile with a larger ntb!\n");
            exit(-1);
        }
    }
    printf("	FDTCC reads %d travel-times\n", k);

    // memory check
    npp = (int)(((wa + wb + 2)) / delta + 1);
    nss = (int)(((was + wbs + 2)) / delta + 1);
    size = ne * ns * 3;
    memory_require = ((size / (1024.0 * 1024.0 * 1024.0)) * npp) * sizeof(float);
    printf("	Memory require > %.2lf GB.\n", memory_require);

    // creat event pairs with theoretical travel times
    printf("	Creating database... \n");
    PT = (PAIR*)malloc(sizeof(PAIR) * np);
    for (i = 0; i < np; i++) {
        PT[i].pk = (PICK*)malloc(sizeof(PICK) * (ns * 2));
    }
    Cal_tt(PO, PT, EVE, ST);
    // update PT based on PO
#pragma omp parallel for shared(PT, PIN, np, ns) private(i, j)
    for (i = 0; i < np; i++) {
        for (j = 0; j < 2 * ns; j++) {
            Replace(PT, PIN, i, j);
        }
    }
#pragma omp barrier
    // read wavefrom
    fp1 = fopen(INPUT1, "w");
    fp2 = fopen(INPUT2, "w");
    fp3 = fopen(INPUT3, "w");
    if (fp1 == NULL || fp2 == NULL || fp3 == NULL) {
        fprintf(stderr, "Can't open INPUT files\n");
        exit(-1);
    }
    for (i = 0; i < ne; i++) {
        for (j = 0; j < ns; j++) {
            if (f == 0) {
                fprintf(fp1, "%s/%s/%s.%s.%c%cZ  %s      %.2lf   %d\n", wavDir,
                    EVE[i].date, ST[j].net, ST[j].sta, ST[j].comp[0], ST[j].comp[1],
                    ST[j].sta, EVE[i].sec, EVE[i].event);
                fprintf(fp2, "%s/%s/%s.%s.%c%cE  %s      %.2lf   %d\n", wavDir,
                    EVE[i].date, ST[j].net, ST[j].sta, ST[j].comp[0], ST[j].comp[1],
                    ST[j].sta, EVE[i].sec, EVE[i].event);
                fprintf(fp3, "%s/%s/%s.%s.%c%cN  %s     %.2lf   %d\n", wavDir,
                    EVE[i].date, ST[j].net, ST[j].sta, ST[j].comp[0], ST[j].comp[1],
                    ST[j].sta, EVE[i].sec, EVE[i].event);
            } else if (f == 1) {
                fprintf(fp1, "%s/%d/%s.%s.%c%cZ  %s     0.0   %d\n", wavDir,
                    EVE[i].event, ST[j].net, ST[j].sta, ST[j].comp[0], ST[j].comp[1],
                    ST[j].sta, EVE[i].event);
                fprintf(fp2, "%s/%d/%s.%s.%c%cE  %s     0.0   %d\n", wavDir,
                    EVE[i].event, ST[j].net, ST[j].sta, ST[j].comp[0], ST[j].comp[1],
                    ST[j].sta, EVE[i].event);
                fprintf(fp3, "%s/%d/%s.%s.%c%cN  %s     0.0   %d\n", wavDir,
                    EVE[i].event, ST[j].net, ST[j].sta, ST[j].comp[0], ST[j].comp[1],
                    ST[j].sta, EVE[i].event);
            }
        }
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    ptriger = (float*)malloc(ne * ns * sizeof(float));
    s1triger = (float*)malloc(ne * ns * sizeof(float));
    s2triger = (float*)malloc(ne * ns * sizeof(float));
    staP = (char**)malloc(sizeof(char*) * ne * ns);
    staS1 = (char**)malloc(sizeof(char*) * ne * ns);
    staS2 = (char**)malloc(sizeof(char*) * ne * ns);
    la_staP = (char**)malloc(sizeof(char*) * ne * ns);
    la_staS1 = (char**)malloc(sizeof(char*) * ne * ns);
    la_staS2 = (char**)malloc(sizeof(char*) * ne * ns);
    for (i = 0; i < ne * ns; i++) {
        staP[i] = (char*)malloc(sizeof(char) * 256);
        staS1[i] = (char*)malloc(sizeof(char) * 256);
        staS2[i] = (char*)malloc(sizeof(char) * 256);
        la_staP[i] = (char*)malloc(sizeof(char) * 10);
        la_staS1[i] = (char*)malloc(sizeof(char) * 10);
        la_staS2[i] = (char*)malloc(sizeof(char) * 10);
    }
    waveP = (float**)calloc(ne * ns, sizeof(float*));
    for (i = 0; i < ne * ns; i++)
        waveP[i] = (float*)calloc(npp, sizeof(float));
    waveS1 = (float**)calloc(ne * ns, sizeof(float*));
    for (i = 0; i < ne * ns; i++)
        waveS1[i] = (float*)calloc(nss, sizeof(float));
    waveS2 = (float**)calloc(ne * ns, sizeof(float*));
    for (i = 0; i < ne * ns; i++)
        waveS2[i] = (float*)calloc(nss, sizeof(float));
    labelP = (int*)malloc(ne * ns * sizeof(int));
    labelS1 = (int*)malloc(ne * ns * sizeof(int));
    labelS2 = (int*)malloc(ne * ns * sizeof(int));
    markP = (int*)malloc(ne * ns * sizeof(int));
    markS1 = (int*)malloc(ne * ns * sizeof(int));
    markS2 = (int*)malloc(ne * ns * sizeof(int));
    fp1 = fopen(INPUT1, "r");
    fp2 = fopen(INPUT2, "r");
    fp3 = fopen(INPUT3, "r");
    if (fp1 == NULL || fp2 == NULL || fp3 == NULL) {
        fprintf(stderr, "Can't open INPUT files\n");
        exit(-1);
    }
    for (i = 0; i < ne * ns; i++) {
        fscanf(fp1, "%s %s %f %d", staP[i], la_staP[i], &ptriger[i], &labelP[i]);
        fscanf(fp2, "%s %s %f %d", staS1[i], la_staS1[i], &s1triger[i],
            &labelS1[i]);
        fscanf(fp3, "%s %s %f %d", staS2[i], la_staS2[i], &s2triger[i],
            &labelS2[i]);
    }
    Correct_Pshift(PT, ptriger, la_staP, labelP);
    Correct_Sshift(PT, s1triger, la_staS1, labelS1);
    Correct_Sshift(PT, s2triger, la_staS2, labelS2);
#pragma omp parallel for shared(waveP, waveS1, waveS2, was, wbs, staP, ptriger, s1triger, s2triger, low, high, timezone, wb, wa, markP, markS1, markS2, np, ns) private(hd1, hd2, hd3, i, j)
    for (i = 0; i < ne * ns; i++) {
	markP[i] = 1;
        markS1[i] = 1;
        markS2[i] = 1;
        if ((waveP[i] = read_sac2(staP[i], &hd1, -3, ptriger[i] - timezone - wb - 1,
                 ptriger[i] - timezone + wa + 1))
            == NULL) {
            markP[i] = 0;
            //fprintf(stderr,"no station %s\n",staP[i]);
        } else if (low > 0 && high > 0) {
	    bpcc(waveP[i], hd1, low, high);
	}
	//in case user want to check waveform
	//char tmp[100]; 
	//sprintf(tmp,"%d/%s.%s.%c%cZ",EVE[i/ns].event,ST[i%ns].net, ST[i%ns].sta, ST[i%ns].comp[0],ST[i%ns].comp[1]);
	//printf("%d/%s.%s.%c%cZ\n",EVE[i/ns].event,ST[i%ns].net, ST[i%ns].sta, ST[i%ns].comp[0],ST[i%ns].comp[1]);
	//if(markP[i]==1)write_sac(tmp,hd1,waveP[i]);
        if ((waveS1[i] = read_sac2(staS1[i], &hd2, -3, s1triger[i] - timezone - wbs - 1,
                 s1triger[i] - timezone + was + 1))
            == NULL) {
            markS1[i] = 0;
            //fprintf(stderr,"no station %s\n",staS1[i]);
        } else if (low > 0 && high > 0) {
            bpcc(waveS1[i], hd2, low, high);
        }
        if ((waveS2[i] = read_sac2(staS2[i], &hd3, -3, s2triger[i] - wbs - timezone - 1,
                 s2triger[i] + was + timezone + 1))
            == NULL) {
            markS2[i] = 0;
            //fprintf(stderr,"no station %s\n",staS2[i]);
        } else if (low > 0 && high > 0) {
            bpcc(waveS2[i], hd3, low, high);
        }
    }
#pragma omp barrier
    Cal_sSNR(waveS1, waveS2, markS1);
    Cal_pSNR(waveP, markP);
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    // calculate cc
    printf("	FDTCC starts to calculate ccv\n");
    int pair_point[2];
#pragma omp parallel for shared(PT, EVE, waveP) private(i, j, pair_point)
    for (i = 0; i < np; i++) {
        Search_event(PT, EVE, pair_point, i);
        for (j = 0; j < ns; j++) {
            if (markP[pair_point[0] * ns + j] == 0 || markP[pair_point[1] * ns + j] == 0 || EVE[pair_point[0]].pSNR[j] <= thre_SNR || EVE[pair_point[1]].pSNR[j] <= thre_SNR) {
                PT[i].pk[2 * j].quality = 0;
                continue;
            }
            PT[i].pk[2 * j].quality = 1;
            SubccP(PT, waveP, pair_point, i, j);
        }
    }
#pragma omp barrier
#pragma omp parallel for shared(PT, EVE, waveS1, waveS2) private(i, j, \
    pair_point)
    for (i = 0; i < np; i++) {
        Search_event(PT, EVE, pair_point, i);
        for (j = 0; j < ns; j++) {
            if (markS1[pair_point[0] * ns + j] == 0 || markS1[pair_point[1] * ns + j] == 0 || EVE[pair_point[0]].sSNR[j] <= thre_SNR || EVE[pair_point[1]].sSNR[j] <= thre_SNR) {
                PT[i].pk[2 * j + 1].quality = 0;
                continue;
            }
            PT[i].pk[2 * j + 1].quality = 1;

            SubccS(PT, waveS1, waveS2, pair_point, i, j);
        }
    }
#pragma omp barrier
    // output dt.cc
    fp1 = fopen(OUTPUT, "w");
    for (i = 0; i < np; i++) {
        fprintf(fp1, "#	%d	%d	0\n", PT[i].event1, PT[i].event2);
        for (j = 0; j < 2 * ns; j++) {
            if (PT[i].pk[j].quality == 1 && PT[i].pk[j].ccv >= threshold && PT[i].pk[j].ccv > 0) {
                fprintf(fp1, "%5s %10.4f %10.2lf %3s\n", PT[i].pk[j].sta,
                    PT[i].pk[j].arr1 - PT[i].pk[j].arr2 + PT[i].pk[j].shift,
                    PT[i].pk[j].ccv, PT[i].pk[j].phase);
            }
        }
    }
    fclose(fp1);
    printf("	Results were written in dt.cc\n");
    // free memory
    for (i = 0; i < NP; i++) {
        free(PO[i].pk);
    }
    for (i = 0; i < np; i++) {
        free(PT[i].pk);
    }
    free(PT);
    free(PO);
    free(ST);
    free(EVE);
    free(TB);
    for (i = 0; i < ne * ns; i++) {
        free(staP[i]);
        free(staS1[i]);
        free(staS2[i]);
        free(waveP[i]);
        free(waveS1[i]);
        free(waveS2[i]);
        free(la_staP[i]);
        free(la_staS1[i]);
        free(la_staS2[i]);
    }
    free(la_staP);
    free(la_staS1);
    free(la_staS2);
    free(staP);
    free(staS1);
    free(staS2);
    free(waveP);
    free(waveS1);
    free(waveS2);
    free(ptriger);
    free(s1triger);
    free(s2triger);
    free(labelP);
    free(labelS1);
    free(labelS2);
    for (i = 0; i < NE; i++) {
        free(PIN[i].pa);
        free(PIN[i].time);
    }
    free(PIN);
}

// transfer hms to sec
void Transfer_sec(EVENT* EVE, int ne)
{
    int i;
    int h, m, s, ms;
    for (i = 0; i < ne; i++) {
        switch (strlen(EVE[i].time)) {
        case 1:
            sscanf(EVE[i].time, "%1d", &ms);
            EVE[i].sec = ((double)((ms))) / 100.0;
            break;
        case 2:
            sscanf(EVE[i].time, "%2d", &ms);
            EVE[i].sec = ((double)(ms)) / 100.0;
            break;
        case 3:
            sscanf(EVE[i].time, "%1d%2d", &s, &ms);
            EVE[i].sec = ((double)((s)*100 + ms)) / 100.0;
            break;
        case 4:
            sscanf(EVE[i].time, "%2d%2d", &s, &ms);
            EVE[i].sec = ((double)((s)*100 + ms)) / 100;
            break;
        case 5:
            sscanf(EVE[i].time, "%1d%2d%2d", &m, &s, &ms);
            EVE[i].sec = ((double)((m * 60 + s) * 100 + ms)) / 100;
            break;
        case 6:
            sscanf(EVE[i].time, "%2d%2d%2d", &m, &s, &ms);
            EVE[i].sec = ((double)((m * 60 + s) * 100 + ms)) / 100;
            break;
        case 7:
            sscanf(EVE[i].time, "%1d%2d%2d%2d", &h, &m, &s, &ms);
            EVE[i].sec = ((double)((h * 3600 + m * 60 + s) * 100 + ms)) / 100;
            break;
        case 8:
            sscanf(EVE[i].time, "%2d%2d%2d%2d", &h, &m, &s, &ms);
            EVE[i].sec = ((double)((h * 3600 + m * 60 + s) * 100 + ms)) / 100;
            break;
        default:
            printf("Wrong time point\n");
            break;
        }
    }
}

void Cal_tt(PAIR* PO, PAIR* PT, EVENT* EVE, STATION* ST)
{
    int i, j, ih, ig, k;
    extern float trx, tdx, tdh;
    extern int np, ns;
    int event[2];
    extern TTT* TB;
    double GCarc1, GCarc2;
    #pragma omp parallel for shared(PO, PT, EVE, ST, TB, np, ns, trx, tdx, tdh) private(i, j, ih, ig, k, event, GCarc1, GCarc2)
    for (i = 0; i < np; i++) {
        Search_event(PO, EVE, event, i);
        PT[i].event1 = PO[i].event1;
        PT[i].event2 = PO[i].event2;
        k = 0;
        for (j = 0; j < ns; j++) {
            strcpy(PT[i].pk[k].sta, ST[j].sta);
            strcpy(PT[i].pk[k + 1].sta, ST[j].sta);
            strcpy(PT[i].pk[k].phase, "P");
            strcpy(PT[i].pk[k + 1].phase, "S");
            ddistaz(ST[j].stla, ST[j].stlo, EVE[event[0]].evla, EVE[event[0]].evlo, &GCarc1);
            ih = rint(EVE[event[0]].evdp / tdh);
            ig = ih * rint(trx / tdx) + rint(GCarc1 / tdx);
	    PT[i].pk[k].arr1 = TB[ig].ptime + (GCarc1 - TB[ig].gdist) * TB[ig].prayp + (EVE[event[0]].evdp - TB[ig].dep) * TB[ig].phslow;
            PT[i].pk[k + 1].arr1 = TB[ig].stime + (GCarc1 - TB[ig].gdist) * TB[ig].srayp + (EVE[event[0]].evdp - TB[ig].dep) * TB[ig].shslow;
            ddistaz(ST[j].stla, ST[j].stlo, EVE[event[1]].evla, EVE[event[1]].evlo, &GCarc2);
            ih = rint(EVE[event[1]].evdp / tdh);
            ig = ih * rint(trx / tdx) + rint(GCarc2 / tdx);
            PT[i].pk[k].arr2 = TB[ig].ptime + (GCarc2 - TB[ig].gdist) * TB[ig].prayp + (EVE[event[1]].evdp - TB[ig].dep) * TB[ig].phslow;
            PT[i].pk[k + 1].arr2 = TB[ig].stime + (GCarc2 - TB[ig].gdist) * TB[ig].srayp + (EVE[event[1]].evdp - TB[ig].dep) * TB[ig].shslow;
            PT[i].pk[k].diff = PT[i].pk[k].arr2 - PT[i].pk[k].arr1;
            PT[i].pk[k + 1].diff = PT[i].pk[k + 1].arr2 - PT[i].pk[k + 1].arr1;
	    k = k + 2;
        }
    }
    #pragma omp barrier
}

void Search_event(PAIR* PO, EVENT* EVE, int* serial, int n)
{
    int i;
    extern int ne;
    for (i = 0; i < ne; i++) {
        if (EVE[i].event == PO[n].event1) {
            serial[0] = i;
            break;
        }
    }
    for (i = 0; i < ne; i++) {
        if (EVE[i].event == PO[n].event2) {
            serial[1] = i;
            break;
        }
    }
}

/* * Modified by M. Zhang
c Subroutine to calculate the Great Circle Arc distance
c    between two sets of geographic coordinates
c
c Given:  stalat => Latitude of first point (+N, -S) in degrees
c         stalon => Longitude of first point (+E, -W) in degrees
c         evtlat => Latitude of second point
c         evtlon => Longitude of second point
c
c Returns:  delta => Great Circle Arc distance in degrees
c           az    => Azimuth from pt. 1 to pt. 2 in degrees
c           baz   => Back Azimuth from pt. 2 to pt. 1 in degrees
c
c If you are calculating station-epicenter pairs, pt. 1 is the station
c
c Equations take from Bullen, pages 154, 155
c
c T. Owens, September 19, 1991
c           Sept. 25 -- fixed az and baz calculations
c
  P. Crotwell, Setember 27, 1994
            Converted to c to fix annoying problem of fortran giving wrong
               answers if the input doesn't contain a decimal point.
*/
void ddistaz(double stalat, double stalon, double evtlat, double evtlon,
    double* delta)
{
    // double stalat, stalon, evtlat, evtlon;
    // double delta, az, baz;
    double scolat, slon, ecolat, elon;
    double a, b, c, d, e, aa, bb, cc, dd, ee, g, gg, h, hh, k, kk;
    double rhs1, rhs2, sph, rad, del, daz, az, dbaz, pi, piby2;
    /*
     stalat = atof(argv[1]);
     stalon = atof(argv[2]);
     evtlat = atof(argv[3]);
     evtlon = atof(argv[4]);
  */
    pi = 3.141592654;
    piby2 = pi / 2.0;
    rad = 2. * pi / 360.0;
    sph = 1.0 / 298.257;

    scolat = piby2 - atan((1. - sph) * (1. - sph) * tan(stalat * rad));
    ecolat = piby2 - atan((1. - sph) * (1. - sph) * tan(evtlat * rad));
    slon = stalon * rad;
    elon = evtlon * rad;
    a = sin(scolat) * cos(slon);
    b = sin(scolat) * sin(slon);
    c = cos(scolat);
    d = sin(slon);
    e = -cos(slon);
    g = -c * e;
    h = c * d;
    k = -sin(scolat);
    aa = sin(ecolat) * cos(elon);
    bb = sin(ecolat) * sin(elon);
    cc = cos(ecolat);
    dd = sin(elon);
    ee = -cos(elon);
    gg = -cc * ee;
    hh = cc * dd;
    kk = -sin(ecolat);
    del = acos(a * aa + b * bb + c * cc);
    *delta = del / rad; // delta
}

void Correct_Pshift(PAIR* PT, float* a, char** b, int* c)
{
    int i, j, k;
    #pragma omp parallel for shared(PT, a, b, c, ne, ns, np) private(i, j, k)
    for (i = 0; i < ne * ns; i++) {
        for (j = 0; j < np; j++) {
            if (PT[j].event1 == c[i]) {
                for (k = 0; k < 2 * ns; k++) {
                    if (strcmp(PT[j].pk[k].sta, b[i]) == 0 && strcmp(PT[j].pk[k].phase, "P") == 0) {
                        a[i] = a[i] + PT[j].pk[k].arr1;
			break;
                    }
                }
                break;
            }
            if (PT[j].event2 == c[i]) {
                for (k = 0; k < 2 * ns; k++) {
                    if (strcmp(PT[j].pk[k].sta, b[i]) == 0 && strcmp(PT[j].pk[k].phase, "P") == 0) {
                        a[i] = a[i] + PT[j].pk[k].arr2;
                        break;
                    }
                }
                break;
            }
        }
    }
    #pragma omp barrier
}

void Correct_Sshift(PAIR* PT, float* a, char** b, int* c)
{
    int i, j, k;
    #pragma omp parallel for shared(PT, a, b, c, ne, ns, np) private(i, j, k)
    for (i = 0; i < ne * ns; i++) {
        for (j = 0; j < np; j++) {
            if (PT[j].event1 == c[i]) {
                for (k = 0; k < 2 * ns; k++) {
                    if (strcmp(PT[j].pk[k].sta, b[i]) == 0 && strcmp(PT[j].pk[k].phase, "S") == 0) {
                        a[i] = a[i] + PT[j].pk[k].arr1;
                        break;
                    }
                }
                break;
            }
            if (PT[j].event2 == c[i]) {
                for (k = 0; k < 2 * ns; k++) {
                    if (strcmp(PT[j].pk[k].sta, b[i]) == 0 && strcmp(PT[j].pk[k].phase, "S") == 0) {
                        a[i] = a[i] + PT[j].pk[k].arr2;
                        break;
                    }
                }
                break;
            }
        }
    }
    #pragma omp barrier
}

void SubccP(PAIR* PT, float** waveP, int* a, int i, int j)
{
    int k, kk, Npoint, Wpoint, ref_shift;
    float s_p;
    double cc, norm, normMaster, tmp;
    extern int ns;
    extern float delta, wa, wb, wf, thre_shift;
    float w, wa1;
    int t_shift = (int)(1/delta);
    s_p = PT[i].pk[2 * j + 1].arr1 - PT[i].pk[2 * j].arr1;
    if (s_p <= 0)
        PT[i].pk[2 * j].quality = 0;
    else {
        wa1 = wa;
        w = wf;
        if (wa > 0.9 * s_p)
            wa1 = 0.9 * s_p;
        if (w > 0.5 * (wa1 + wb))
            w = 0.5 * (wa1 + wb);
        ref_shift = (int)(w / delta);
        w = ref_shift * delta;
        Npoint = (int)(2 * w / delta - 0.5);
        Wpoint = (int)((wa1 + wb) / delta - 0.5);
        PT[i].pk[2 * j].ccv = 0;
        PT[i].pk[2 * j].shift = 0;
        normMaster = 0.0;
        norm = 0.0;
        for (k = 0; k <= Wpoint; k++) {
            norm += waveP[a[1] * ns + j][k + t_shift] * waveP[a[1] * ns + j][k + t_shift];
            normMaster += waveP[a[0] * ns + j][k + t_shift] * waveP[a[0] * ns + j][k + t_shift];
        }
        for (k = 0; k <= Npoint; k++) {
            cc = 0.0;
            if (k <= ref_shift) {
                for (kk = ref_shift - k; kk <= Wpoint; kk++) {
                    cc += waveP[a[0] * ns + j][kk - ref_shift + k + t_shift] * waveP[a[1] * ns + j][kk + t_shift];
                }
            } else {
                for (kk = 0; kk <= Wpoint - (k - ref_shift); kk++) {
                    cc += waveP[a[0] * ns + j][kk + k - ref_shift + t_shift] * waveP[a[1] * ns + j][kk + t_shift];
                }
            }
            tmp = cc / (sqrt(norm) * sqrt(normMaster));
            if (fabs(tmp) > PT[i].pk[2 * j].ccv) {
                PT[i].pk[2 * j].ccv = fabs(tmp);
                PT[i].pk[2 * j].shift = k * delta - w;
            }
        }
        if (fabs(PT[i].pk[2 * j].arr1 - PT[i].pk[2 * j].arr2 + PT[i].pk[2 * j].shift) > thre_shift) {
            PT[i].pk[2 * j].quality = 0;
        }
    }
}

void SubccS(PAIR* PT, float** waveS1, float** waveS2, int* a, int i, int j)
{
    int Npoint, Wpoint, k, kk, ref_shift;
    float s_p;
    double cc1, norm1, normMaster1, cc2, norm2, normMaster2, tmp;
    extern int ns;
    extern float wbs, was, delta, wfs, thre_shift;
    int tt;
    int t_shift = (int)(1/delta);
    float w, wbs1;
    s_p = PT[i].pk[2 * j + 1].arr1 - PT[i].pk[2 * j].arr1;
    if (s_p <= 0)
        PT[i].pk[2 * j + 1].quality = 0;
    else {
        wbs1 = wbs;
        w = wfs;
        tt = 0;
        if (wbs1 > 0.5 * s_p) {
            wbs1 = 0.5 * s_p;
            tt = (int)((wbs - 0.5 * s_p) / delta - 0.5);
        }
        if (w > 0.5 * (was + wbs1)) {
            w = 0.5 * (was + wbs1);
        }
        ref_shift = (int)(w / delta);
        w = ref_shift * delta;
        Npoint = (int)(2 * w / delta - 0.5);
        Wpoint = (int)((was + wbs1) / delta - 0.5);
        PT[i].pk[2 * j + 1].ccv = 0;
        PT[i].pk[2 * j + 1].shift = 0;
        normMaster1 = 0.0;
        norm1 = 0.0;
        normMaster2 = 0.0;
        norm2 = 0.0;
        for (k = 0; k <= Wpoint; k++) {
            norm1 += waveS1[a[1] * ns + j][tt + k +t_shift] * waveS1[a[1] * ns + j][tt + k +t_shift];
            normMaster1 += waveS1[a[0] * ns + j][tt + k + t_shift] * waveS1[a[0] * ns + j][tt + k + t_shift];
            norm2 += waveS2[a[1] * ns + j][tt + k + t_shift] * waveS2[a[1] * ns + j][k + tt + t_shift];
            normMaster2 += waveS2[a[0] * ns + j][tt + k + t_shift] * waveS2[a[0] * ns + j][tt + k + t_shift];
        }
        for (k = 0; k <= Npoint; k++) {
            cc1 = 0.0;
            cc2 = 0.0;
            if (k <= ref_shift) {
                for (kk = ref_shift - k; kk <= Wpoint; kk++) {
                    cc1 += waveS1[a[0] * ns + j][tt + kk - ref_shift + k + t_shift] * waveS1[a[1] * ns + j][tt + kk + t_shift];
                    cc2 += waveS2[a[0] * ns + j][tt + kk - ref_shift + k + t_shift] * waveS2[a[1] * ns + j][tt + kk + t_shift];
                }
            } else {
                for (kk = 0; kk <= Wpoint - (k - ref_shift); kk++) {
                    cc1 += waveS1[a[0] * ns + j][tt + kk - ref_shift + k + t_shift] * waveS1[a[1] * ns + j][tt + kk + t_shift];
                    cc2 += waveS2[a[0] * ns + j][tt + kk - ref_shift + k + t_shift] * waveS2[a[1] * ns + j][tt + kk + t_shift];
                }
            }
            tmp = ((cc1 / (sqrt(norm1) * sqrt(normMaster1))) + (cc2 / (sqrt(norm2) * sqrt(normMaster2)))) / 2;
            if (fabs(tmp) > PT[i].pk[2 * j + 1].ccv) {
                PT[i].pk[2 * j + 1].ccv = fabs(tmp);
                PT[i].pk[2 * j + 1].shift = k * delta - w;
            }
        }
        if (fabs(PT[i].pk[2 * j + 1].arr1 - PT[i].pk[2 * j + 1].arr2 + PT[i].pk[2 * j + 1].shift) > thre_shift)
            PT[i].pk[2 * j + 1].quality = 0;
    }
}

void Cal_pSNR(float** wave, int* mark)
{
    extern float delta, wa, wb;
    extern int ns, ne;
    int i, j, k;
    double s, n;
    int t_shift = (int)(1/delta);
    int spoint, npoint;
    extern EVENT* EVE;
    spoint = (int)(wa / delta - 0.5);
    npoint = (int)(wb / delta - 0.5);
    #pragma omp parallel for shared(EVE, wave, mark, t_shift, spoint, npoint, ne, ns) private(i, j, k, s, n)
    for (i = 0; i < ne; i++) {
        for (k = 0; k < ns; k++) {
            if (mark[k + i * ns] == 0) {
                EVE[i].pSNR[k] = 0;
                continue;
            }
            s = 0;
            n = 0;
            for (j = 0; j < spoint; j++)
                s += wave[i * ns + k][j + npoint + t_shift] * wave[i * ns + k][j + npoint + t_shift];
            for (j = 0; j < npoint; j++)
                n += wave[i * ns + k][j + t_shift] * wave[i * ns + k][j + t_shift];
            EVE[i].pSNR[k] = (s / spoint) / (n / npoint);
        }
    }
    #pragma omp barrier
}

void Cal_sSNR(float** wave1, float** wave2, int* mark)
{
    extern float delta, was, wbs;
    extern int ns, ne;
    int i, j, k;
    int t_shift = (int)(1/delta);
    double s1, n1, s2, n2;
    int spoint, npoint;
    extern EVENT* EVE;
    spoint = (int)(was / delta - 0.5);
    npoint = (int)(wbs / delta - 0.5);
    #pragma omp parallel for shared(EVE, wave1, wave2, mark, t_shift, spoint, npoint, ne, ns) private(i, j, k, s1, s2, n1, n2)
    for (i = 0; i < ne; i++) {
        for (k = 0; k < ns; k++) {
            if (mark[i * ns + k] == 0) {
                EVE[i].sSNR[k] = 0;
                continue;
            }
            s1 = 0;
            n1 = 0;
            s2 = 0;
            n2 = 0;
            for (j = 0; j < spoint; j++)
                s1 += wave1[i * ns + k][j + npoint + t_shift] * wave1[i * ns + k][j + npoint + t_shift];
            for (j = 0; j < npoint; j++)
                n1 += wave1[i * ns + k][j + t_shift] * wave1[i * ns + k][j + t_shift];
            for (j = 0; j < spoint; j++)
                s2 += wave2[i * ns + k][j + npoint + t_shift] * wave2[i * ns + k][j + npoint + t_shift];
            for (j = 0; j < npoint; j++)
                n2 += wave2[i * ns + k][j + t_shift] * wave2[i * ns + k][j + t_shift];
            EVE[i].sSNR[k] = 0.5 * ((s1 / spoint) / (n1 / npoint) + (s2 / spoint) / (n2 / npoint));
        }
    }
    #pragma omp barrier
}

void Replace(PAIR* PT, PHASEINF* PIN, int a, int b)
{
    int i, j;
    extern int np, ns;
    for (i = 0; i < 2 * ns; i++) {
        if (strcmp(PT[a].pk[b].sta, PIN[PT[a].event1].pa[i].sta) == 0 && strcmp(PT[a].pk[b].phase, PIN[PT[a].event1].pa[i].phase) == 0) {
            PT[a].pk[b].arr1 = PIN[PT[a].event1].time[i];
            break;
        }
    }
    for (i = 0; i < 2 * ns; i++) {
        if (strcmp(PT[a].pk[b].sta, PIN[PT[a].event2].pa[i].sta) == 0 && strcmp(PT[a].pk[b].phase, PIN[PT[a].event2].pa[i].phase) == 0) {
            PT[a].pk[b].arr2 = PIN[PT[a].event2].time[i];
            break;
        }
    }
}

//hanning taper before bandpass filter
void taper(float* yarray, int nlen, float start, float end)
{
    float ang, cs;
    int m1, m2, m3, m4, m5;
    int i, j, k, xi;

    m1 = (int)(nlen * start + 0.5);
    m2 = m1 + 1;

    ang = 3.1415926 / (float)(m1);

    for (i = 0; i <= m1; i++) {
        xi = i;
        cs = (1 - cos(xi * ang)) / 2.0;
        yarray[i] = yarray[i] * cs;
    }

    m3 = (int)(nlen * end + 0.5);
    m5 = nlen - m3 - 1;
    m4 = m5 + 1;
    ang = 3.1415926 / (float)(m3);

    for (k = m2; k <= m5; k++) {
        yarray[k] = yarray[k];
    }
    for (j = m4; j < nlen; j++) {
        xi = j + 1 - nlen;
        cs = (1 - cos(xi * ang)) / 2.0;
        yarray[j] = yarray[j] * cs;
    }
}

/* remove trend a*i + b */
void rtrend(float* y, int n)
{
    int i;
    double a, b, a11, a12, a22, y1, y2;
    y1 = y2 = 0.;
    for (i = 0; i < n; i++) {
        y1 += i * y[i];
        y2 += y[i];
    }
    a12 = 0.5 * n * (n - 1);
    a11 = a12 * (2 * n - 1) / 3.;
    a22 = n;
    b = a11 * a22 - a12 * a12;
    a = (a22 * y1 - a12 * y2) / b;
    b = (a11 * y2 - a12 * y1) / b;
    for (i = 0; i < n; i++) {
        y[i] = y[i] - a * i - b;
    }
}

//do bandpass filtering for templates and traces
void bpcc(float* yarray, SACHEAD hd, float low, float high)
{
    /* Local variables */
    //float low, high;
    double attenuation, transition_bandwidth;
    int nlen;
    //SACHEAD hd;
    double delta_d;
    int order;
    //float *yarray;
    int passes;
    float total, sum, mean, taperb;
    int j;

    sum = 0.0;
    //rmean
    for (j = 0; j < hd.npts; j++) {
        sum += yarray[j];
    }
    for (j = 0; j < hd.npts; j++) {
        yarray[j] = yarray[j] - sum / hd.npts;
    }
    delta_d = hd.delta;
    nlen = hd.npts;

    //rtrend
    rtrend(yarray, nlen);
    /*taper function*/
    //taper(yarray,hd.npts,0.05,0.05); //sac default hanning window 0.05
    taperb = 0.0001;
    if (hd.npts < 20000) {
        taperb = 0.01;
    }
    taper(yarray, hd.npts, taperb, taperb);
    passes = 2;
    order = 4;
    transition_bandwidth = 0.0;
    attenuation = 0.0;
    xapiir(yarray, nlen, SAC_BUTTERWORTH, transition_bandwidth, attenuation, order, SAC_BANDPASS, low, high, delta_d, passes);
    //    write_sac("test.sac",hd,yarray);
    //    exit(-1);
}
