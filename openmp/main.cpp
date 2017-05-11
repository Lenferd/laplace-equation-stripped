#include <iostream>
#include <fstream>
#include "omp.h"
#include <cmath>

using std::string;

struct Settings {
    int dim;
    double epsilon;
    double yStart;
    double yEnd;
    double xStart;
    double xEnd;
    int vectSize;

};

double fy1(double y) {
    return exp(sin(M_PI * y));
}

double fy2(double y) {
    return sin(M_PI * y) + y;
}

double fx1(double x) {
    return 0;
}

double fx2(double x) {
    return 0;
}

int main(int argc, char** argv) {

    const int NUM_THREAD = 1;


    FILE *infile = fopen("../../initial/settings.ini", "r");

    if (infile == NULL) {
        std::cout << "File open error" << std::endl;
        exit(-1);
    }
    Settings settings;

    fscanf(infile, "DIM=%d\n", &settings.dim);
    fscanf(infile, "EPS=%lf\n", &settings.epsilon);
    fscanf(infile, "XSTART=%lf\n", &settings.xStart);
    fscanf(infile, "XEND=%lf\n", &settings.xEnd);
    fscanf(infile, "YSTART=%lf\n", &settings.yStart);
    fscanf(infile, "YEND=%lf\n", &settings.yEnd);
    settings.vectSize = settings.dim * settings.dim;


    double time_S, time_E;
    double *vect;
    vect = new double[settings.vectSize];

    for (int k = 0; k < settings.vectSize; ++k) {
        vect[k] = 0;
    }

    // boundaries fill
    double h = fabs(settings.xEnd - settings.xStart)/ (settings.vectSize - 1);
    double xPos = settings.xStart;
    double yPos = settings.yStart;
    for (int i = 0; i < settings.dim; i++) {
        vect[i * settings.dim] = fx1(xPos);
        vect[i] = fy1(yPos);

        vect[i * settings.dim + settings.dim - 1] = fx2(xPos);
        vect[(settings.dim - 1) * settings.dim + i] = fy2(yPos);

        xPos += h;
        yPos += h;
    }


    string filename = "../../result/result_openmp_original.txt";
    FILE *fp;
    fp = fopen(filename.c_str(), "w");
//
    for (int i = 0; i < settings.dim; i++) {
        for (int j = 0; j < settings.dim; j++)
            fprintf(fp, "%.15le ", vect[i*(settings.dim) + j]);
        fprintf(fp, "\n");
    }

    int stepCounter = 0;

    double globChange;
    double locChange;

    double tempPrevVal;
    double tempChange;

    omp_set_num_threads(NUM_THREAD);
    omp_lock_t globChangeLock;
    omp_init_lock(&globChangeLock);

    omp_lock_t *rowLock = new omp_lock_t[settings.dim];
    for (int l = 0; l < settings.dim; ++l) {
        omp_init_lock(&rowLock[l]);
    }
    time_S = omp_get_wtime();

    do {
        globChange = 0;
        #pragma omp parallel for shared(vect, settings, globChange) private(tempPrevVal, tempChange, locChange)
        for (int j = 1; j < settings.dim - 1; ++j) {
            locChange = 0;

            omp_set_lock(&rowLock[j+1]);
            omp_set_lock(&rowLock[j]);
            omp_set_lock(&rowLock[j-1]);

            for (int i = 1; i < settings.dim - 1; ++i) {

                tempPrevVal = vect[i * settings.dim + j];
                vect[i * settings.dim + j] = 0.25 * (vect[i * settings.dim + j + 1] + vect[i * settings.dim + j - 1] +
                                                    vect[(i + 1) * settings.dim + j] + vect[(i - 1) * settings.dim + j]);
                tempChange = fabs(vect[settings.dim * i + j] - tempPrevVal);
                if (locChange < tempChange) {
                    locChange = tempChange;
                }
            }

            omp_set_lock(&globChangeLock);
            if (globChange < locChange) {
                globChange = locChange;
            }
            omp_unset_lock(&globChangeLock);

            omp_unset_lock(&rowLock[j-1]);
            omp_unset_lock(&rowLock[j]);
            omp_unset_lock(&rowLock[j+1]);
        }
        ++stepCounter;
    } while ( globChange > settings.epsilon);


    time_E = omp_get_wtime();

    string filename2 = "../../result/result_openmp.txt";
    FILE *fp2;
    fp2 = fopen(filename2.c_str(), "w");
//
    for (int i = 0; i < settings.dim; i++) {
        for (int j = 0; j < settings.dim; j++)
            fprintf(fp2, "%.15le ", vect[i * settings.dim + j]);
        fprintf(fp2, "\n");
    }

    fclose(fp2);

    printf("Proc count:\t %d\n", NUM_THREAD);
    printf("Epsilon:\t %lf\n", settings.epsilon);
    printf("Dim size:\t %d\n", settings.dim);
    printf("Step calc:\t %d\n", stepCounter);
    printf("Run time:\t %.15lf\n", time_E-time_S);


    return 0;
}