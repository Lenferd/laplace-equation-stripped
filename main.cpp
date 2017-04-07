#include <iostream>
#include "omp.h"
#include <cmath>

struct Settings {
    int dim;
    double epsilon;
    double yStartTime;
    double yEndTime;
    double xStartTime;
    double xEndTime;

};

double funct1(double x, double y) {

}

int main(int argc, char** argv) {

    const int NUM_THREAD = 4;

    Settings settings;
    settings.dim = 1000 + 2;  // 2 - boundaries
    settings.epsilon = 1e-3;

    settings.xStartTime = 0;
    settings.xEndTime = 1;

    settings.yStartTime = 0;
    settings.yEndTime = 1;


    double time_S, time_E;
    double *vect;
    vect = new double[settings.dim * settings.dim];

    for (int k = 0; k < settings.dim * settings.dim; ++k) {
        vect[k] = 0;
    }

    // File vector?
    for (int i = 1; i < settings.dim - 1; ++i) {
        for (int j = 1; j < settings.dim - 1; ++j) {
            vect[settings.dim * i + j] = 0;
        }
    }

    int stepCounter = 0;

    double globChange;
    double locChange;

    double tempPrevVal;
    double tempChange;

    omp_set_num_threads(NUM_THREAD);
    omp_lock_t globChangeLock;
    omp_init_lock(&globChangeLock);

    time_S = omp_get_wtime();

    do {
        globChange = 0;

        #pragma omp parallel for shared(vect, settings, globChange) private(tempPrevVal, tempChange, locChange)
        for (int j = 1; j < settings.dim - 1; ++j) {
            locChange = 0;
            for (int i = 1; i < settings.dim - 1; ++i) {

                tempPrevVal = vect[i * settings.dim + j];
                vect[i * settings.dim + j] = 0.25 * (vect[i * settings.dim + j + 1] + vect[i * settings.dim + j - 1] +
                                                    vect[(i + 1) * settings.dim + j + 1] + vect[(i + 1) * settings.dim + j - 1]);
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
        }
        ++stepCounter;
    } while ( globChange > settings.epsilon);

    time_E = omp_get_wtime();

    printf("Run time:\t %.15lf\n", time_E-time_S);
}