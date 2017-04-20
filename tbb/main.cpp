#include <iostream>
#include <cmath>
#include <tbb/tbb.h>

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
    return sin(M_PI * y);
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

    const int NUM_THREAD = 2;

    Settings settings;
    settings.dim = 10000 + 2;  // 2 - boundaries
    settings.vectSize = settings.dim * settings.dim;
    settings.epsilon = 1e-4;

    settings.xStart = 0;
    settings.xEnd = 1;

    settings.yStart = 0;
    settings.yEnd= 1;


    tbb::tick_count t_start, t_finish;
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

    int stepCounter = 0;

    double globChange;
    double locChange;

    double tempPrevVal;
    double tempChange;

    double* loc_Change = new double[settings.vectSize];

    double glob_change;

    t_start = tbb::tick_count::now();

    do {
        glob_change = 0;
        for (int k = 0; k < settings.vectSize; ++k) {
            loc_Change[k] = 0;
        }
        tbb::parallel_for(tbb::blocked_range<int>(1, settings.dim - 1, 1),
                          [vect, loc_Change, settings](const tbb::blocked_range<int>& r) {
                              for (int j = r.begin(); j < r.end(); ++j) {
                                  for (int i = 1; i < settings.dim - 1; ++i) {
                                      double tempPrevVal = vect[i * settings.dim + j];
                                      vect[i * settings.dim + j] = 0.25 * (vect[i * settings.dim + j + 1] + vect[i * settings.dim + j - 1] +
                                                                                vect[(i + 1) * settings.dim + j] + vect[(i - 1) * settings.dim + j]);
                                      double tempChange = fabs(vect[settings.dim * i + j] - tempPrevVal);
                                      if (loc_Change[j] < tempChange) {
                                          loc_Change[j] = tempChange;
                                      }
                                  }
                              }
                          }, tbb::static_partitioner());

        ++stepCounter;

        for (int l = 0; l < settings.vectSize; ++l) {
            if (loc_Change[l] > glob_change) {
                glob_change = loc_Change[l];
            }
        }
        std::cout << glob_change << std::endl;
    } while ( glob_change > settings.epsilon);

    t_finish = tbb::tick_count::now();


    printf("Proc count:\t %d\n", NUM_THREAD);
    printf("Epsilon:\t %lf\n", settings.epsilon);
    printf("Dim size:\t %d\n", settings.dim);
    printf("Step calc:\t %d\n", stepCounter);
    printf("Run time:\t %.15lf\n", (t_finish-t_start).seconds());


    return 0;
}