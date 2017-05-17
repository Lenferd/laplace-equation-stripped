#include <iostream>
#include <cmath>
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"

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


    tbb::task_scheduler_init init(1);

    FILE *infile = fopen("../../../initial/settings.ini", "r");

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

    string filename = "../../../result/result_tbb_original.txt";
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
                          }, tbb::simple_partitioner());

        ++stepCounter;

        for (int l = 0; l < settings.vectSize; ++l) {
            if (loc_Change[l] > glob_change) {
                glob_change = loc_Change[l];
            }
        }
//        std::cout << glob_change << std::endl;
    } while ( glob_change > settings.epsilon);

    t_finish = tbb::tick_count::now();

    string filename2 = "../../../result/result_tbb.txt";
    FILE *fp2;
    fp2 = fopen(filename2.c_str(), "w");

    for (int i = 0; i < settings.dim; i++) {
        for (int j = 0; j < settings.dim; j++)
            fprintf(fp2, "%.15le ", vect[i * settings.dim + j]);
        fprintf(fp2, "\n");
    }

    fclose(fp2);

//    printf("Proc count:\t %d\n", );
    printf("Epsilon:\t %lf\n", settings.epsilon);
    printf("Dim size:\t %d\n", settings.dim);
    printf("Step calc:\t %d\n", stepCounter);
    printf("Run time:\t %.15lf\n", (t_finish-t_start).seconds());


    return 0;
}