#include <iostream>
#include <fstream>
#include <cmath>
#include "mpi.h"
#include "omp.h"

const int ROOT = 0;
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

int main(int argc, char **argv) {

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
    settings.vectSize = settings.dim * settings.dim;    // with +2 boundaries


    double time_S, time_E;
    double *vect;
    vect = new double[settings.vectSize];

    for (int k = 0; k < settings.vectSize; ++k) {
        vect[k] = 0;
    }

    // boundaries fill
    double h = fabs(settings.xEnd - settings.xStart) / (settings.vectSize - 1);
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


    string filename = "../../result/result_mpi_original.txt";
    FILE *fp;
    fp = fopen(filename.c_str(), "w");
//
    for (int i = 0; i < settings.dim; i++) {
        for (int j = 0; j < settings.dim; j++)
            fprintf(fp, "%.15le ", vect[i * (settings.dim) + j]);
        fprintf(fp, "\n");
    }


    int sizeP, rankP;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

//    int proc_dim = settings.dim/ sizeP;
    if (settings.dim % sizeP != 0) {
        std::cout << "Wrong size" << std::endl;
        exit(-1);
    }
    int proc_row_size = settings.dim / sizeP;
    int proc_column_size = settings.dim / sizeP + 2;

    int proc_full_size = proc_row_size * proc_column_size;
    double *proc_vect = new double[proc_full_size];

    int stepCounter = 0;

    double procChange;
    double locChange;

    double tempPrevVal;
    double tempChange;


    // we should scatter our data

    int *displs = new int[sizeP];
    int *sendcounts = new int[sizeP];

    displs[0] = 0;
    sendcounts[0] = proc_full_size;

    for (int l = 1; l < sizeP; ++l) {
        displs[l] = proc_full_size * l - proc_row_size;
        sendcounts[l] = proc_full_size;
    }

    MPI_Scatterv(vect, sendcounts, displs, MPI_DOUBLE,
                 proc_vect, proc_full_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);


//    omp_set_num_threads(NUM_THREAD);
//    omp_lock_t globChangeLock;
//    omp_init_lock(&globChangeLock);
//
//    omp_lock_t *rowLock = new omp_lock_t[settings.dim];
//    for (int l = 0; l < settings.dim; ++l) {
//        omp_init_lock(&rowLock[l]);
//    }
//    time_S = omp_get_wtime();
//

    do {
//        std::cout << rankP << std::endl << std::flush;
//         swap data
        std::cout << "kek1" << std::endl;
        if (sizeP == 1) {
            std::cout << "kek2" << std::endl;
            ;
        } else if (rankP == ROOT) {
            std::cout << "kek3 " << rankP << std::endl;
            int offset = proc_row_size * proc_column_size - 1;
            MPI_Send(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP, 0, MPI_COMM_WORLD, &status);
//            MPI_Sendrecv(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP + 1, 0,
//                         proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP, 0, MPI_COMM_WORLD, &status);
            std::cout << "kek32 " << rankP << std::endl;
        } else if (rankP == sizeP - 1) {
            std::cout << "kek4 " << rankP << std::endl;
            MPI_Send(proc_vect, proc_row_size, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(proc_vect, proc_row_size, MPI_DOUBLE, rankP, 0, MPI_COMM_WORLD, &status);
//            MPI_Sendrecv(proc_vect, proc_row_size, MPI_DOUBLE, rankP - 1, 0,
//                         proc_vect, proc_row_size, MPI_DOUBLE, rankP, 0, MPI_COMM_WORLD, &status);
            std::cout << "kek42 " << rankP << std::endl;
        } else {
            std::cout << "kek5" << std::endl;
            // first row
            MPI_Sendrecv(proc_vect, proc_row_size, MPI_DOUBLE, rankP - 1, 0,
                         proc_vect, proc_row_size, MPI_DOUBLE, rankP, 0, MPI_COMM_WORLD, &status);
            // last row
            int offset = proc_row_size * proc_column_size - 1;
            MPI_Sendrecv(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP + 1, 0,
                         proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP, 0, MPI_COMM_WORLD, &status);
        }
        std::cout << "kek" << std::endl;

        procChange = 0;
////        #pragma omp parallel for shared(vect, settings, globChange) private(tempPrevVal, tempChange, locChange)
        for (int j = 1; j < proc_column_size - 1; ++j) {    // rows
            locChange = 0;
//
//            omp_set_lock(&rowLock[j+1]);
//            omp_set_lock(&rowLock[j]);
//            omp_set_lock(&rowLock[j-1]);
//
            for (int i = 1; i < proc_row_size - 1; ++i) {    // colms
//
                tempPrevVal = proc_vect[i * proc_row_size + j];
                proc_vect[i * proc_row_size + j] = 0.25 * (proc_vect[i * proc_row_size + j + 1] +
                                                          proc_vect[i * proc_row_size + j - 1] +
                                                          proc_vect[(i + 1) * proc_row_size + j] +
                                                          proc_vect[(i - 1) * proc_row_size + j]);
                tempChange = fabs(vect[proc_row_size * i + j] - tempPrevVal);
                if (locChange < tempChange) {
                    locChange = tempChange;
                }
            }
//
//            omp_set_lock(&globChangeLock);
            if (procChange < locChange) {
                procChange = locChange;
            }
//            omp_unset_lock(&globChangeLock);
//
//            omp_unset_lock(&rowLock[j-1]);
//            omp_unset_lock(&rowLock[j]);
//            omp_unset_lock(&rowLock[j+1]);
        }
        ++stepCounter;
        std::cout << "procChange: " << procChange << std::endl;
    } while ( procChange > settings.epsilon);
//
//
//    time_E = omp_get_wtime();
//
//    string filename2 = "../../result/result_openmp.txt";
//    FILE *fp2;
//    fp2 = fopen(filename2.c_str(), "w");
////
//    for (int i = 0; i < settings.dim; i++) {
//        for (int j = 0; j < settings.dim; j++)
//            fprintf(fp2, "%.15le ", vect[i * settings.dim + j]);
//        fprintf(fp2, "\n");
//    }
//
//    fclose(fp2);
//
//    printf("Proc count:\t %d\n", NUM_THREAD);
//    printf("Epsilon:\t %lf\n", settings.epsilon);
//    printf("Dim size:\t %d\n", settings.dim);
//    printf("Step calc:\t %d\n", stepCounter);
//    printf("Run time:\t %.15lf\n", time_E-time_S);

                MPI_Finalize();

                return 0;
            }