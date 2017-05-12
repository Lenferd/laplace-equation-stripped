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

    int sizeP, rankP;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    Settings settings;
    double *vect = new double[settings.vectSize];
    if (rankP == ROOT) {
        FILE *infile = fopen("../../initial/settings.ini", "r");

        if (infile == NULL) {
            std::cout << "File open error" << std::endl;
            exit(-1);
        }

        fscanf(infile, "DIM=%d\n", &settings.dim);
        fscanf(infile, "EPS=%lf\n", &settings.epsilon);
        fscanf(infile, "XSTART=%lf\n", &settings.xStart);
        fscanf(infile, "XEND=%lf\n", &settings.xEnd);
        fscanf(infile, "YSTART=%lf\n", &settings.yStart);
        fscanf(infile, "YEND=%lf\n", &settings.yEnd);
        settings.vectSize = settings.dim * settings.dim;    // with +2 boundaries


        double time_S, time_E;

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
            for (int j = 0; j < settings.dim; j++) {
                fprintf(fp, "%.15le ", vect[i * (settings.dim) + j]);
            }
            fprintf(fp, "\n");
        }
    }

    MPI_Bcast(&settings.vectSize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&settings.dim, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&settings.epsilon, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&settings.xStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&settings.xEnd, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&settings.yStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&settings.yEnd, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);


//    int proc_dim = settings.dim/ sizeP;
    if (settings.dim % sizeP != 0) {
        std::cout << "Wrong size" << std::endl;
        exit(-1);
    }
    int proc_row_size = settings.dim;
    int proc_column_size = settings.dim / sizeP + 2;

    int proc_full_size = proc_row_size * proc_column_size;
    double *proc_vect = new double[proc_full_size];

    int stepCounter = 0;

    double globChange;
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
        displs[l] = proc_full_size * l - proc_row_size * 2;
        sendcounts[l] = proc_full_size;
    }

    MPI_Scatterv(vect, sendcounts, displs, MPI_DOUBLE,
                 proc_vect, proc_full_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    std::cout << "proc_full_size: " << proc_full_size << std::endl;
    std::cout << "settings.VectSize: " << settings.vectSize << std::endl;
    std::cout << "row size: " << proc_row_size << std::endl;
    std::cout << "column size: " << proc_column_size << std::endl;

    if (sizeP == 1) {
        exit(0);
    }
    do {
//        int dest = (rankP + 1) % sizeP;
//        int source = (rankP - 1 + sizeP) % sizeP;

        if (rankP == ROOT) {
            int offset = proc_full_size - proc_row_size * 2;    // send penult row
            MPI_Send(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD);
            offset = proc_full_size - proc_row_size;        // get las row
            MPI_Recv(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, &status);
        } else if (rankP == sizeP - 1) {
            int offset = proc_row_size;
            MPI_Send(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(proc_vect, proc_row_size, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, &status);
        } else {
            int offset = proc_full_size - proc_row_size * 2;    // send penult row
            MPI_Send(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD);
            offset = proc_row_size;
            MPI_Send(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD);

            offset = proc_full_size - proc_row_size;        // get las row
            MPI_Recv(proc_vect + offset, proc_row_size, MPI_DOUBLE, rankP + 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(proc_vect, proc_row_size, MPI_DOUBLE, rankP - 1, 0, MPI_COMM_WORLD, &status);
        }

//        std::cout << "kek e" << std::endl;

        procChange = 0;
        for (int j = 1; j < proc_column_size - 1; ++j) {    // rows
            locChange = 0;

            for (int i = 1; i < proc_row_size - 1; ++i) {    // colms
//
                tempPrevVal = proc_vect[j * proc_row_size + i];
                proc_vect[j * proc_row_size + i] = 0.25 * (proc_vect[j * proc_row_size + i + 1] +
                                                           proc_vect[j * proc_row_size + i - 1] +
                                                           proc_vect[(j + 1) * proc_row_size + i] +
                                                           proc_vect[(j - 1) * proc_row_size + i]);
//                if (j == 2 and i == 5) {
//                    std::cout << "proc_vect[j * proc_row_size + i]: " << rankP << "   "
//                              << proc_vect[j * proc_row_size + i] << std::endl;
//                }
                tempChange = fabs(proc_vect[proc_row_size * j + i] - tempPrevVal);
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
        MPI_Reduce(&procChange, &globChange, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);

        if (rankP == ROOT) {
            std::cout << "globChange: " << globChange << std::endl;
//        std::cout << "===" << std::endl;
        }
        MPI_Bcast(&globChange, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    } while (globChange > settings.epsilon);
//
    displs = new int[sizeP];
    int *recvcounts = new int[sizeP];

    if (rankP == 1) {
        for (int i = 0; i < proc_column_size; i++) {
            for (int j = 0; j < proc_row_size; j++)
                std::cout << proc_vect[i * proc_row_size + j] << " ";
            std::cout << std::endl;
        }
    }

    displs[0] = 0;
    recvcounts[0] = proc_full_size;

    for (int l = 1; l < sizeP; ++l) {
        displs[l] = proc_full_size * l - proc_row_size * 2;
        recvcounts[l] = proc_full_size;
    }

    double *recv_vect = new double[proc_full_size * sizeP];
//    MPI_Gather(proc_vect, proc_full_size, MPI_DOUBLE, vect, proc_full_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    MPI_Gatherv(proc_vect, proc_full_size, MPI_DOUBLE, recv_vect, recvcounts, displs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
//    MPI_Gatherv(proc_vect, sendcounts, displs, MPI_DOUBLE, vect, proc_full_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD)
//
//    time_E = omp_get_wtime();
    std::cout << "EXIT rankP: " << rankP << std::endl;
    if (rankP == ROOT) {

        string filename2 = "../../result/result_" + std::to_string(sizeP)+ "_mpi.txt";
        FILE *fp2;
        fp2 = fopen(filename2.c_str(), "w");


        for (int i = 0; i < settings.dim; i++) {
            for (int j = 0; j < settings.dim; j++)
                fprintf(fp2, "%.15le ", recv_vect[i * settings.dim + j]);
            fprintf(fp2, "\n");
        }

        fclose(fp2);

        printf("Proc count:\t %d\n", NUM_THREAD);
        printf("Epsilon:\t %lf\n", settings.epsilon);
        printf("Dim size:\t %d\n", settings.dim);
        printf("Step calc:\t %d\n", stepCounter);
//    printf("Run time:\t %.15lf\n", time_E-time_S);

    }
    MPI_Finalize();

    return 0;
}