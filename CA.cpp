#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NCOLS 2000 //Dummy values
#define NROWS 2000 //Dummy values
#define NGENS 2000 //Dummy values

#define NCELLS ((NCOLS*NROWS)/nProc)+4 // +4 for ghost cells???

int* readM;
int* writeM;

int nProc, rank, rankLeft, rankRight, rankUp, rankDown; 
MPI_Datatype columnType, rowType, cornerType; //Is cornerType needed?

//Todo: Create all the side functions

void init(){}

void initCA(){}

void swap(){}

inline void sendBorders(){}

inline void recvBorders(){}

inline void transFuncInside(){}

inline void transFuncBorders(){}

inline void transFunc(){}


int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    readM = new int[NCELLS];
    writeM = new int[NCELLS];

    MPI_Comm cartComm;
    int dims[2] = {0, 0}; 
    MPI_Dims_create(nProc, 2, dims); 
    int periods[2] = {1,1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cartComm);
    MPI_Cart_shift(cartComm, 0, 1, &rankLeft, &rankRight);
    MPI_Cart_shift(cartComm, 1, 1, &rankUp, &rankDown);

    //Todo: Create the columnType, rowType and cornerType

    //Todo: Initialize the matrix
    init();

    //Todo: Create the initial configuration
    initCA();

    //Todo: Create the main loop
    for(int i = 0; i < NGENS; i++){
        sendBorders();
        transFuncInside();
        recvBorders();
        transFuncBorders();
        //Allegro here(?)
        swap();
    }

    MPI_Finalize();
    return 0;
}

