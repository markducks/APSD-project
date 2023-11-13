#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "CA.h"

#define NCOLS 9 //Dummy values
#define NROWS 9 //Dummy values
#define NGENS 1 //Dummy values
 
#define v(r,c) ((r)*(colsPerProc+2) + (c))

int nProc, rank, rankLeft, rankRight, rankUp, rankDown; 
int rankDiagonalUpLeft, rankDiagonalUpRight, rankDiagonalDownLeft, rankDiagonalDownRight;
MPI_Datatype columnType, rowType, cornerType; 
MPI_Comm cartComm;
int* readM;
int* writeM;
int dims[2];
int rowsPerProc, colsPerProc;

void testInit(){
    for(int i=0; i<rowsPerProc+2; i++){
        for(int j=0; j<colsPerProc+2; j++){
            readM[v(i,j)] = rank;
        }
    }
}


void init(){
    for(int i=0; i<rowsPerProc+2; i++){
        for(int j=0; j<colsPerProc+2; j++){
            readM[v(i,j)] = 0;
        }
    }
}

void initRandomCA(){
    //Random initialization of Game Of Life
    for(int i=1; i<rowsPerProc+1; i++){
        for(int j=1; j<colsPerProc+1; j++){
            readM[v(i,j)] = rand()%2;
        }
    }
}


inline void send(){
    MPI_Request req;

    MPI_Isend(&readM[v(1,0)], 1, rowType, rankUp, 17, cartComm, &req);
    MPI_Isend(&readM[v(rowsPerProc,0)], 1, rowType, rankDown, 18, cartComm,  &req);
    MPI_Isend(&readM[v(0,1)], 1, columnType, rankLeft, 19, cartComm, &req);
    MPI_Isend(&readM[v(0,colsPerProc)], 1, columnType, rankRight, 20, cartComm, &req);

    
    //Send the corners to the diagonal neighbors
    MPI_Isend(&readM[v(1,1)], 1, cornerType, rankDiagonalUpLeft, 21, cartComm, &req);
    MPI_Isend(&readM[v(1,colsPerProc)], 1, cornerType, rankDiagonalUpRight, 22, cartComm, &req);
    MPI_Isend(&readM[v(rowsPerProc,1)], 1, cornerType, rankDiagonalDownLeft, 23, cartComm, &req);
    MPI_Isend(&readM[v(rowsPerProc,colsPerProc)], 1, cornerType, rankDiagonalDownRight, 24, cartComm, &req);
    
}

inline void recv(){
    MPI_Status status;
    //Receive the borders from the neighbors
    MPI_Recv(&readM[v(rowsPerProc+1,0)],1,rowType,rankDown,17,cartComm, &status);
    MPI_Recv(&readM[v(0,0)],1,rowType,rankUp,18,cartComm,&status);
    MPI_Recv(&readM[v(0,colsPerProc+1)],1,columnType,rankRight,19,cartComm,&status);
    MPI_Recv(&readM[v(0,0)],1,columnType,rankLeft,20,cartComm,&status);
    
    
    //Receive the diagonals from the neighbors
    MPI_Recv(&readM[v(rowsPerProc+1,colsPerProc+1)],1,cornerType,rankDiagonalDownRight,21,cartComm,&status);
    MPI_Recv(&readM[v(rowsPerProc+1,0)],1,cornerType,rankDiagonalDownLeft,22,cartComm,&status);
    MPI_Recv(&readM[v(0,colsPerProc+1)],1,cornerType,rankDiagonalUpRight,23,cartComm,&status);
    MPI_Recv(&readM[v(0,0)],1,cornerType,rankDiagonalUpLeft,24,cartComm,&status);
    
}

inline void transFunc(int i, int j){
    int nNeighbors = 0;
    nNeighbors += readM[v(i-1,j-1)];
    nNeighbors += readM[v(i-1,j)];
    nNeighbors += readM[v(i-1,j+1)];
    nNeighbors += readM[v(i,j-1)];
    nNeighbors += readM[v(i,j+1)];
    nNeighbors += readM[v(i+1,j-1)];
    nNeighbors += readM[v(i+1,j)];
    nNeighbors += readM[v(i+1,j+1)];
    if(readM[v(i,j)] == 1){
        if(nNeighbors < 2 || nNeighbors > 3){
            writeM[v(i,j)] = 0;
        }
        else{
            writeM[v(i,j)] = 1;
        }
    }
    else{
        if(nNeighbors == 3){
            writeM[v(i,j)] = 1;
        }
        else{
            writeM[v(i,j)] = 0;
        }
    }
}

inline void transFuncInside(){
    for(int i=1; i<rowsPerProc; i++){  //is correct?
        for(int j=1; j<colsPerProc; j++){  //is correct?
            transFunc(i,j);
        }
    }
}

inline void transFuncBorders(){
    for(int i=0; i<rowsPerProc; i++){  //is correct?
        transFunc(i,0);
        transFunc(i,colsPerProc+1);
    }
    for(int j=0; j<colsPerProc; j++){  //is correct?
        transFunc(0,j);
        transFunc(rowsPerProc+1,j);
    }
}

void swap(){
    int* temp = readM;
    readM = writeM;
    writeM = temp;
}

void printCoords(MPI_Comm cartComm){
    int myCoords[2] = {0, 0};
    MPI_Cart_coords(cartComm, rank, 2, myCoords);
    // Print all the ranks and their coordinates
    printf("Rank: %d, coords: (%d,%d)\n rankLeft: %d, rankRight: %d, rankUp: %d, rankDown: %d\n rankDiagonalUpLeft: %d, rankDiagonalUpRight: %d, rankDiagonalDownLeft: %d, rankDiagonalDownRight: %d\n\n", rank, myCoords[0], myCoords[1], rankLeft, rankRight, rankUp, rankDown, rankDiagonalUpLeft, rankDiagonalUpRight, rankDiagonalDownLeft, rankDiagonalDownRight);
}

inline void printGather(){
    int* gatherBuffer = NULL;
    if(rank == 0){
        gatherBuffer = new int[nProc * (rowsPerProc+2) * (colsPerProc+2)];
    }

    MPI_Gather(readM, (rowsPerProc+2) * (colsPerProc+2), MPI_INT, gatherBuffer, (rowsPerProc+2) * (colsPerProc+2), MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < nProc; i++) {
            printf("Matrix from rank %d:\n", i);
            for (int j = 0; j < (rowsPerProc+2); j++) {
                for (int k = 0; k < (colsPerProc+2); k++) {
                    printf("%d ", gatherBuffer[i * (rowsPerProc+2) * (colsPerProc+2) + j * (colsPerProc+2) + k]);
                }
                printf("\n");
            }
            printf("\n");
        }
        delete[] gatherBuffer;
    }
}

inline void rankDiscovering(){
    dims[0] = dims[1] = 0;
    MPI_Dims_create(nProc, 2, dims); 
    int periods[2] = {1,1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cartComm);
    MPI_Cart_shift(cartComm, 0, 1, &rankUp, &rankDown);
    MPI_Cart_shift(cartComm, 1, 1, &rankLeft, &rankRight);

    int coords[2];
    MPI_Cart_coords(cartComm, rank, 2, coords); 
    coords[0] = (coords[0] - 1 + dims[0]) % dims[0];
    coords[1] = (coords[1] - 1 + dims[1]) % dims[1];
    MPI_Cart_rank(cartComm, coords, &rankDiagonalUpLeft);
    coords[1] = (coords[1] + 2) % dims[1];
    MPI_Cart_rank(cartComm, coords, &rankDiagonalUpRight);
    coords[0] = (coords[0] + 2) % dims[0];
    MPI_Cart_rank(cartComm, coords, &rankDiagonalDownRight);
    coords[1] = (coords[1] - 2 + dims[1]) % dims[1];
    MPI_Cart_rank(cartComm, coords, &rankDiagonalDownLeft);
}

//MAIN
int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
    rankDiscovering();

    rowsPerProc = NROWS/dims[0];
    colsPerProc = NCOLS/dims[1];
    readM = new int[(rowsPerProc+2)*(colsPerProc+2)];
    writeM = new int[(rowsPerProc+2)*(colsPerProc+2)];

    MPI_Type_vector(rowsPerProc+2, 1, colsPerProc+2, MPI_INT, &columnType);
    MPI_Type_commit(&columnType);
    MPI_Type_vector(colsPerProc+2, 1, 1, MPI_INT, &rowType);
    MPI_Type_commit(&rowType);
    MPI_Type_vector(1, 1, 1, MPI_INT, &cornerType);
    MPI_Type_commit(&cornerType);

    init();

    initRandomCA();

    for(int i = 0; i < NGENS; i++){
        send();
        transFuncInside();
        recv();
        printGather();
        transFuncBorders();
        //Allegro here(?)
        swap();
    }

    MPI_Type_free(&columnType);
    MPI_Type_free(&rowType);
    MPI_Type_free(&cornerType);
    delete[] readM;
    delete[] writeM;

    MPI_Finalize();
    return 0;
}

