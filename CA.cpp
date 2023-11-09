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
int rowsPerProc, colsPerProc;

void testInit(){
    for(int i=0; i<rowsPerProc+2; i++){
        for(int j=0; j<colsPerProc+2; j++){
            readM[v(i,j)] = rank;
        }
    }
}

void printTest(){
    printf("RANK : %d\n", rank);
    for(int i=0; i<rowsPerProc+2; i++){
        for(int j=0; j<colsPerProc+2; j++){
            printf("%d ", readM[v(i,j)]);
        }
        printf("\n");
    }
    printf("\n");
}

void init(){
    for(int i=0; i<rowsPerProc+2; i++){
        for(int j=0; j<colsPerProc+2; j++){
            readM[v(i,j)] = 0;
        }
    }
}

void initCA(){
    //Random initialization of Game Of Life
    for(int i=1; i<rowsPerProc+1; i++){
        for(int j=1; j<colsPerProc+1; j++){
            readM[v(i,j)] = rand()%2;
        }
    }
}

void swap(){}

inline void send(){
    MPI_Request req;
    //Send the borders to the neighbors
    MPI_Send(&readM[v(1,1)], 1, rowType, rankUp, 17, cartComm);
    MPI_Send(&readM[v(rowsPerProc,1)], 1, rowType, rankDown, 18, cartComm);
    MPI_Send(&readM[v(1,1)], 1, columnType, rankLeft, 19, cartComm);
    MPI_Send(&readM[v(1,colsPerProc)], 1, columnType, rankRight, 20, cartComm);

    //Send the corners to the diagonal neighbors
    MPI_Send(&readM[v(1,1)], 1, cornerType, rankDiagonalUpLeft, 21, cartComm);
    MPI_Send(&readM[v(1,colsPerProc)], 1, cornerType, rankDiagonalUpRight, 22, cartComm);
    MPI_Send(&readM[v(rowsPerProc,1)], 1, cornerType, rankDiagonalDownLeft, 23, cartComm);
    MPI_Send(&readM[v(rowsPerProc,colsPerProc)], 1, cornerType, rankDiagonalDownRight, 24, cartComm);

}

inline void recv(){
    MPI_Status status;
    //Receive the borders from the neighbors
    MPI_Recv(&readM[v(rowsPerProc+1,0)],1,rowType,rankDown,17,cartComm, &status);
    MPI_Recv(&readM[v(0,0)],1,rowType,rankUp,18,cartComm,&status);
    MPI_Recv(&readM[v(0,colsPerProc+1)],1,columnType,rankRight,19,cartComm,&status);
    MPI_Recv(&readM[v(0,0)],1,columnType,rankLeft,20,cartComm,&status);

    //Receive the diagonals from the neighbors
    MPI_Recv(&readM[v(rowsPerProc+2,colsPerProc+2)],1,cornerType,rankDiagonalDownRight,21,cartComm,&status);
    MPI_Recv(&readM[v(rowsPerProc+2,0)],1,cornerType,rankDiagonalDownLeft,22,cartComm,&status);
    MPI_Recv(&readM[v(rowsPerProc+2,colsPerProc+2)],1,cornerType,rankDiagonalUpRight,23,cartComm,&status);
    MPI_Recv(&readM[v(0,0)],1,cornerType,rankDiagonalUpLeft,24,cartComm,&status);
}

inline void transFuncInside(){}

inline void transFuncBorders(){}

inline void transFunc(){}

void printCoords(MPI_Comm cartComm){
    int myCoords[2] = {0, 0};
    MPI_Cart_coords(cartComm, rank, 2, myCoords);
    // Print all the ranks and their coordinates
    printf("Rank: %d, coords: (%d,%d)\n rankLeft: %d, rankRight: %d, rankUp: %d, rankDown: %d\n rankDiagonalUpLeft: %d, rankDiagonalUpRight: %d, rankDiagonalDownLeft: %d, rankDiagonalDownRight: %d\n\n", rank, myCoords[0], myCoords[1], rankLeft, rankRight, rankUp, rankDown, rankDiagonalUpLeft, rankDiagonalUpRight, rankDiagonalDownLeft, rankDiagonalDownRight);
}

//MAIN
int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int dims[2] = {0, 0}; 
    MPI_Dims_create(nProc, 2, dims); 
    int periods[2] = {1,1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cartComm);
    MPI_Cart_shift(cartComm, 0, 1, &rankUp, &rankDown);
    MPI_Cart_shift(cartComm, 1, 1, &rankLeft, &rankRight);

    //Calculate the ranks of the processes in the diagonal
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

    //printCoords(cartComm); //Todo: Remove

    rowsPerProc = NROWS/dims[0];
    colsPerProc = NCOLS/dims[1];
    printf("rowsPerProc: %d, colsPerProc: %d\n", rowsPerProc, colsPerProc);
    readM = new int[(rowsPerProc+2)*(colsPerProc+2)];
    writeM = new int[(rowsPerProc+2)*(colsPerProc+2)];

    MPI_Type_vector(NROWS+2, 1, NCOLS+2, MPI_INT, &columnType); //NROWS+4 or NROWS?
    MPI_Type_commit(&columnType);
    MPI_Type_vector(NCOLS+2, 1, 1, MPI_INT, &rowType);
    MPI_Type_commit(&rowType);
    MPI_Type_vector(1, 1, 1, MPI_INT, &cornerType); //Useless(?) but created for consistency
    MPI_Type_commit(&cornerType);

    //init();
    testInit();

    //initCA();

    for(int i = 0; i < NGENS; i++){
        send();
        //transFuncInside();
        recv();
        printTest();
        //transFuncBorders();
        //Allegro here(?)
        //swap();
    }

    MPI_Type_free(&columnType);
    MPI_Type_free(&rowType);
    MPI_Type_free(&cornerType);
    delete[] readM;
    delete[] writeM;

    MPI_Finalize();
    return 0;
}

