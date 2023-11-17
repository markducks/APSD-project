#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <allegro.h> //  /usr/include/allegro.h   /usr/lib/x86_64-linux-gnu/liballeg.so
#include "CA.h"

#define NCOLS 120 
#define NROWS 180
#define NGENS 3500 
#define CELL_SIZE 2 
 
#define v(r,c) ((r)*(colsPerProc+2) + (c))

#define g(i, j) ((i)*NCOLS + (j))


int nProc, rank, rankLeft, rankRight, rankUp, rankDown; 
int rankDiagonalUpLeft, rankDiagonalUpRight, rankDiagonalDownLeft, rankDiagonalDownRight;
MPI_Datatype columnType, rowType, cornerType, gatherType; 
MPI_Comm cartComm;
int* readM;
int* writeM;
int dims[2];
int rowsPerProc, colsPerProc;

int black, white;
int dir = 1;
BITMAP *buffer;
int* gatherBuffer;

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

inline void initLWSS(){
    readM[v(1,2)] = 1;
    readM[v(1,3)] = 1;
    readM[v(1,4)] = 1;
    readM[v(1,5)] = 1;
    readM[v(2,1)] = 1;
    readM[v(2,5)] = 1;
    readM[v(3,5)] = 1;
    readM[v(4,4)] = 1;
}

void initGlider(){
    // Clear the grid
    for(int i=0; i<rowsPerProc+2; i++){
        for(int j=0; j<colsPerProc+2; j++){
            readM[v(i,j)] = 0;
        }
    }

    // Set the glider pattern
    readM[v(1,2)] = 1;
    readM[v(2,3)] = 1;
    readM[v(3,1)] = 1;
    readM[v(3,2)] = 1;
    readM[v(3,3)] = 1;
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
    for(int i=2; i<rowsPerProc; i++){  //is correct?
        for(int j=2; j<colsPerProc; j++){  //is correct?
            transFunc(i,j);
        }
    }
}

inline void transFuncBorders(){
    for(int i=1; i<rowsPerProc+1; i++){
        transFunc(i,1);
        transFunc(i,colsPerProc);
    }
    for(int j=1; j<colsPerProc+1;j++){
        transFunc(1,j);
        transFunc(rowsPerProc,j);
    }
}

void swap(){
    int* temp = readM;
    readM = writeM;
    writeM = temp;
}

void allegroInit(){
   
	allegro_init();
	set_color_depth(8);
	buffer = create_bitmap(NCOLS, NROWS);
	set_gfx_mode( GFX_AUTODETECT_WINDOWED, NCOLS, NROWS, 0, 0);

    white = makecol(255, 255, 255);
	

}

void initGather(){
    MPI_Type_vector(rowsPerProc, colsPerProc, colsPerProc+2, MPI_INT, &gatherType);
    MPI_Type_commit(&gatherType);
    gatherBuffer = new int[NROWS*NCOLS];
}

void gatherData(){
    int myCoords[2] = {0, 0};
    MPI_Cart_coords(cartComm, rank, 2, myCoords);
    // Print all the ranks and their coordinates
    printf("Rank: %d, coords: (%d,%d)\n rankLeft: %d, rankRight: %d, rankUp: %d, rankDown: %d\n rankDiagonalUpLeft: %d, rankDiagonalUpRight: %d, rankDiagonalDownLeft: %d, rankDiagonalDownRight: %d\n\n", rank, myCoords[0], myCoords[1], rankLeft, rankRight, rankUp, rankDown, rankDiagonalUpLeft, rankDiagonalUpRight, rankDiagonalDownLeft, rankDiagonalDownRight);
}

inline void print(){

    for (int i = 0; i < NROWS; i++) {
        for (int j = 0; j < NCOLS; j++) {
            switch (gatherBuffer[g(i,j)]) {
                case 0:
                    putpixel(buffer, i, j, black);
                    break;
                case 1:
                    putpixel(buffer, i, j, white);
                    break;
        };
    };

    blit(buffer, screen, 0, 0, 0, 0, NROWS, NCOLS);
    }
}

inline void printWithCells(){
    for (int i = 0; i < NROWS; i++) {
        for (int j = 0; j < NCOLS; j++) {
            if(gatherBuffer[g(i,j)] == 1){
                rectfill(buffer, j * CELL_SIZE, i * CELL_SIZE, j*CELL_SIZE+CELL_SIZE, i*CELL_SIZE+CELL_SIZE, white);
            }
            else{
                rectfill(buffer, j * CELL_SIZE, i * CELL_SIZE, j*CELL_SIZE+CELL_SIZE, i*CELL_SIZE+CELL_SIZE, black);
            }
        };
    };
    blit(buffer, screen, 0, 0, 0, 0, NROWS, NCOLS);
}

inline void gatherAndPrint(){
    MPI_Gather(&readM[v(1,1)], 1, gatherType, gatherBuffer, 1, gatherType, 0, MPI_COMM_WORLD);
        if(rank == 0){
            printWithCells();
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

void allegroDraw(int step){
   	//draw every pixel with a color depending on the state
	for (int row = 0; row < rowsPerProc; row++)
		for (int col = 0; col < colsPerProc; col++)
			switch (readM[v(row,col)]) {
			case 0:
				putpixel(buffer, col, row, black);
				break;
			case 1:
				putpixel(buffer, col, row, white);
				break;
			}

	textprintf_ex(buffer, font, 0, 0, white, black, "Step: %d", step);

	blit(buffer, screen, 0, 0, 0, 0, rowsPerProc, colsPerProc);
}

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

    //initRandomCA();
    initLWSS();

    /*
    initGather();

    if(rank == 0){
        allegroInit();
    }*/

    double startTime = MPI_Wtime();

    for(int i = 0; i < NGENS; i++){
        send();
        transFuncInside();
        recv();
        transFuncBorders();
        swap();
    }

    double endTime = MPI_Wtime();
    double elapsedTime = endTime - startTime;
    if (rank == 0) {  // Only print the time on the root process
        printf("Time elapsed: %f seconds\n", elapsedTime);
    }


    MPI_Type_free(&columnType);
    MPI_Type_free(&rowType);
    MPI_Type_free(&cornerType);
    delete[] readM;
    delete[] writeM;

    MPI_Finalize();
    return 0;
}

