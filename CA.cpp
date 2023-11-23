#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <allegro.h> //  /usr/include/allegro.h   /usr/lib/x86_64-linux-gnu/liballeg.so
#include "CA.h"

#define NCOLS 600 
#define NROWS 600
#define NGENS 1500
#define CELL_SIZE 3
 
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
int totalSubmatrixSize;

int black, white, red;
int dir = 1;
BITMAP *buffer;
int* gatherBuffer;

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


void init(){
    for(int i=0; i<rowsPerProc+2; i++){
        for(int j=0; j<colsPerProc+2; j++){
            readM[v(i,j)] = 0;
        }
    }
}



void initRandomCA(){

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



void allegroInit(){
   
	allegro_init();
    set_color_depth(8);
    set_gfx_mode( GFX_AUTODETECT_WINDOWED, NCOLS, NROWS, 0, 0);
	buffer = create_bitmap(NCOLS, NROWS);
	
    white = makecol(255, 255, 255);
	black = makecol(0, 0, 0);

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



void allegroDraw(int step){
   	int coords[2];
    MPI_Cart_coords(cartComm, rank, 2, coords);

    int globalRow = coords[0] * (rowsPerProc);
    int globalCol = coords[1] * (colsPerProc);



    for(int i = globalRow; i < globalRow+rowsPerProc; i++){
        for(int j = globalCol; j < globalCol+colsPerProc; j++){
			switch (readM[v(i,j)]) {
			case 0:
                rectfill(buffer, i*CELL_SIZE, j*CELL_SIZE, (i+1)*CELL_SIZE, (j+1)*CELL_SIZE, black);
                break;
            case 1:
                rectfill(buffer, i*CELL_SIZE, j*CELL_SIZE, (i+1)*CELL_SIZE, (j+1)*CELL_SIZE, white);
                break;
			}
        }
    }

	textprintf_ex(buffer, font, 0, 0, white, black, "Step: %d", step);

	blit(buffer, screen, 0, 0, 0, 0, rowsPerProc, colsPerProc);
}



inline void printWithCells(){
    

    for (int rank = 0; rank < nProc; rank++) {
   
        int coords[2];
        MPI_Cart_coords(cartComm, rank, 2, coords); 

        int globalRow = coords[0] * (rowsPerProc);
        int globalCol = coords[1] * (colsPerProc);


        for(int i = globalRow; i < globalRow+rowsPerProc; i++){
            for(int j = globalCol; j < globalCol+colsPerProc; j++){
                switch(gatherBuffer[g(i,j)]){
                case 0:
                    rectfill(buffer, i*CELL_SIZE, j*CELL_SIZE, (i+1)*CELL_SIZE, (j+1)*CELL_SIZE, black);
                    break;
                case 1:
                    rectfill(buffer, i*CELL_SIZE, j*CELL_SIZE, (i+1)*CELL_SIZE, (j+1)*CELL_SIZE, white);
                    break;
			}
            }
        }

    }

    blit(buffer, screen, 0, 0, 0, 0, NCOLS, NROWS );

}






void swap(){
    int* temp = readM;
    readM = writeM;
    writeM = temp;
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


void printColorsBasedOnRank(){
    int procOnRows = dims[0];
    int procOnCols = dims[1];

    for (int rank = 0; rank < nProc; rank++) {
   
        int coords[2];
        MPI_Cart_coords(cartComm, rank, 2, coords);

        int globalRow = coords[0] * (rowsPerProc+2);
        int globalCol = coords[1] * (colsPerProc+2);


        for(int i = globalRow; i < globalRow+rowsPerProc+2; i++){
            for(int j = globalCol; j < globalCol+colsPerProc+2; j++){
               switch(rank){
                    case 0:
                        putpixel(buffer, i, j, makecol(255, 0, 0));
                        break;
                    case 1:
                        putpixel(buffer, i, j, makecol(0, 255, 0));
                        break;
                    case 2:
                        putpixel(buffer, i, j, makecol(0, 0, 255));
                        break;
                    case 3:
                        putpixel(buffer, i, j, makecol(255, 255, 0));
                        break;
                    case 4:
                        putpixel(buffer, i, j, makecol(255, 0, 255));
                        break;
                    case 5:
                        putpixel(buffer, i, j, makecol(0, 255, 255));
                        break;
                    case 6:
                        putpixel(buffer, i, j, makecol(255, 255, 255));
                        break;
                    case 7:
                        putpixel(buffer, i, j, makecol(0, 0, 0));
                        break;
                    case 8:
                        putpixel(buffer, i, j, makecol(255, 128, 0));
                        break;
                    case 9:
                        putpixel(buffer, i, j, makecol(255, 0, 128));
                        break;
                    case 10:
                        putpixel(buffer, i, j, makecol(128, 255, 0));
                        break;
               }
			}
        }

        blit(buffer, screen, 0, 0, 0, 0, NCOLS, NROWS );
    }
}

inline void gatherAndPrint(int i){
    MPI_Gather(&readM[v(1,1)], 1, gatherType, gatherBuffer, totalSubmatrixSize, MPI_INT, 0, cartComm);
        if(rank == 0){
            printWithCells();
            //printColorsBasedOnRank();
        }
}


int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
    rankDiscovering();

    rowsPerProc = NROWS/dims[0];
    colsPerProc = NCOLS/dims[1];
    totalSubmatrixSize = (rowsPerProc)*(colsPerProc);

    readM = new int[(rowsPerProc+2)*(colsPerProc+2)];
    writeM = new int[(rowsPerProc+2)*(colsPerProc+2)];

    MPI_Type_vector(rowsPerProc+2, 1, colsPerProc+2, MPI_INT, &columnType);
    MPI_Type_commit(&columnType);
    MPI_Type_vector(colsPerProc+2, 1, 1, MPI_INT, &rowType);
    MPI_Type_commit(&rowType);
    MPI_Type_vector(1, 1, 1, MPI_INT, &cornerType);
    MPI_Type_commit(&cornerType);

    init(); //Inizializza la matrice a 0

    initRandomCA(); //Inizializza la matrice con valori random
    //initLWSS();

    //initGlider();

    MPI_Type_vector(rowsPerProc, colsPerProc, colsPerProc+2, MPI_INT, &gatherType);
    MPI_Type_commit(&gatherType);

    int procOnRows = dims[0];
    int procOnCols = dims[1];

    gatherBuffer = new int[NROWS*NCOLS];
    
    if(rank == 0){
        allegroInit();
    }

    /*
    double startTime = MPI_Wtime();
    */

    for(int i = 0; i < NGENS; i++){
        send();
        transFuncInside();
        recv();
        transFuncBorders();
        gatherAndPrint(i);
        swap();
        printf("Step %d\n",i);
    }

    /*
    double endTime = MPI_Wtime();
    double elapsedTime = endTime - startTime;
    if (rank == 0) {  
        printf("Time elapsed: %f seconds\n", elapsedTime);
    }
    */

    MPI_Type_free(&columnType);
    MPI_Type_free(&rowType);
    MPI_Type_free(&cornerType);
    delete[] readM;
    delete[] writeM;
    if(rank == 0){
        delete[] gatherBuffer;
    }

    
    MPI_Finalize();
    return 0;
}

