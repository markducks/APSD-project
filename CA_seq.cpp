#include <stdio.h>
#include <stdlib.h>
#include <chrono>
using namespace std::chrono;

#define NCOLS 120 
#define NROWS 180
#define NGENS 3500  

#define v(r,c) ((r)*(NCOLS) + (c))

int* readM;
int* writeM;

void init(){
    for(int i=0; i<NROWS; i++){
        for(int j=0; j<NCOLS; j++){
            readM[v(i,j)] = 0;
        }
    }
}

void initRandomCA(){
    //Random initialization of Game Of Life
    for(int i=0; i<NROWS; i++){
        for(int j=0; j<NCOLS; j++){
            readM[v(i,j)] = rand()%2;
        }
    }
}

void initLWSS(){
    for(int i=0; i<NROWS; i++){
        for(int j=0; j<NCOLS; j++){
            readM[v(i,j)] = 0;
        }
    }

    readM[v(1,2)] = 1;
    readM[v(1,3)] = 1;
    readM[v(1,4)] = 1;
    readM[v(1,5)] = 1;
    readM[v(2,1)] = 1;
    readM[v(2,5)] = 1;
    readM[v(3,5)] = 1;
    readM[v(4,4)] = 1;
}

inline void transFuncCell(int i, int j){
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

inline void transFunc(){
    for(int i=1; i<NROWS; i++){
        for(int j=1; j<NCOLS; j++){
            transFuncCell(i,j);
        }
    }
}

void swap(){
    int* temp = readM;
    readM = writeM;
    writeM = temp;
}


int main(int argc, char *argv[]){
    
    readM = new int[NCOLS*NROWS];
    writeM = new int[NCOLS*NROWS];

    init();

    initLWSS();

    auto start = high_resolution_clock::now();

    for(int i = 0; i < NGENS; i++){
        transFunc();    
        swap();
    }

    auto stop = high_resolution_clock::now();
    auto elapsedTime = duration_cast<microseconds>(stop - start).count() / 1000000.0;
    printf("Time elapsed: %f seconds\n", elapsedTime);

    return 0;
}

