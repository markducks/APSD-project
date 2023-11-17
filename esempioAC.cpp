#include <allegro.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
using namespace std;

// File di prova per Automi Cellulari
// Usa dati di setup dal file Configuration.txt

int **readMatrix;
int **writeMatrix;

int numStep;
int numRows;
int numCols;
int rowCentro;
int colCentro;
int raggio;

int black, white;
int dir = 1;
BITMAP *buffer;

void initAllegro() {


	allegro_init();
	install_keyboard();
	set_color_depth(24);
	buffer = create_bitmap(numCols, numRows);
	set_gfx_mode( GFX_AUTODETECT_WINDOWED, numCols, numRows, 0, 0);

	black = makecol(0, 0, 0);
	white = makecol(255, 255, 255);


}

void readConfigFile() {
	char str[20];
	FILE *file;
	file = fopen("Configuration.txt", "r");
	if (file) {
		fscanf(file, "%s", str);
		numRows = atoi(str);
		fscanf(file, "%s", str);
		numCols = atoi(str);
		fscanf(file, "%s", str);
		rowCentro = atoi(str);
		fscanf(file, "%s", str);
		colCentro = atoi(str);
		fscanf(file, "%s", str);
		raggio = atoi(str);
		fscanf(file, "%s", str);
		numStep = atoi(str);
		fclose(file);
	}

}

void createMatrix(int **&m, int numRows, int numCols) {
	m = new int*[numRows];
	for (int row = 0; row < numRows; ++row) {
		m[row] = new int[numCols];
	}
}

void initModel() {
	createMatrix(readMatrix, numRows, numCols);
	createMatrix(writeMatrix, numRows, numCols);
	for (int row = 0; row < numRows; row++) {
		for (int col = 0; col < numCols; col++) {
			if (sqrt(pow((row - rowCentro), 2) + pow((col - colCentro), 2))
					<= raggio) {
				readMatrix[row][col] = 1;
				writeMatrix[row][col] = 1;
			} else {
				readMatrix[row][col] = 0;
				writeMatrix[row][col] = 0;
			}
		}

	}
}

void drawWithAllegro(int step) {
	//draw every pixel with a color depending on the state
	for (int row = 0; row < numRows; row++)
		for (int col = 0; col < numCols; col++)
			switch (readMatrix[row][col]) {
			case 0:
				putpixel(buffer, col, row, black);
				break;
			case 1:
				putpixel(buffer, col, row, white);
				break;
			}

	textprintf_ex(buffer, font, 0, 0, white, black, "step: %i", step);
	//Move the modified buffer to the screen
//	printf("numRows=%d numCols=%d\n",numRows,numCols);
	blit(buffer, screen, 0, 0, 0, 0, numRows, numCols);

}

int incMod(int num, int inc, int dim){
	return (num+inc+dim)%dim;
}

void transitionFunction(int row, int col) {
	writeMatrix[row][col]=  readMatrix[row][incMod(col,dir,numCols)];

}

void executeTransitionFunction() {
	for (int row = 0; row < numRows; row++)
		for (int col = 0; col < numCols; col++)
			transitionFunction(row, col);
}

void swapMatrices() {
	int **p = readMatrix;
	readMatrix = writeMatrix;
//	for (int row = 0; row < numRows; row++)
//			for (int col = 0; col < numCols; col++)
//				readMatrix[row][col]=writeMatrix[row][col];
	writeMatrix = p;
}

void controlLoop() {
	readConfigFile();
	initAllegro();

	initModel();
	for (int step = 0; step < numStep; ++step) {
		drawWithAllegro(step);
		readkey();
		if (key[KEY_LEFT])
			dir = 1;
		if (key[KEY_RIGHT])
			dir = -1;
		executeTransitionFunction();
		swapMatrices();

	}

	readkey();
}

int main() {
	controlLoop();
	return 0;
}
END_OF_MAIN()
