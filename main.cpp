#include <iostream>
#include "SiliconOxideMatrix.h"

/* run this program using the console pauser or add your own getch, system("pause") or input loop */

void configurator (void);

int main(int argc, char** argv) {
	
	configurator();
	 
    system("pause");
	return 0;
}

void configurator (void)
{
	int numRows = _NUM_ROWS;
    int numCols = _NUM_COLS;
	
	std::cout << "Matrix has numRows " << numRows << std::endl;
	std::cout << "Matrix has numCols " << numCols << std::endl;
	
    SiliconOxideMatrix matrix;

    matrix.fillMatrix();
    matrix.printMatrixToFile("map.txt");
    matrix.printMatrixToImage("map.ppm");
    //matrix.initializeArrayO();
}
