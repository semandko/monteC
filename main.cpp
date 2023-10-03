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
    int percentage = 0;;
	
	std::cout << "Matrix has numRows " << numRows << std::endl;
	std::cout << "Matrix has numCols " << numCols << std::endl;

    std::cout << "Please enter the integer Oxygen percentage (from 0 to 100): ";
    std::cin >> percentage;

    if (percentage >= 0 && percentage <= 100)
    {
        SiliconOxideMatrix matrix;

        matrix.fillMatrix();
        matrix.calculateNumberOfOxygenAtomsByInputPercentage(percentage);
        matrix.putOxygenAtomsOnRdmPlaces();
        matrix.debug_getNumber();
        matrix.printMatrixToFile("map.txt");
        //matrix.printMatrixToImage("map.ppm");
        //matrix.initializeArrayO();
    }
    else
    {
        std::cout << "Invalid input percentage: " << percentage << std::endl;
    }
    
}
