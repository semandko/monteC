#include <iostream>
#include <chrono>
#include "SiliconOxideMatrix.h"

/* run this program using the console pauser or add your own getch, system("pause") or input loop */

void configurator (void);

int main(int argc, char** argv)
{
	
	configurator();
	 
    system("pause");
	return 0;
}

void configurator (void)
{
    int percentage = 0;

    std::cout << "Please enter the integer Oxygen percentage (from 0 to 100): ";
    std::cin >> percentage;

    if (percentage >= 0 && percentage <= 100)
    {
        SiliconOxideMatrix matrix;

        matrix.fillMatrix();
        matrix.calculateNumberOfOxygenAtomsByInputPercentage(percentage);
        matrix.putOxygenAtomsOnRdmPlaces();
        matrix.printMatrixToFile("map_old.txt");

        auto start = std::chrono::high_resolution_clock::now();
        matrix.evolutionFindOxygen();
        auto stop = std::chrono::high_resolution_clock::now();

        matrix.printMatrixToFile("map_evo.txt");

	
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

        std::cout << "Evolution took " << duration.count() << " seconds" << std::endl;
        matrix.debug_getNumber();

        //matrix.debug_getNumber();
        //matrix.printMatrixToFile("map.txt");
        //matrix.printMatrixToImage("map.ppm");
        //matrix.initializeArrayO();
    }
    else
    {
        std::cout << "Invalid input percentage: " << percentage << std::endl;
    }
    
}
