#include <iostream>
#include <thread>
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

        std::vector<std::thread> threads(4);
        threads[0] = std::thread(&SiliconOxideMatrix::evolution, &matrix);
        threads[1] = std::thread(&SiliconOxideMatrix::evolution, &matrix);
        threads[2] = std::thread(&SiliconOxideMatrix::evolution, &matrix);
        threads[3] = std::thread(&SiliconOxideMatrix::evolution, &matrix);

        auto start = std::chrono::high_resolution_clock::now();
        threads[0].join();
        threads[1].join();
        threads[2].join();
        threads[3].join();

        auto stop = std::chrono::high_resolution_clock::now();

        matrix.printMatrixToFile("map_evo.txt");

	
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

        std::cout << "Evolution took " << duration.count() << " seconds" << std::endl;
        matrix.debug_getNumber();

        //matrix.printMatrixToFile("map.txt");
        //matrix.printMatrixToImage("map.ppm");
        //matrix.initializeArrayO();
    }
    else
    {
        std::cout << "Invalid input percentage: " << percentage << std::endl;
    }
    
}
