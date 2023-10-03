// SiliconOxideMatrix.cpp
#include "SiliconOxideMatrix.h"
#include <iostream>
#include <map>
#include <cmath>
#include <random>

SiliconOxideMatrix::SiliconOxideMatrix()
{
	maxi = 10;
	maxj = 10;
    generalNumberOfAtoms = maxi * maxj;
    numberOfOxygenAtoms = 0; // zero by default
	data.resize(maxi, std::vector<double>(maxj, 1)); // Initialize the data vector
}

SiliconOxideMatrix::~SiliconOxideMatrix()
{
	// at the moment empty
}

void SiliconOxideMatrix::fillMatrix()
{	
    for (int i = 0; i < maxi; ++i)
    {
        for (int j = 0; j < maxj; ++j)
        {
            // according to the mathematics theory - array starts from index 1, so decide to adapt C++ array indexing to mathematics
            if ((1 == i % 2) && (1 == j % 2)) // Kremniy places
            {
                data[i][j] = 0.0;
            }
            else if ((0 == i % 2) && (0 == j % 2)) // ??? N/A places
            {
                data[i][j] = 9.0;
            }
            else if (((1 == i % 2) && (0 == j % 2)) || ((0 == i % 2) && (1 == j % 2))) // Oxygen places
            {
                allPossibleOxygenPlaces.push_back({i, j});
            }
            else
            {
                // empty at the moment
            }
        }
    } 
}

void SiliconOxideMatrix::printMatrixToFile(const std::string& fileName)
{
    std::ofstream outputFile(fileName); // Create an output file stream

    if (outputFile.is_open())
    {
        for (int i = 0; i < maxi; i++)
        {
            for (int j = 0; j < maxj; j++)
            {
                outputFile << data[i][j] << " "; // Write data to the file
            }
            outputFile << std::endl; // Add a newline at the end of each row
        }
        outputFile.close(); // Close the file
        std::cout << "Matrix has been written to " << fileName << std::endl;
    }
    else
    {
        std::cerr << "Unable to open the file " << fileName << std::endl;
    }
}

void SiliconOxideMatrix::calculateNumberOfOxygenAtomsByInputPercentage(const int percentage)
{
    numberOfOxygenAtoms = generalNumberOfAtoms / 2 * (percentage / 100.0);
    std::cout << "generalNumberOfAtoms: " << generalNumberOfAtoms << std::endl;
    std::cout << "numberOfOxygenAtoms: " << numberOfOxygenAtoms << std::endl;
    std::cout << "allPossibleOxygenPlaces: " << allPossibleOxygenPlaces.size() << std::endl;
}

void SiliconOxideMatrix::putOxygenAtomsOnRdmPlaces()
{
    std::vector<std::pair<int, int>> oxygenPlacesBuff = allPossibleOxygenPlaces;
    std::random_device rd;
    std::mt19937 mercenne_algorith (rd());

    while (0 < numberOfOxygenAtoms)
    {
        std::uniform_real_distribution<double> dist(0, oxygenPlacesBuff.size());
        int rnd_index = dist(mercenne_algorith);

        std::cout << "remainOxygenAtomsToPlace: " << numberOfOxygenAtoms << " rnd_index: " << rnd_index << " possition: " << oxygenPlacesBuff[rnd_index].first << ", " << oxygenPlacesBuff[rnd_index].second << std::endl;
        data[oxygenPlacesBuff[rnd_index].first][oxygenPlacesBuff[rnd_index].second] = 2;
        allOxygenPossitions.push_back({true, oxygenPlacesBuff[rnd_index].first, oxygenPlacesBuff[rnd_index].second});
        oxygenPlacesBuff.erase(oxygenPlacesBuff.begin() + rnd_index);

        // todo Critical section for multithreading
        --numberOfOxygenAtoms;
    }
}

void SiliconOxideMatrix::debug_getNumber()
{
    int count = 0;
    for (int i = 0; i < maxi; ++i)
    {
        for (int j = 0; j < maxj; ++j)
        {
            if (data[i][j] == 2)
            {
                ++count;
            }
        }
    }
    std::cout << "count: " << count << std::endl;
}

// Method to calculate the number of oxygen atoms and handle boundary transitions
int SiliconOxideMatrix::calculationNumOxigen(int i, int j) { 
    int value = 0;
    // Calculate the number of adjacent cells with value -1
    if (data[i - 1][j] == -1) {
        ++value;
    }   
    if (data[i + 1][j] == -1) {
        ++value;
    }
    if (data[i][j - 1] == -1) {
        ++value;
    }
    if (data[i][j + 1] == -1) {
        ++value;
    }
    return value; // Return the calculated value
}

void SiliconOxideMatrix::printMatrixToImage(const std::string& fileName) {
    std::ofstream outputFile(fileName);
    
    if (outputFile.is_open()) {
        outputFile << "P3\n" << maxj << " " << maxi << "\n255\n";
        // Define a map to associate values with colors
        std::map<int, std::string> colorMap;
        colorMap[0] = "0 0 0 ";        // Black
        colorMap[1] = "255 0 0 ";      // Red
        colorMap[2] = "0 0 255 ";      // Blue
        colorMap[3] = "0 255 0 ";      // Green
        colorMap[4] = "255 255 0 ";    // Yellow

        for (int i = 0; i < maxi - 1; i++) {
            for (int j = 0; j < maxj - 1; j++) {
                if ((i + 1) % 2 && (j + 1) % 2) { // [even][even] position calculation
                    unsigned int value = calculationNumOxigen(i, j); // Number of oxygen atoms calculation
                    // Use the color map to set the color based on the value
                    if (colorMap.find(value) != colorMap.end()) {
                        outputFile << colorMap[value];
                    } else {
                        outputFile << "255 255 255 "; // Default to white for unknown values
                    }
                }
            }
            outputFile << "\n";
        }
        outputFile.close();
        std::cout << "Matrix has been saved as image map.ppm " << fileName << std::endl;
    } else {
        std::cerr << "Unable to save the image map.ppm" << fileName << std::endl;
    }
}

double SiliconOxideMatrix::randomGenerator(unsigned int first_interval, unsigned int last_interval) {
    std::random_device rd; // Create a random number generator engine
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(first_interval, last_interval); // Create a uniform distribution in the range [first_interval, last_interval]
    return dist(gen);
}

bool SiliconOxideMatrix::metropolisCondition(double a, double b) {
	
	double randomValue = randomGenerator(0, 1); // Generate a random number in the range [0, 1]
	double metropolisValue = exp(-((b - a) / kT_eV)); // Calculate the Metropolis condition value

    // Compare the random value with the Metropolis value
    if (randomValue > metropolisValue) {
        return false; // Reject the move
    } else {
        return true; // Accept the move
    }
}

bool SiliconOxideMatrix::isOxygenInCell(int i, int j) {
    if ((i % 2 == 0 && j % 2 == 1) || (i % 2 == 1 && j % 2 == 0)) {
        // Check if it's an [even][odd] or [odd][even] cell
        if (data[i][j] == -1) {
            return true; // Element in [even][odd] or [odd][even] cell is -1
        }
    }
    return false; // Element is not in [even][odd] or [odd][even] cell or is not -1
}

bool SiliconOxideMatrix::findingOxigeninCeil(int &index_i, int & index_j) {
    while (!isOxygenInCell(index_i, index_j)) {
    	index_i = rand() % maxi;
    	index_j = rand() % maxj;
	}
	return true;
}

void SiliconOxideMatrix::evolution(void) {
    unsigned int iteration = 1000000000;

    while (iteration--) {

        int index_i = rand() % maxi;
        int index_j = rand() % maxj;

		findingOxigeninCeil(index_i, index_j);
        
        // Process the cell with data[index_i][index_j]
        // ...
    }
}

double SiliconOxideMatrix::Penalty(int i, int j) {
    return Delta[calculationNumOxigen(i, j)];
}



// Oxigen init procedurec
void SiliconOxideMatrix::updateDataMatrix(std::vector<std::pair<int, int>>& ArrayO) {
    int index = rand() % ArrayO.size();
    data[ArrayO[index].first][ArrayO[index].second] = 1;
    ArrayO.erase(ArrayO.begin() + index);
}

void SiliconOxideMatrix::updateArrayOLoop(std::vector<std::pair<int, int>>& ArrayO) {
    
	int count = 0;   
    while (count < maxIterations && !ArrayO.empty()) {
        updateDataMatrix(ArrayO);
        count++;
    }
}

void SiliconOxideMatrix::initializeArrayO() {	
    // Create ArrayO by iterating over coordinates
    for (int n = 1; n <= maxi / 2; n++) {
        for (int m = 1; m <= maxj / 2; m++) {
            ArrayO.push_back({2 * n - 1, 2 * m});
            ArrayO.push_back({2 * n, 2 * m - 1});
        }
    }
    updateArrayOLoop(ArrayO);
}


