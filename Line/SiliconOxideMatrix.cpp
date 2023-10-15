// SiliconOxideMatrix.cpp
#include "SiliconOxideMatrix.h"
#include <iostream>
#include <map>
#include <cmath>
#include <random>
#include <cfloat> // DBL_MAX

SiliconOxideMatrix::SiliconOxideMatrix()
{
    N = 500;
    generalNumberOfAtoms = N * N;
    numberOfOxygenAtoms = 0; // zero by default
	data.resize(N, std::vector<MatrixPossition>(N, {1, true})); // Initialize the data vector

    horizontalKremniyPlacesOfset =  {{0, -1}, {0, -3}, {-2, -1}, {-2, 1}, {0, 1}, {0, 3}, {2, 1}, {2, -1}};
    horizontalOxygenPlacesOfset =   {{0, -2}, {-1, -1}, {-1, 1}, {0, 2}, {1, 1}, {1, -1}};
    verticalKremniyPlacesOfset =    {{-1, -2}, {-1, 0}, {-3, 0}, {-1, 2}, {1, 2}, {1, 0}, {3, 0}, {1, -2}};
    verticalOxygenPlacesOfset =     {{-1, -1}, {-2, 0}, {-1, 1}, {1, 1}, {2, 0}, {1, -1}};
}

SiliconOxideMatrix::~SiliconOxideMatrix()
{
	// at the moment empty
}

void SiliconOxideMatrix::fillMatrix()
{	
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            // according to the mathematics theory - array starts from index 1, so decide to adapt C++ array indexing to mathematics
            if ((1 == i % 2) && (1 == j % 2)) // Kremniy places
            {
                data[i][j].value = 0.0;
            }
            else if ((0 == i % 2) && (0 == j % 2)) // N/A places
            {
                data[i][j].value = 9.0;
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
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                outputFile << data[i][j].value << " "; // Write data to the file
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

int SiliconOxideMatrix::getRdmIntNumber(int start, int notIncluydedEnd)
{
    std::random_device rd;
    std::mt19937 mercenne_algorith (rd());
    std::uniform_real_distribution<double> dist(start, notIncluydedEnd);

    return dist(mercenne_algorith);
}

double SiliconOxideMatrix::getRdmDoubleNumber(int start, int notIncludedEnd)
{
    std::random_device rd;
    std::mt19937 mercenne_algorith (rd());
    std::uniform_real_distribution<> dist(start, std::nextafter(notIncludedEnd, DBL_MAX));

    return dist(mercenne_algorith);
}

void SiliconOxideMatrix::putOxygenAtomsOnRdmPlaces()
{
    std::vector<std::pair<int, int>> oxygenPlacesBuff = allPossibleOxygenPlaces;

    while (0 < numberOfOxygenAtoms)
    {
        int rnd_index = getRdmIntNumber(0, oxygenPlacesBuff.size());

        // std::cout << "remainOxygenAtomsToPlace: " << numberOfOxygenAtoms << " rnd_index: " << rnd_index << " possition: " << oxygenPlacesBuff[rnd_index].first << ", " << oxygenPlacesBuff[rnd_index].second << std::endl;
        data[oxygenPlacesBuff[rnd_index].first][oxygenPlacesBuff[rnd_index].second].value = 2;
        allOxygenPossitions.push_back({oxygenPlacesBuff[rnd_index].first, oxygenPlacesBuff[rnd_index].second});
        oxygenPlacesBuff.erase(oxygenPlacesBuff.begin() + rnd_index);

        // todo Critical section for multithreading
        --numberOfOxygenAtoms;
    }
}

bool SiliconOxideMatrix::checkingHorizontalCellOccupationAndJumping(int row, int column)
{
    bool placesavailability = false;

    // check Oxygen places occupation and places availability
    for (const auto OxygenPlaceOfset : horizontalOxygenPlacesOfset)
    {
        int indexI = (row + OxygenPlaceOfset.first + N) % N;
        int indexJ = (column + OxygenPlaceOfset.second + N) % N;

        // no need at the moment
        // if (!data[indexI][indexJ].valid)
        // {
        //     return false;
        // }

        if (data[indexI][indexJ].value == 1.0)
        {
            placesavailability = true;
            break; // to be removed
        }
    }

    return placesavailability;
}

bool SiliconOxideMatrix::checkingVerticalCellOccupationAndJumping(int row, int column)
{
    bool placesavailability = false;

    // check Oxygen places occupation and places availability
    for (const auto OxygenPlaceOfset : verticalOxygenPlacesOfset)
    {
        int indexI = (row + OxygenPlaceOfset.first + N) % N;
        int indexJ = (column + OxygenPlaceOfset.second + N) % N;

        // no need at the moment
        // if (!data[indexI][indexJ].valid)
        // {
        //     return false;
        // }

        if (data[indexI][indexJ].value == 1.0)
        {
            placesavailability = true;
            break; // to be removed
        }
    }

    return placesavailability;
}

void SiliconOxideMatrix::evolution()
{
    unsigned long long n = 10000000;
    while (n--)
    {
        // std::cout << "Remain number of iteration: " << n << std::endl;
        int oxygenRndPossitionIndex = getRdmIntNumber(0, allOxygenPossitions.size());
        int oxygenRndRow = allOxygenPossitions[oxygenRndPossitionIndex].first;
        int oxygenRndColumn = allOxygenPossitions[oxygenRndPossitionIndex].second;

        // std::cout << "RND oxygen row: " << oxygenRndRow << ", column: " << oxygenRndColumn << std::endl;

        // no need at the moment
        // if (data[oxygenRndRow][oxygenRndColumn].valid)
        // {
            if (oxygenRndRow % 2 == 1)
            {
                // std::cout << "Horizontal type" << std::endl;
                bool ret = checkingHorizontalCellOccupationAndJumping(oxygenRndRow, oxygenRndColumn);

                if (ret)
                {
                    // std::cout << "Found available places to jump" << std::endl;

                    int toJumpPossitionIndex, toJumpRow, toJumpColumn;
                    std::vector<std::pair<int, int>> possitionsToJump;
                    double A = calculatePenaltyForKremniy(oxygenRndRow, (oxygenRndColumn - 1 + N) % N) + calculatePenaltyForKremniy(oxygenRndRow, (oxygenRndColumn + 1 + N) % N);
                    double nA = 0;

                    // std::cout << "Energy A: " << A << std::endl;

                    for (const auto placeOfset : horizontalOxygenPlacesOfset)
                    {
                        int buffRow = (oxygenRndRow + placeOfset.first + N) % N;
                        int buffColumn = (oxygenRndColumn + placeOfset.second + N) % N;
                        if (data[buffRow][buffColumn].value == 1.0)
                        {
                            possitionsToJump.push_back({buffRow, buffColumn});
                        }
                    }

                    toJumpPossitionIndex = getRdmIntNumber(0, possitionsToJump.size());
                    toJumpRow = possitionsToJump[toJumpPossitionIndex].first;
                    toJumpColumn = possitionsToJump[toJumpPossitionIndex].second;

                    // std::cout << "RND toJumpPossition row: " << toJumpRow << ", column: " << toJumpColumn << std::endl;

                    if (data[(toJumpRow - 1 + N) % N][toJumpColumn].value == 0.0)
                    {
                        nA = calculatePenaltyForKremniy((toJumpRow - 1 + N) % N, oxygenRndColumn) + calculatePenaltyForKremniy((oxygenRndRow + 1 + N) % N, oxygenRndColumn);
                    }
                    else
                    {
                        nA = calculatePenaltyForKremniy(toJumpRow, (oxygenRndColumn - 1 + N) % N) + calculatePenaltyForKremniy(oxygenRndRow, (oxygenRndColumn + 1 + N) % N);
                    }

                    // std::cout << "Energy nA: " << nA << std::endl;

                    if (nA < A)
                    {
                        allOxygenPossitions[oxygenRndPossitionIndex].first = toJumpRow;
                        allOxygenPossitions[oxygenRndPossitionIndex].second = toJumpColumn;
                        data[toJumpRow][toJumpColumn].value = 2.0;
                        data[oxygenRndRow][oxygenRndColumn].value = 1.0;
                        // std::cout << "JUMPED!" << std::endl;
                    }
                    else
                    {
                        // std::cout << "nA >= A, trying metropolis condition" << std::endl;

                        if (metropolisCondition(A, nA))
                        {
                            allOxygenPossitions[oxygenRndPossitionIndex].first = toJumpRow;
                            allOxygenPossitions[oxygenRndPossitionIndex].second = toJumpColumn;
                            data[toJumpRow][toJumpColumn].value = 2.0;
                            data[oxygenRndRow][oxygenRndColumn].value = 1.0;
                            // std::cout << "JUMPED!" << std::endl;
                        }
                        else
                        {
                            // std::cout << "Metropolis condition to aged" << std::endl;
                        }
                    }

                    // no need at the moment
                    // for (const auto ofsetValue : horizontalOxygenPlacesOfset)
                    // {
                    //     data[oxygenRndRow + ofsetValue.first][oxygenRndColumn + ofsetValue.second].valid = true;
                    // }
                }
                else
                {
                    // std::cout << "No places to jump" << std::endl;
                }
            }
            else
            {
                // std::cout << "Vertical type" << std::endl;
                bool ret = checkingVerticalCellOccupationAndJumping(oxygenRndRow, oxygenRndColumn);

                if (ret)
                {
                    // std::cout << "Found available places to jump" << std::endl;

                    int toJumpPossitionIndex, toJumpRow, toJumpColumn;
                    std::vector<std::pair<int, int>> possitionsToJump;
                    double A = calculatePenaltyForKremniy((oxygenRndRow - 1 + N) % N, oxygenRndColumn) + calculatePenaltyForKremniy((oxygenRndRow + 1 + N) % N, oxygenRndColumn);
                    double nA;

                    // std::cout << "Energy A: " << A << std::endl;

                    for (const auto placeOfset : verticalOxygenPlacesOfset)
                    {
                        int buffRow = (oxygenRndRow + placeOfset.first + N) % N;
                        int buffColumn = (oxygenRndColumn + placeOfset.second + N) % N;
                        if (data[buffRow][buffColumn].value == 1.0)
                        {
                            possitionsToJump.push_back({buffRow, buffColumn});
                        }
                    }

                    toJumpPossitionIndex = getRdmIntNumber(0, possitionsToJump.size());
                    toJumpRow = possitionsToJump[toJumpPossitionIndex].first;
                    toJumpColumn = possitionsToJump[toJumpPossitionIndex].second;

                    // std::cout << "RND toJumpPossition row: " << toJumpRow << ", column: " << toJumpColumn << std::endl;

                    if (data[(toJumpRow - 1 + N) % N][toJumpColumn].value == 0.0)
                    {
                        nA = calculatePenaltyForKremniy((toJumpRow - 1 + N) % N, oxygenRndColumn) + calculatePenaltyForKremniy((oxygenRndRow + 1 + N) % N, oxygenRndColumn);
                    }
                    else
                    {
                        nA = calculatePenaltyForKremniy(toJumpRow, (oxygenRndColumn - 1 + N) % N) + calculatePenaltyForKremniy(oxygenRndRow, (oxygenRndColumn + 1 + N) % N);
                    }

                    // std::cout << "Energy nA: " << nA << std::endl;

                    if (nA < A)
                    {
                        allOxygenPossitions[oxygenRndPossitionIndex].first = toJumpRow;
                        allOxygenPossitions[oxygenRndPossitionIndex].second = toJumpColumn;
                        data[toJumpRow][toJumpColumn].value = 2.0;
                        data[oxygenRndRow][oxygenRndColumn].value = 1.0;
                        // std::cout << "JUMPED!" << std::endl;
                    }
                    else
                    {
                        // std::cout << "nA >= A, trying metropolis condition" << std::endl;
                        if (metropolisCondition(A, nA))
                        {
                            allOxygenPossitions[oxygenRndPossitionIndex].first = toJumpRow;
                            allOxygenPossitions[oxygenRndPossitionIndex].second = toJumpColumn;
                            data[toJumpRow][toJumpColumn].value = 2.0;
                            data[oxygenRndRow][oxygenRndColumn].value = 1.0;
                            // std::cout << "JUMPED!" << std::endl;
                        }
                        else
                        {
                            // std::cout << "Metropolis condition not aged" << std::endl;
                        }
                    }

                    // no need at the moment
                    // for (const auto ofsetValue : horizontalOxygenPlacesOfset)
                    // {
                    //     data[oxygenRndRow + ofsetValue.first][oxygenRndColumn + ofsetValue.second].valid = true;
                    // }
                }
                else
                {
                    // std::cout << "No places to jump" << std::endl;
                }
            }
        // }
        // else
        // {
        //     // std::cout << "Possition is not valid";
        // }
        // std::cout << "End iteration" << std::endl <<std::endl;
    }
}

void SiliconOxideMatrix::debug_getNumber()
{
    int count = 0;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (data[i][j].value == 2)
            {
                ++count;
            }
        }
    }
    std::cout << "count: " << count << std::endl;
}

// Method to calculate the number of oxygen atoms and handle boundary transitions
int SiliconOxideMatrix::calculationNumberOfOxigens(int i, int j)
{ 
    int value = 0;
    // Calculate the number of adjacent cells with value -1
    if (data[(i - 1 + N) % N][j].value == 2) 
    {
        ++value;
    }   
    if (data[(i + 1 + N) % N][j].value == 2) 
    {
        ++value;
    }
    if (data[i][(j - 1 + N) % N].value == 2) 
    {
        ++value;
    }
    if (data[i][(j + 1 + N) % N].value == 2) 
    {
        ++value;
    }
    return value; // Return the calculated value
}

bool SiliconOxideMatrix::metropolisCondition(double a, double b)
{
	
	double randomValue = getRdmDoubleNumber(0, 1); // Generate a random number in the range [0, 1]
	double metropolisValue = exp(-((b - a) / kT_eV)); // Calculate the Metropolis condition value

    // std::cout << "randomValue: " << randomValue << " metropolisValue: " << metropolisValue << std::endl;
    // Compare the random value with the Metropolis value
    if (randomValue > metropolisValue)
    {
        return false; // Reject the move
    }
    else
    {
        return true; // Accept the move
    }
}

double SiliconOxideMatrix::calculatePenaltyForKremniy(int i, int j) 
{
    return Delta[calculationNumberOfOxigens(i, j)];
}

// void SiliconOxideMatrix::printMatrixToImage(const std::string& fileName) {
//     std::ofstream outputFile(fileName);
    
//     if (outputFile.is_open()) {
//         outputFile << "P3\n" << maxj << " " << maxi << "\n255\n";
//         // Define a map to associate values with colors
//         std::map<int, std::string> colorMap;
//         colorMap[0] = "0 0 0 ";        // Black
//         colorMap[1] = "255 0 0 ";      // Red
//         colorMap[2] = "0 0 255 ";      // Blue
//         colorMap[3] = "0 255 0 ";      // Green
//         colorMap[4] = "255 255 0 ";    // Yellow

//         for (int i = 0; i < maxi - 1; i++) {
//             for (int j = 0; j < maxj - 1; j++) {
//                  if ((i + 1) % 2 && (j + 1) % 2) { // [even][even] position calculation
                //     unsigned int value = calculationNumberOfOxigens(i, j); // Number of oxygen atoms calculation
                //     // Use the color map to set the color based on the value
                //     if (colorMap.find(value) != colorMap.end()) {
                //         outputFile << colorMap[value];
                //     } else {
                //         outputFile << "255 255 255 "; // Default to white for unknown values
                //     }
                // }
//             }
//             outputFile << "\n";
//         }
//         outputFile.close();
//         std::cout << "Matrix has been saved as image map.ppm " << fileName << std::endl;
//     } else {
//         std::cerr << "Unable to save the image map.ppm" << fileName << std::endl;
//     }
// }

// bool SiliconOxideMatrix::isOxygenInCell(int i, int j) {
//     if ((i % 2 == 0 && j % 2 == 1) || (i % 2 == 1 && j % 2 == 0)) {
//         // Check if it's an [even][odd] or [odd][even] cell
//         if (data[i][j] == -1) {
//             return true; // Element in [even][odd] or [odd][even] cell is -1
//         }
//     }
//     return false; // Element is not in [even][odd] or [odd][even] cell or is not -1
// }

// bool SiliconOxideMatrix::findingOxigeninCeil(int &index_i, int & index_j) {
//     while (!isOxygenInCell(index_i, index_j)) {
//     	index_i = rand() % maxi;
//     	index_j = rand() % maxj;
// 	}
// 	return true;
// }

// void SiliconOxideMatrix::evolution(void) {
//     unsigned int iteration = 1000000000;

//     while (iteration--) {

//         int index_i = rand() % maxi;
//         int index_j = rand() % maxj;

// 		findingOxigeninCeil(index_i, index_j);
        
//         // Process the cell with data[index_i][index_j]
//         // ...
//     }
// }

// Oxigen init procedurec
// void SiliconOxideMatrix::updateDataMatrix(std::vector<std::pair<int, int>>& ArrayO) {
//     int index = rand() % ArrayO.size();
//     data[ArrayO[index].first][ArrayO[index].second] = 1;
//     ArrayO.erase(ArrayO.begin() + index);
// }

// void SiliconOxideMatrix::updateArrayOLoop(std::vector<std::pair<int, int>>& ArrayO) {
    
// 	int count = 0;   
//     while (count < maxIterations && !ArrayO.empty()) {
//         updateDataMatrix(ArrayO);
//         count++;
//     }
// }

// void SiliconOxideMatrix::initializeArrayO() {	
//     // Create ArrayO by iterating over coordinates
//     for (int n = 1; n <= maxi / 2; n++) {
//         for (int m = 1; m <= maxj / 2; m++) {
//             ArrayO.push_back({2 * n - 1, 2 * m});
//             ArrayO.push_back({2 * n, 2 * m - 1});
//         }
//     }
//     updateArrayOLoop(ArrayO);
// }


