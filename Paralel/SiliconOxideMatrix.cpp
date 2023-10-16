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
    evo_iterations = 10000000;
	data.resize(N, std::vector<MatrixPossition>(N, {1, true})); // Initialize the data vector

    // horizontalKremniyPlacesOfset =  {{0, -1}, {0, -3}, {-2, -1}, {-2, 1}, {0, 1}, {0, 3}, {2, 1}, {2, -1}};
    horizontalOxygenPlacesOfset =   {{0, -2}, {-1, -1}, {-1, 1}, {0, 2}, {1, 1}, {1, -1}};
    // verticalKremniyPlacesOfset =    {{-1, -2}, {-1, 0}, {-3, 0}, {-1, 2}, {1, 2}, {1, 0}, {3, 0}, {1, -2}};
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
                allPossibleOxygenPlaces.push_back({i, j, true});
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
    std::cout << "Iterations number: " << evo_iterations << std::endl;
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
    std::vector<OxygenPos> oxygenPlacesBuff = allPossibleOxygenPlaces;

    int buff = numberOfOxygenAtoms;
    while (0 < buff)
    {
        int rnd_index = getRdmIntNumber(0, oxygenPlacesBuff.size());

        // std::cout << "remainOxygenAtomsToPlace: " << numberOfOxygenAtoms << " rnd_index: " << rnd_index << " possition: " << oxygenPlacesBuff[rnd_index].first << ", " << oxygenPlacesBuff[rnd_index].second << std::endl;
        data[oxygenPlacesBuff[rnd_index].row][oxygenPlacesBuff[rnd_index].column].value = 2;
        allOxygenPossitions.push_back({oxygenPlacesBuff[rnd_index].row, oxygenPlacesBuff[rnd_index].column, true});
        oxygenPlacesBuff.erase(oxygenPlacesBuff.begin() + rnd_index);

        --buff;
    }
}

bool SiliconOxideMatrix::checkingCellOccupationAndJumping(int row, int column, bool method)
{
    bool placesavailability = false;
    std::vector<std::pair<int, int>> OxygenPlacesOfset;

    if (method)
    {
        OxygenPlacesOfset = horizontalOxygenPlacesOfset;
    }
    else
    {
        OxygenPlacesOfset = verticalOxygenPlacesOfset;
    }

    // check Oxygen places occupation and places availability
    for (const auto OxygenPlaceOfset : OxygenPlacesOfset)
    {
        int indexI = (row + OxygenPlaceOfset.first + N) % N;
        int indexJ = (column + OxygenPlaceOfset.second + N) % N;

        if (!data[indexI][indexJ].valid)
        {
            return false;
        }

        if (data[indexI][indexJ].value == 1.0)
        {
            placesavailability = true;
        }
    }

    return placesavailability;
}

std::vector<std::pair<int, int>> SiliconOxideMatrix::getJumpPossitionsAndChangeValidationToNV(int row, int column, bool method)
{
    std::vector<std::pair<int, int>> OxygenPlacesOfset;
    std::vector<std::pair<int, int>> possitionsToJump;

    if (method)
    {
        OxygenPlacesOfset = horizontalOxygenPlacesOfset;
    }
    else
    {
        OxygenPlacesOfset = verticalOxygenPlacesOfset;
    }

    data[row][column].valid = false;
    for (const auto placeOfset : OxygenPlacesOfset)
    {
        int buffRow = (row + placeOfset.first + N) % N;
        int buffColumn = (column + placeOfset.second + N) % N;
        if (data[buffRow][buffColumn].value == 1.0)
        {
            possitionsToJump.push_back({buffRow, buffColumn});
        }
        data[buffRow][buffColumn].valid = false;
    }

    return possitionsToJump;
}

void SiliconOxideMatrix::changeValidationToV(int row, int column, bool method)
{
    std::vector<std::pair<int, int>> OxygenPlacesOfset;

    if (method)
    {
        OxygenPlacesOfset = horizontalOxygenPlacesOfset;
    }
    else
    {
        OxygenPlacesOfset = verticalOxygenPlacesOfset;
    }

    data[row][column].valid = true;
    for (const auto ofsetValue : OxygenPlacesOfset)
    {
        int buffRow = (row + ofsetValue.first + N) % N;
        int buffColumn = (column + ofsetValue.second + N) % N;
        data[buffRow][buffColumn].valid = true;
    }
}

void SiliconOxideMatrix::evolution()
{
    while (evo_iterations > 0)
    {   
        // std::cout << "Remain number of iteration: " << evo_iterations << std::endl;
        int oxygenRndPossitionIndex = getRdmIntNumber(0, allOxygenPossitions.size());

        valid_m.lock();
        if (allOxygenPossitions[oxygenRndPossitionIndex].valid)
        {
            allOxygenPossitions[oxygenRndPossitionIndex].valid = false;
            valid_m.unlock();

            int oxygenRndRow = allOxygenPossitions[oxygenRndPossitionIndex].row;
            int oxygenRndColumn = allOxygenPossitions[oxygenRndPossitionIndex].column;
            bool method = (oxygenRndRow % 2 == 1); // horizontal == 1, vertical == 0

            // std::cout << "RND oxygen row: " << oxygenRndRow << ", column: " << oxygenRndColumn << std::endl;

            valid_m.lock();
            if (checkingCellOccupationAndJumping(oxygenRndRow, oxygenRndColumn, method))
            {
                valid_m.unlock();

                // std::cout << "Found available places to jump" << std::endl;

                std::vector<std::pair<int, int>> possitionsToJump;
                int toJumpPossitionIndex, toJumpRow, toJumpColumn;
                double A = 0;
                double nA = 0;

                valid_m.lock();
                possitionsToJump = getJumpPossitionsAndChangeValidationToNV(oxygenRndRow, oxygenRndColumn, method);
                valid_m.unlock();

                if (method)
                {
                    A = calculatePenaltyForKremniy(oxygenRndRow, (oxygenRndColumn - 1 + N) % N) + calculatePenaltyForKremniy(oxygenRndRow, (oxygenRndColumn + 1 + N) % N);
                }
                else
                {
                    A = calculatePenaltyForKremniy((oxygenRndRow - 1 + N) % N, oxygenRndColumn) + calculatePenaltyForKremniy((oxygenRndRow + 1 + N) % N, oxygenRndColumn);
                }

                // std::cout << "Energy A: " << A << std::endl;

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
                    valid_m.lock();
                    allOxygenPossitions[oxygenRndPossitionIndex].row = toJumpRow;
                    allOxygenPossitions[oxygenRndPossitionIndex].column = toJumpColumn;
                    data[toJumpRow][toJumpColumn].value = 2.0;
                    data[oxygenRndRow][oxygenRndColumn].value = 1.0;
                    valid_m.unlock();

                    // std::cout << "JUMPED!" << std::endl;
                }
                else
                {
                    // std::cout << "nA >= A, trying metropolis condition" << std::endl;

                    if (metropolisCondition(A, nA))
                    {
                        valid_m.lock();
                        allOxygenPossitions[oxygenRndPossitionIndex].row = toJumpRow;
                        allOxygenPossitions[oxygenRndPossitionIndex].column = toJumpColumn;
                        data[toJumpRow][toJumpColumn].value = 2.0;
                        data[oxygenRndRow][oxygenRndColumn].value = 1.0;
                        valid_m.unlock();

                        // std::cout << "JUMPED!" << std::endl;
                    }
                    else
                    {
                        // std::cout << "Metropolis condition to aged" << std::endl;
                    }
                }

                valid_m.lock();
                changeValidationToV(oxygenRndRow, oxygenRndColumn, method);
                valid_m.unlock();
            }
            else
            {
                valid_m.unlock();
                // std::cout << "Ofset position(s) not valid or no available places" << std::endl;
            }

            valid_m.lock();
            allOxygenPossitions[oxygenRndPossitionIndex].valid = true;
            valid_m.unlock();


            // valid_m.lock();
            // if (debug_getNumber_test() <= 2)
            // {   
            //     std::cout << "Something went wrong:" << std::endl;
            //     std::cout << "  Oxygen number:" << debug_getNumber_test() << std::endl;
            //     std::cout << "  Remain iterations: " << evo_iterations << std::endl;
            //     std::cout << "  Row: " << oxygenRndRow << ", Column: " << oxygenRndColumn << std::endl;

            //     for (int i = 0; i < allOxygenPossitions.size(); ++i)
            //     {
            //         std::cout << "  Oxygen number: " << i << ", Row: " << allOxygenPossitions[i].row << ", Column: " << allOxygenPossitions[i].column << std::endl;
            //     }
                
            //     printMatrixToFile(std::to_string(evo_iterations));

            //     numberOfOxygenAtoms = debug_getNumber_test();
            //     exit(0);
            // }
            // valid_m.unlock();
        }
        else
        {
            valid_m.unlock();
            // std::cout << "Oxygen possition not valid" << std::endl;
        }
        
        // std::cout << "End iteration" << std::endl <<std::endl;

        //printMatrixToFile("map_evo.txt");
        --evo_iterations;
    }
}

int SiliconOxideMatrix::debug_getNumber_test()
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

    return count;
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


