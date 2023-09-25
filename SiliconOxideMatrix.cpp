// SiliconOxideMatrix.cpp
#include "SiliconOxideMatrix.h"
#include <iostream>

SiliconOxideMatrix::SiliconOxideMatrix() {
	
}

SiliconOxideMatrix::~SiliconOxideMatrix() {
	
}

SiliconOxideMatrix::SiliconOxideMatrix(int numRows, int numCols) : maxi(numRows), maxj(numCols) {
    data.resize(maxi, std::vector<double>(maxj, 0.0)); // Initialize the data vector with zeros
}

void SiliconOxideMatrix::fillMatrix() {
    for (int i = 0; i < maxi; i++) {
        for (int j = 0; j < maxj; j++) {
            data[i][j] = 0.0;
        }
    }
}

void SiliconOxideMatrix::printMatrixToFile(const std::string& fileName) {
	
    std::ofstream outputFile(fileName); // Create an output file stream

    if (outputFile.is_open()) {
        for (int i = 0; i < maxi; i++) {
            for (int j = 0; j < maxj; j++) {
                outputFile << data[i][j] << " "; // Write data to the file
            }
            outputFile << std::endl; // Add a newline at the end of each row
        }
        outputFile.close(); // Close the file
        std::cout << "Matrix has been written to " << fileName << std::endl;
    } else {
        std::cerr << "Unable to open the file " << fileName << std::endl;
    }
}

void SiliconOxideMatrix::printMatrixToImage(const std::string& fileName) {
    std::ofstream outputFile(fileName);

    if (outputFile.is_open()) {
        outputFile << "P3\n" << maxj << " " << maxi << "\n255\n";

        for (int i = 0; i < maxi; i++) {
            for (int j = 0; j < maxj; j++) {
                double value = data[i][j];

                if (value == 0.0) {
                    // Set red color
                    outputFile << "255 0 0 ";
                } else {
                    // Set other color (e.g., white)
                    outputFile << "255 255 255 ";
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

double SiliconOxideMatrix::Penalty(int iota, int phi) {
    int x = static_cast<int>(data[iota - 1][phi]) +
            static_cast<int>(data[iota][((phi + 1) % maxj)]) +
            static_cast<int>(data[(iota + 1) % maxi][phi]) +
            static_cast<int>(data[iota][phi - 1]);

    return Delta[x + 1];
}

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


