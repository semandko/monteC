#include "static.h"

int main(int argc, char** argv) {
    // Output the value of matrix dimension n
    std::cout << "matrix n = " << n << std::endl;

    // Resize the grid matrix to dimensions n x n and initialize it with zeros
    grid.resize(n, std::vector<int>(n, 1));

    // Fill the grid matrix with values using the fillMatrix function
    fillMatrix();

    // Prompt the user to input the percentage of Oxygen (O) atoms in the matrix
    int percentO;
    std::cout << "Input percentage of Oxygen: ";
    std::cin >> percentO;

    // Calculate the number of Oxygen (O) atoms based on the percentage
    numberO = (n * n / 2) * percentO / 100;

    // Record the start time for measuring the duration of the fillOxigen function
    auto start = std::chrono::high_resolution_clock::now();

    // Fill the grid with Oxygen (O) atoms using the fillOxigen function
    fillOxigen();

    // Record the stop time to calculate the duration of the fillOxigen function
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // Output the duration of the fillOxigen function execution
    std::cout << "Time taken by fillOxigen: " << duration.count() << " microseconds" << std::endl;

    // Generate random coordinates i and j within the matrix dimensions
    int i = rand() % n;
    int j = rand() % n;

    // Output the randomly selected coordinates i and j
    std::cout << "i = " << i ;
    std::cout << "; j = " << j << std::endl;

    // Calculate and output the number of Oxygen (O) atoms surrounding the selected coordinates
    std::cout << "Number of Oxygen atoms: " << countSurroundingOs(i, j) << std::endl;

    // Output the entire grid matrix to the console
    for (const std::vector<int>& row : grid) {
        for (int value : row) {
            std::cout << value << ' ';
        }
        std::cout << std::endl; // Move to a new line after each row of the matrix
    }

    // Open an output file named "grid.txt" for writing
    std::ofstream outputFile("grid.txt");

    // Check if the file was successfully opened
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open output file!" << std::endl;
        return 1; // Return an error code
    }

    // Write the entire grid matrix to the "grid.txt" file
    for (const std::vector<int>& row : grid) {
        for (int value : row) {
            outputFile << value << ' ';
        }
        outputFile << std::endl; // Move to a new line in the file after each row of the matrix
    }

    // Close the output file
    outputFile.close();

    // Wait for user input before exiting the program
    std::cin.get();

    return 0; // Return a success code
}

