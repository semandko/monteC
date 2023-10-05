#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <map>
#include <cmath>
#include <random>


#define _CONST_A 1.38e-23
#define _CONST_B 1.602e-19
#define _K_B (_CONST_A / _CONST_B)
#define _X_0 0.1
#define _T 1000
#define _kT_EV(_T) (1.0 / 11604 * (_T + 273.15)) // 0.11 eV for 1000 C

// typedefs
typedef struct {
    int i;
    int j;
} tcell;

typedef struct {
    tcell cell;
    int value;
    bool valid;
} tRomb;

// variables
int n = 6;
int numberO;

std::vector<std::vector<tRomb>> grid(n, std::vector<tRomb>(n, { {0, 0}, 1, 1 }));
std::vector<double> Delta{0.0, 0.5, 0.51, 0.22, 0.0}; // penalty energy
double penaltyValue;
double x0 = _X_0;
double kB = _K_B;
double kT_eV = _kT_EV(_T);
int T = _T;

// prototypes
void fillMatrix();
void fillOxigen();
int countSurroundingOs(int i, int j);
double randomGenerator(unsigned int first_interval, unsigned int last_interval);
double penalty(int i, int j);
bool metropolisCondition(double a, double b);
void evolution(void);

void testOcupationO(void);

bool outputFile(void);
void configurator(void);

// code section *************************************************
void fillMatrix() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((i % 2 == 1) && (j % 2 == 1)) {
                grid[i][j].value = 0.0; // Set value to 1 for [even][even] cells
            }
        }
    }
}

void fillOxigen() {
    while (numberO > 0) {
        int i = rand() % n;
        int j = rand() % n;
		
		// cells type O [even][odd] or [odd][even]
        if (((i % 2 == 0) && (j % 2 == 1)) || ((i % 2 == 1) && (j % 2 == 0))) {
            if (grid[i][j].value != 8) {
                grid[i][j].value = 8;
                numberO--;
            }
        }
    }
}

int countSurroundingOs(int i, int j) {
    int count = 0; // number of cells type O

    if (i % 2 == 1 && j % 2 == 1) { // [even][even] cell of type K
        int di[] = {-1, 1, 0, 0}; // vertical
        int dj[] = {0, 0, -1, 1}; // horisontal

        for (int dir = 0; dir < 4; dir++) {
            int ni = (i + di[dir]) % n;
            int nj = (j + dj[dir]) % n;

            if (grid[ni][nj].value == 8) {
                count++;
            }
        }
    } else {
    	std::cout << "not valid ";
	}

    return count;
}

double randomGenerator(unsigned int first_interval, unsigned int last_interval) {
    std::random_device rd; // Create a random number generator engine
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(first_interval, last_interval);
    return dist(gen);
}

double pPenalty(int i, int j) {
    return Delta[countSurroundingOs(i, j)];
}

bool metropolisCondition(double a, double b) {
	
	double randomValue = randomGenerator(0, 1); // Generate a random number in the range [0, 1]
	double metropolisValue = exp(-((b - a) / kT_eV)); // Calculate the Metropolis condition value

    // Compare the random value with the Metropolis value
    if (randomValue > metropolisValue) {
        return false; // Reject the jump
    } else {
        return true; // Accept the jump
    }
}

void evolution(void) {
	
}

// testing block start
void testOcupationO(void) {
    for (int k = 0; k < 5; ++k) {
        int i = rand() % n;
        int j = rand() % n;

        std::cout << k << ": [" << i << "][" << j << "] = " << countSurroundingOs(i, j) << " cells" << std::endl;
    }
}
// testing block stop

bool outputFile(void) {
    std::ofstream outputFile("grid.txt");

    // Check if the file was successfully opened
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open output file!" << std::endl;
        return 1; // Return an error code
    }

    // Write the entire grid matrix to the "grid.txt" file
    for (const std::vector<tRomb>& row : grid) {
        for (const tRomb& cell : row) {
            outputFile << cell.value << ' ';
        }
        outputFile << std::endl;
    }

    outputFile.close();
    return 0;
}

void printMatrixToImage(const std::string& fileName) {
    std::ofstream outputFile(fileName);
    
    if (outputFile.is_open()) {
        outputFile << "P3\n" << n << " " << n << "\n255\n";
        std::map<int, std::string> colorMap; // Define a map to associate values with colors
        colorMap[0] = "0 0 0 ";        // Black
        colorMap[1] = "255 0 0 ";      // Red
        colorMap[2] = "0 0 255 ";      // Blue
        colorMap[3] = "0 255 0 ";      // Green
        colorMap[4] = "255 255 0 ";    // Yellow

        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - 1; j++) {
                if ((i % 2 == 1) && (j % 2 == 1)) { // [even][even] position calculation
                    unsigned int value = countSurroundingOs(i, j); // Number of oxygen atoms calculation surround [even][even] cell
                    // Use the color map to set the color based on the value
                    if (colorMap.find(value) != colorMap.end()) {
                        outputFile << colorMap[value];
                    }
                }
            }
            outputFile << "\n";
        }
        outputFile.close();
        std::cout << "Matrix has been saved as image " << fileName << std::endl;
    } else {
        std::cerr << "Unable to save the image " << fileName << std::endl;
    }
}

void configurator(void) {

    std::cout << "Number of all elements = " << n << std::endl; // Output the value of matrix dimension n

    grid.resize(n, std::vector<tRomb>(n, { {0, 0}, 1, 1 })); // Resize the grid matrix to dimensions n x n and initialize it with 1s

    fillMatrix();

    int percentO; // Prompt the user to input the percentage of Oxygen (type O cell) in the matrix
    std::cout << "Input percentage of Oxygen: ";
    std::cin >> percentO;

    numberO = (n * n / 2) * percentO / 100; // Calculate the number of O cells
    std::cout << "Number of Oxygen cells = " << numberO << std::endl;

    auto start = std::chrono::high_resolution_clock::now(); // Record the start time for measuring the duration of the fillOxigen function

    fillOxigen(); // Fill the grid with O cells

    auto stop = std::chrono::high_resolution_clock::now(); // Record the stop time to calculate the duration of the fillOxigen function
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Time taken by fillOxigen = " << duration.count() << " microseconds" << std::endl;

}

//************************************************************************************************
int main(int argc, char** argv) {

    configurator();
    testOcupationO();

    outputFile();
    printMatrixToImage("map.ppm");

    std::cin.get(); // Wait for user input before exiting the program

    return 0; // Return a success code
}

