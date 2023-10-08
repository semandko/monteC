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

typedef struct {
	tcell cellTypeO;
	double penalty;
} tJumpCell;


// variables
int n = 100;
int numberO;

std::vector<std::vector<tRomb>> grid(n, std::vector<tRomb>(n, { {0, 0}, 1, 1 }));
std::vector<std::vector<tJumpCell>> jumpFreeSpaces; // 

std::vector<double> Delta{0.0, 0.5, 0.51, 0.22, 0.0}; // penalty energy

double penaltyValue;

double x0 = _X_0;
double kB = _K_B;
double kT_eV = _kT_EV(_T);
int T = _T;
int jumpingCounts;

// prototypes
void fillMatrix();
void fillOxigen();
int countSurroundingOs(int i, int j, bool isPresent);
void countSurroundingRomb(int i, int j);
double randomGenerator(unsigned int first_interval, unsigned int last_interval);
double penalty(int i, int j);
bool metropolisCondition(double oldPenaltySum, double newPenaltySum);
bool rombPenaltyCalculation(tJumpCell& jumpCell);
bool findOForJumping(tJumpCell& jumpCell);
bool initJumping(tJumpCell& jumpCell);
void jumping(tJumpCell& jumpCell);
void evolution(void);

bool printMatrixToTxt(void);
void printMatrixToImage(const std::string& fileName);

// void printMatrixToImagePNG(const std::string& fileName);

void configurator(void);

// testing block
void testOcupationO(void);
// testing block


// code section *************************************************
void fillMatrix() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((i % 2 == 1) && (j % 2 == 1)) {
                grid[i][j].value = 0.0; // Set value to 1 for [odd][odd] cells
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
    int countPresent = 0;

    if (i % 2 == 1 && j % 2 == 1) { // [odd][odd] cell of type K
        int di[] = {-1, 1, 0, 0}; // vertical
        int dj[] = {0, 0, -1, 1}; // horizontal

        for (int dir = 0; dir < 4; dir++) {
            int ni = (i + di[dir] + n) % n; // type O index
            int nj = (j + dj[dir] + n) % n; // type O index

        	// is type O presented
            if (grid[ni][nj].value == 8) {
                countPresent++;
            }
        }
    }
	else {
        std::cout << "not valid " << std::endl;
    }

    return countPresent;
}

void countSurroundingRomb(int i, int j) {
    if (i % 2 == 1 && j % 2 == 1) { // [odd][odd] cell of type K
        int di[] = {-1, 1, 0, 0}; // vertical
        int dj[] = {0, 0, -1, 1}; // horizontal

        for (int dir = 0; dir < 4; dir++) {
            int ni = (i + di[dir] + n) % n; // type O index
            int nj = (j + dj[dir] + n) % n; // type O index

            // Check if it's a free space
            if (grid[ni][nj].value != 8) {
            	
            	// structure tJumpCell
                tJumpCell jumpCell;
                jumpCell.cellTypeO.i = ni;
                jumpCell.cellTypeO.j = nj;
                jumpCell.penalty = 0.0;
                
                // vector of structure tJumpCell
                std::vector<tJumpCell> jumpCellVector; // Create a vector to hold the jumpCell
                jumpCellVector.push_back(jumpCell);
                // vector of vectors of structure tJumpCell
                jumpFreeSpaces.push_back(jumpCellVector); // Push the vector containing jumpCell
            }
        }
    } 
	else {
        std::cout << "not valid " << std::endl;
    }
}

bool rombPenaltyCalculation(tJumpCell& jumpCell) {
    bool ret = false;
    int i = jumpCell.cellTypeO.i; // for type O
    int j = jumpCell.cellTypeO.j; // for type O
	int i_a, j_a;
	int i_b, j_b;
	
    if (i % 2 == 1) { // Index for type O
        // Horizontal case for cells type K
		i_a = i; j_a = (j - 1 + n) % n; // left
		i_b = i; j_b = (j + 1) % n; // right
    }
	else {
        // Vertical case for cells type K 
        i_a = (i - 1 + n) % n; 	j_a = j; // up
		i_b = (i + 1) % n; 		j_b = j; // down
    }
    
    countSurroundingRomb(i_a, j_a); // is free other 3 spaces surrounding type K cell a within the romb
	countSurroundingRomb(i_b, j_b);  // is free other 3 spaces surrounding type K cell b within the romb
		
    if (!jumpFreeSpaces.empty()) {
        jumpCell.penalty = penalty(i_a, j_a) + penalty(i_b, j_b); // Calculate the penalty for the new cell
        ret = true;
    } 
	else {
        ret = false; // No free space to jump
    }
        
    return ret;
}

double randomGenerator(unsigned int first_interval, unsigned int last_interval) {
    std::random_device rd; // Create a random number generator engine
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(first_interval, last_interval);
    return dist(gen);
}

double penalty(int i, int j) {
	int count = countSurroundingOs(i, j); // true to calculate O type is presented
	double penalty = 0;
	
	if (i % 2 == 1 && j % 2 == 1) { // [odd][odd] cell of type K
		penalty = Delta[count];
		// grid[i][j].value = penalty;
	}
    return penalty;
}

bool metropolisCondition(double oldPenaltySum, double newPenaltySum) {
	
	bool ret = false;
	
	double randomValue = randomGenerator(0, 1); // Generate a random number in the range [0, 1]
		
	double metropolisValue = exp(-((newPenaltySum - oldPenaltySum) / kT_eV)); // Calculate the Metropolis condition value
	
    // Compare the random value with the Metropolis criteria
    if (randomValue > metropolisValue) {
        ret = false; // Reject the jump
    } else {
        ret =  true; // Accept the jump
    }
    return ret;
}
	
// finding free place to jump cell type O, if possible calculate penalty energy and fix current [i][j]
bool findOForJumping(tJumpCell& jumpCell) {
    bool ret = false;

    while (true) {
		int k = rand() % n;
        int m = rand() % n;
	
        // Cells type O [even][odd] or [odd][even]
        if (((k % 2 == 0) && (m % 2 == 1)) || ((k % 2 == 1) && (m % 2 == 0))) {

            if (grid[k][m].value == 8) {
				jumpCell.cellTypeO.i = k;
                jumpCell.cellTypeO.j = m;
                ret = rombPenaltyCalculation(jumpCell);
				
				if (ret == true) {
					break;
				}     
            }
        }
    }
    return ret;
}

bool initJumping(tJumpCell& jumpCell) {
	bool ret = false;
	
    if (true == findOForJumping(jumpCell)) {	
		ret = true;   
    }
	else {
		ret = false;
    }
    return ret;
}

void jumping(tJumpCell& jumpCell) {
    bool isJump = false;

	// Check if jumpFreeSpaces is not empty
	if (!jumpFreeSpaces.empty()) {
		
	    int randomIndex = static_cast<int>(randomGenerator(0, jumpFreeSpaces.size() - 1)); // Generate a random index within the range of jumpFreeSpaces
		
	    // Get the corresponding new cell to jump
	    tJumpCell jumpNewCell;
	    jumpNewCell.cellTypeO.i = jumpFreeSpaces[randomIndex][0].cellTypeO.i;
	    jumpNewCell.cellTypeO.j = jumpFreeSpaces[randomIndex][0].cellTypeO.j;
	    jumpNewCell.penalty = 0.0;

        rombPenaltyCalculation(jumpNewCell); // Calculate the penalty for the new cell

        if (jumpNewCell.penalty < jumpCell.penalty) {
            isJump = true;
        }
		else {
            isJump = metropolisCondition(jumpCell.penalty, jumpNewCell.penalty); // Implement metropolis condition and assign the result to isJump
        }

        if (isJump) {
			// Rewrite 8 to the new position in the grid vector
			grid[jumpNewCell.cellTypeO.i][jumpNewCell.cellTypeO.j].value = 8;
			grid[jumpCell.cellTypeO.i][jumpCell.cellTypeO.j].value = 1;
			jumpingCounts++;
        }

        // Clear fields of jumpCell
        jumpCell.cellTypeO.i = -1;
        jumpCell.cellTypeO.j = -1;
        jumpCell.penalty = 0.0;

        // Clear all elements of jumpFreeSpaces
        jumpFreeSpaces.clear();
    }
}

void evolution(void) {

    // Double iteration 1000000000; // Working
    double iteration = 100000; // Testing
    
    while (iteration--) {
    	
    	tJumpCell jumpCell; // current cell
    	
    	if (initJumping(jumpCell)) {
    		jumping(jumpCell); // checking local cell inside
		}
    }

    std::cout << "Evolution is done" << std::endl;
}

//***************************************************************************
bool printMatrixToTxt(const std::string& fileName) {
    std::ofstream outputFile(fileName);

    if (!outputFile.is_open()) {
        std::cerr << "Unable to open output file: " << fileName << std::endl;
        return false;
    }

    for (const std::vector<tRomb>& row : grid) {
        for (const tRomb& cell : row) {
            outputFile << cell.value << ' ';
        }
        outputFile << std::endl;
    }

    outputFile.close();
    return true;
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

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i % 2 == 1) && (j % 2 == 1)) { // [odd][odd] position calculation
                    unsigned int value = countSurroundingOs(i, j); // Number of oxygen atoms calculation surround [odd][odd] cell
                    // Use the color map to set the color based on the value
                    if (colorMap.find(value) != colorMap.end()) {
                        outputFile << colorMap[value];
                    }
                } else {
                    // Set white color for cells that do not meet the condition
                    outputFile << "255 255 255 ";
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

    std::cout << "Number of all elements = " << n*n << std::endl;

    grid.resize(n, std::vector<tRomb>(n, { {0, 0}, 1, 1 }));

    fillMatrix();

    int percentO; // Prompt the user to input the percentage of Oxygen (type O cell) in the matrix
    std::cout << "Input percentage of Oxygen: ";
    std::cin >> percentO;
	
	std::cout << "Number of Oxygen cells = " << (n * n / 2) << std::endl;
    numberO = (n * n / 2) * percentO / 100; // Calculate the number of O cells
    std::cout << "Number of Oxygen cells = " << numberO << std::endl;
    
    fillOxigen(); // Fill the grid with O cells
}

//************************************************************************************************
int main(int argc, char** argv) {

    configurator();
    
    // testing block
    //testOcupationO();
    // testing block

    if (!printMatrixToTxt("init.txt")) {
        std::cerr << "Failed to create 'init.txt'." << std::endl;
    }
    printMatrixToImage("init_map.ppm");    
    
    auto start = std::chrono::high_resolution_clock::now(); // Record the start time for measuring the duration of the fillOxigen function
    
	evolution();
	
    auto stop = std::chrono::high_resolution_clock::now(); // Record the stop time to calculate the duration of the fillOxigen function
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Time taken by fillOxigen = " << duration.count() << " microseconds" << std::endl;
    
    std::cout << "jumping counts = " << jumpingCounts << std::endl;
    
    // Print the grid matrix after evolution to a new file
    if (!printMatrixToTxt("evolution.txt")) {
        std::cerr << "Failed to create 'evolution.txt'." << std::endl;
    }
    printMatrixToImage("evolution_map.ppm");
		
    std::cin.get();

    return 0;
}



// testing block start
void testOcupationO(void) {
	std::cout << "testing start" << std::endl;
    for (int k = 0; k < 5; ++k) {
		int i = static_cast<int>(randomGenerator(0, n)) % n;
        int j = static_cast<int>(randomGenerator(0, n)) % n;

        std::cout << k << ": [" << i << "][" << j << "] = " << countSurroundingOs(i, j) << " cells" << std::endl;
    }
    std::cout << "testing start" << std::endl;
}
// testing block stop

