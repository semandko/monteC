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
#define _SILICON(i,j) (i % 2 == 1) && (j % 2 == 1) // cels type Si [odd][odd]
#define _OXIGEN(i,j) ((i % 2 == 0) && (j % 2 == 1)) || ((i % 2 == 1) && (j % 2 == 0)) // cells type Oxygen [even][odd] or [odd][even]


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
	float penalty;
} tJumpCell;


// variables
int n = 500;
uint32_t iteration = 1000000000; // Testing

float x; // Prompt the user to input the index of Oxygen (type O cell) in the matrix
int numberO;

std::vector<std::vector<tRomb>> grid(n, std::vector<tRomb>(n, { {0, 0}, 1, 1 }));
std::vector<tJumpCell> jumpFreeSpaces; // 

std::vector<float> Delta{0.0, 0.5, 0.51, 0.22, 0.0}; // penalty energy

// surrounding for romb
int di[] = {-1, 1, 0, 0}; // vertical
int dj[] = {0, 0, -1, 1}; // horizontal
		
tJumpCell jumpNewCell; // cell in jumping
tJumpCell jumpCell; // current cell

float penaltyValue;

float x0 = _X_0;
float kB = _K_B;
float kT_eV = _kT_EV(_T);
int T = _T;
int jumpingCounts;

// prototypes
void fillMatrix();
void fillOxigen();
int countSurroundingOs(int i, int j, bool isPresent);
void countSurroundingRomb(int i, int j);
float randomGenerator(unsigned int first_interval, unsigned int last_interval);
float penalty(int i, int j);
bool metropolisCondition(float oldPenaltySum, float newPenaltySum);
bool rombPenaltyCalculation(tJumpCell& jumpCell);
bool findOForJumping(tJumpCell& jumpCell);
void jumping(tJumpCell& jumpCell);
void evolution(void);

bool printMatrixToTxt(void);
void printMatrixToImage(const std::string& fileName);

void configurator(void);

// testing block
void countCells(void);
// testing block


// code section *************************************************

// not optimized because init
void fillMatrix() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (_SILICON(i,j)) {
                grid[i][j].value = 0; // Set value to 1 for [odd][odd] cells
            }
        }
    }
}

// not optimized
void fillOxigen() {
	int nOxigen = numberO;
    while (nOxigen > 0) {
        int i = rand() % n;
        int j = rand() % n;
		
        if (_OXIGEN(i,j)) {
            if (grid[i][j].value != 2) {
                grid[i][j].value = 2;
                nOxigen--;
            }
        }
    }
}

float randomGenerator(unsigned int first_interval, unsigned int last_interval) {
    int randomValue = rand();  // Generate a random integer using rand()
    
    float scaledValue = first_interval + (randomValue / (RAND_MAX + 1.0)) * (last_interval - first_interval);

    return scaledValue;
}


// may be optimized
int countSurroundingOs(int i, int j) {
	int countPresent = 0;

	for (int dir = 0; dir < 4; dir++) {
		int ni = (i + di[dir] + n) % n; // type O index
		int nj = (j + dj[dir] + n) % n; // type O index

		// is type O presented
		if (grid[ni][nj].value == 2) {
			countPresent++;
		}
	}
    return countPresent;
}

// may be optimized
void countSurroundingRomb(int i, int j) {
    if (_SILICON(i,j)) {

		tJumpCell jumpCell;
        for (int dir = 0; dir < 4; dir++) {
            int ni = (i + di[dir] + n) % n; // type O index
            int nj = (j + dj[dir] + n) % n; // type O index

            // Check if it's a free space
            if (grid[ni][nj].value != 2) {
            	
            	// structure tJumpCell     
                jumpCell.cellTypeO.i = ni;
                jumpCell.cellTypeO.j = nj;
                jumpCell.penalty = 0.0;
                
                jumpFreeSpaces.push_back(jumpCell); // Push the vector containing jumpCell
            }
        }
    } 
	else
	{
        std::cout << "not valid " << std::endl;
    }
}

// may be optimized
bool rombPenaltyCalculation(tJumpCell& jumpCell) {
    bool ret = false;
    int i = jumpCell.cellTypeO.i; // for type O
    int j = jumpCell.cellTypeO.j; // for type O
	int i_a, j_a;
	int i_b, j_b;
	jumpCell.penalty = 0;
	
    if (i % 2 == 1) { // Index for type O
        // Horizontal case for cells type K
		i_a = i; j_a = (j - 1 + n) % n; // left
		i_b = i; j_b = (j + 1) % n; // right
    }
	else
	{
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
	else
	{
        ret = false; // No free space to jump
    }
        
    return ret;
}

// may be optimized
float penalty(int i, int j) {
	float penalty = 0;
	
	if (_SILICON(i,j)) {
		int count = countSurroundingOs(i, j);
		penalty = Delta[count];
	}
    return penalty;
}

// may be optimized
bool metropolisCondition(float oldPenaltySum, float newPenaltySum) {
	
	bool ret = false;
	
	float randomValue = randomGenerator(0, 1); // Generate a random number in the range [0, 1]
		
	float metropolisValue = exp(-((newPenaltySum - oldPenaltySum) / kT_eV)); // Calculate the Metropolis condition value
	
    // Compare the random value with the Metropolis criteria
    if (randomValue > metropolisValue) {
        ret = false; // Reject the jump
    }
	else
	{
        ret =  true; // Accept the jump
    }
    return ret;
}

// may be optimized	
// finding free place to jump cell type O, if possible calculate penalty energy and fix current [i][j]
bool findOForJumping(tJumpCell& jumpCell) {
    bool ret = false;

    while (true) {
		int i = rand() % n;
        int j = rand() % n;
	
        if (_OXIGEN(i,j)) {
            if (grid[i][j].value == 2) {
				jumpCell.cellTypeO.i = i;
                jumpCell.cellTypeO.j = j;
                ret = rombPenaltyCalculation(jumpCell);
				
				if (ret == true) {
					break;
				}     
            }
        }
    }
    return ret;
}

// may be optimized
void jumping(tJumpCell& jumpCell) {
    bool isJump = false;

	// Check if jumpFreeSpaces is not empty
	if (!jumpFreeSpaces.empty()) {
		
	    int randomIndex = static_cast<int>(randomGenerator(0, jumpFreeSpaces.size() - 1)); // Generate a random index within the range of jumpFreeSpaces
		
	    // Get the corresponding new cell to jump
	    jumpNewCell.cellTypeO.i = jumpFreeSpaces[randomIndex].cellTypeO.i;
	    jumpNewCell.cellTypeO.j = jumpFreeSpaces[randomIndex].cellTypeO.j;
	    jumpNewCell.penalty = 0.0;

        rombPenaltyCalculation(jumpNewCell); // Calculate the penalty for the new cell

        if (jumpNewCell.penalty < jumpCell.penalty) {
            isJump = true;
        }
		else {
            isJump = metropolisCondition(jumpCell.penalty, jumpNewCell.penalty); // Implement metropolis condition and assign the result to isJump
        }

        if (isJump == true) {
			// Rewrite 2 to the new position in the grid vector
			grid[jumpNewCell.cellTypeO.i][jumpNewCell.cellTypeO.j].value = 2;
			grid[jumpCell.cellTypeO.i][jumpCell.cellTypeO.j].value = 7; // clear old place
			jumpingCounts++;
        }

        jumpFreeSpaces.clear(); // Clear all elements of jumpFreeSpaces
    }
}

// main loop
void evolution(void) {
    
    while (iteration--) {
    	if (findOForJumping(jumpCell) == true) {
    		jumping(jumpCell); // checking local cell inside
		}
    }

    std::cout << "Evolution is done" << std::endl;
}

//***************************************************************************
// output function
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

// output function
void printMatrixToImage(const std::string& fileName) {
    std::ofstream outputFile(fileName);
    
    if (outputFile.is_open()) {
        outputFile << "P3\n" << n << " " << n << "\n255\n";
        std::map<int, std::string> colorMap;
        colorMap[0] = "0 0 0 ";        // Black
        colorMap[1] = "255 0 0 ";      // Red
        colorMap[2] = "0 0 255 ";      // Blue
        colorMap[3] = "0 255 0 ";      // Green
        colorMap[4] = "255 255 0 ";    // Yellow
		
		unsigned int value =0;
		
		int countSi = 0;
		int countO = 0;
		int countGrid = 0;
		
		for (int i = 0; i < n; i++) {
		    for (int j = 0; j < n; j++) {
		    	
		        if (_SILICON(i, j)) { // [odd][odd] position calculation
		            value = countSurroundingOs(i, j);
		            
		            if (colorMap.find(value) != colorMap.end()) {
		                outputFile << colorMap[value];
		                countSi++;
		            }
		        }
		        else
				{
		            if (_OXIGEN(i, j)) {
		                outputFile << "255 192 203 "; // pink for O
		                countO++;
		            }
		            else
					{
		                // [even][even] position calculation
		                outputFile << "255 255 255 "; // white
		                countGrid++;
		            }
		        }
		    }
		    outputFile << "\n";
		}
        
        outputFile.close();
        std::cout << "Matrix has been saved " << fileName << std::endl;
        
        std::cout << "Count Si " << countSi << std::endl;
        std::cout << "Count Ox " << countO << std::endl;
        std::cout << "Count Grid " << countGrid << std::endl;
    }
	else
	{
        std::cerr << "Unable to save the image " << fileName << std::endl;
    }
}

// control output function
void countCells(void) {
	int countSi = 0;
	int countO = 0;
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (_SILICON(i,j)) {
				if (grid[i][j].value == 0) {
					countSi++;
				}
			}
			else {
				if (_OXIGEN(i,j)) {
					if (grid[i][j].value == 2) {
						countO++;
					}
				}
				else {
					if (grid[i][j].value != 1) {
						std::cout << "Not permitted cell = " << grid[i][j].value << std::endl;
					}
				}
			}
		}
	}
	std::cout << "Silicon stop = " << countSi << "; Silicon start = " << n*n / 4 << std::endl;
	std::cout << "Oxigen stop = " << countO << "; Oxigen start = " << numberO << std::endl;
}

// init function
void configurator(void) {

    std::cout << "Elements = " << n*n << std::endl;

    grid.resize(n, std::vector<tRomb>(n, { {0, 0}, 1, 1 }));

    fillMatrix();
	std::cout << "Silicon cells = " << (n * n / 4) << std::endl;
	std::cout << "Oxygen cells = " << (n * n / 2) << std::endl;
	
	while (1) {
		std::cout << "Input stoicheometry index x of Oxygen [0,0.1..2]: ";
		std::cin >> x;
		if ( x < 0 || x > 2) {
			std::cout << "NOTE! x of Oxygen must be in the range [0,0.1..2], you put x = " << x << std::endl;
		}
		else {
			numberO = x*((n * n) / 4); // Calculate the number of O cells depending on stoicheometry SiOx
			std::cout << "Number of Oxygen = " << numberO << std::endl;
			break;
		}
	}

    fillOxigen(); // Fill the grid with O cells
}

//************************************************************************************************
int main(int argc, char** argv) {
	
	int iter = iteration;
    configurator();

    // Count Si; Oxygen

    if (!printMatrixToTxt("init.txt")) {
        std::cerr << "Failed to create 'init.txt'." << std::endl;
    }
    printMatrixToImage("init_map.ppm");

	// estimation
    auto start = std::chrono::high_resolution_clock::now();
    evolution();
    auto stop = std::chrono::high_resolution_clock::now();
    
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    int hours = std::chrono::duration_cast<std::chrono::hours>(duration).count();
    duration -= std::chrono::hours(hours);
    int minutes = std::chrono::duration_cast<std::chrono::minutes>(duration).count();
    duration -= std::chrono::minutes(minutes);
    int seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
    duration -= std::chrono::seconds(seconds);
    int microseconds = duration.count();

    std::cout << "Time taken by fillOxigen = " << hours << " hours, " << minutes << " minutes, " << seconds << " seconds, and " << microseconds << " microseconds" << std::endl;

    // Print the grid matrix after evolution to a new file
    if (!printMatrixToTxt("evolution.txt")) {
        std::cerr << "Failed to create 'evolution.txt'." << std::endl;
    }
    printMatrixToImage("evolution_map.ppm");

	std::cout << "iteration = " << iter << std::endl;
    std::cout << "jumping counts = " << jumpingCounts << std::endl;

    //countCells();

    std::cin.get();

    return 0;
}

