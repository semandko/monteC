#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <map>
#include <cmath>
#include <random>


#define _TESTING 1

#define _CONST_A 1.38e-23
#define _CONST_B 1.602e-19
#define _K_B (_CONST_A / _CONST_B)
#define _X_0 0.1
#define _T (1000)
#define _kT_EV(_T) (1.0 / 11604 * (_T + 273.15)) // 0.11 eV for 1000 C

#define _IS_SILICON 0
#define _IS_GRID_CORE 1
#define _IS_JUMPED 7
#define _IS_OXIGEN_PRESENTED 2
//
#define _SILICON(i,j) 			(i % 2 == 1) && (j % 2 == 1) // cels type Si [odd][odd]
#define _OXIGEN(i,j) 			((i % 2 == 0) && (j % 2 == 1)) || ((i % 2 == 1) && (j % 2 == 0)) // cells type Oxygen [even][odd] or [odd][even]

#define _IS_HORIZONTAL 			(i % 2 == 1)
#define _SI_HORIZONTAL_LEFT 	((j - 1 + n) % (n))
#define _SI_HORIZONTAL_RIGHT 	((j + 1) % (n))
#define _SI_VERTICAL_UP 		((i - 1 + n) % (n))
#define _SI_VERTICAL_DOWN 		((i + 1) % (n))


// typedefs
typedef struct {
    int i;
    int j;
} tcell;

typedef struct {
	tcell cellTypeO;
	float penalty;
} tJumpCell;


// variables
int n = 500;
uint32_t iteration = 1000000000; // Testing

float x; // Prompt the user to input the index of Oxygen (type O cell) in the matrix
int numberO; // Oxigen counter

int** grid; // matrix of elements

std::vector<tJumpCell> jumpFreeSpaces; // 

float Delta[] = {0.0, 0.5, 0.51, 0.22, 0.0}; // penalty energy

// surrounding for romb
int di[] = {-1, 1, 0, 0}; // vertical
int dj[] = {0, 0, -1, 1}; // horizontal

tJumpCell jumpOldCell; // current cell	
tJumpCell jumpNewCell; // cell in jumping
tJumpCell freeCell; // free cell to jump
	
float penaltyValue;

float x0 = _X_0;
float kB = _K_B;
float kT_eV = _kT_EV(_T);
int T = _T;

#ifdef _TESTING
int jumpingCounts = 0;
int jumpingRejected = 0;
int jumpingMetropolis = 0;
int jumpingNoMetropolis = 0;
int wrongSiIndex = 0;
#endif


// prototypes
void fillMatrix();
void fillOxigen();
int isPresentedOxigen(int i, int j, bool isPresent);
void isNotPresentedOxigen(int i, int j);
float randomGenerator(unsigned int first_interval, unsigned int last_interval);
float penalty(int i, int j);
bool metropolisCondition(float oldPenaltySum, float newPenaltySum);
bool rombPenaltyCalculation(tJumpCell& jumpCell);
bool findOxigen(void);
void jumping(void);
void evolution(void);

bool printMatrixToTxt(void);
void printMatrixToImage(const std::string& fileName);

void configurator(void);

#ifdef _TESTING
void countCells(void);
void distributionSi(void);
#endif


// code section *************************************************

// not optimized because init
void fillMatrix() {
	grid = new int*[n];
    for (int i = 0; i < n; ++i) {
        grid[i] = new int[n];
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (_SILICON(i,j)) {
                grid[i][j] = _IS_SILICON; // Set value to 0 for [odd][odd] cells
            }
            else {
            	grid[i][j] = _IS_GRID_CORE;
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
            if (grid[i][j] != _IS_OXIGEN_PRESENTED) {
                grid[i][j] = _IS_OXIGEN_PRESENTED;
                nOxigen--;
            }
        }
    }
}

// generate index for new position within the presented free localized cells
float randomGenerator(unsigned int first_interval, unsigned int last_interval) {
    int randomValue = rand(); 
    float scaledValue = first_interval + (randomValue / (RAND_MAX + 1.0)) * (last_interval - first_interval);

    return scaledValue;
}

// int i, int j for Silicon
int isPresentedOxigen(int i, int j) {
	int countPresent = 0;
	int ni = 0; // ni Oxigen index
	int nj = 0; // nj Oxigen index
	
	for (int dir = 0; dir < 4; dir++) {
		ni = (i + di[dir] + n) % n;
		nj = (j + dj[dir] + n) % n;

		if (grid[ni][nj] == _IS_OXIGEN_PRESENTED) {
			countPresent++;
		}
	}
    return countPresent;
}

// int i, int j for Silicon
void isNotPresentedOxigen(int i, int j) {
	int ni = 0; // ni Oxigen index
	int nj = 0; // nj Oxigen index
	
    for (int dir = 0; dir < 4; dir++) {
        ni = (i + di[dir] + n) % n;
        nj = (j + dj[dir] + n) % n;

        // Check if it's a free space
        if (grid[ni][nj] != _IS_OXIGEN_PRESENTED) {  
            freeCell.cellTypeO.i = ni;
            freeCell.cellTypeO.j = nj;
            freeCell.penalty = 0.0;
            
            jumpFreeSpaces.push_back(freeCell); // collect free spaces
        }
    }
}

// may be optimized
float penalty(int i, int j) {
	float penalty = 0;
	
	if (_SILICON(i,j)) {
		int count = isPresentedOxigen(i, j);
		penalty = Delta[count];
	}
	#ifdef _TESTING
	else {
		wrongSiIndex++;
	}
	#endif
    return penalty;
}

// may be optimized
bool rombPenaltyCalculation(tJumpCell& jumpCell) {
    bool ret = false;
    
    int i = jumpCell.cellTypeO.i; // for Oxigen
    int j = jumpCell.cellTypeO.j; // for Oxigen
	int i_a, j_a; // for Silicon
	int i_b, j_b; // for Silicon
	jumpCell.penalty = 0;
	
    if (_IS_HORIZONTAL) {
        // Horizontal case for Silicon = cell a + cell b
		i_a = i; 					j_a = _SI_HORIZONTAL_LEFT;
		i_b = i; 					j_b = _SI_HORIZONTAL_RIGHT;
    }
	else {
        // Vertical case for Silicon  = cell a + cell b
        i_a = _SI_VERTICAL_UP; 	 	j_a = j;
		i_b = _SI_VERTICAL_DOWN; 	j_b = j;
    }
    
    isNotPresentedOxigen(i_a, j_a); // is free other 3 spaces surrounding type Silicon cell a within the romb
	isNotPresentedOxigen(i_b, j_b);  // is free other 3 spaces surrounding type Silicon cell b within the romb
	
    if (!jumpFreeSpaces.empty()) {
        jumpCell.penalty = penalty(i_a, j_a) + penalty(i_b, j_b); // Calculate the penalty for the new cell
        ret = true;
    } 
	else {
        ret = false; // No free space to jump
    }
        
    return ret;
}

// may be optimized
bool metropolisCondition(float oldPenaltySum, float newPenaltySum) {
	
	bool ret = false;
	
	float randomValue = randomGenerator(0, 1); // Generate a random number in the range [0, 1]
		
	float metropolisValue = exp(-((newPenaltySum - oldPenaltySum) / kT_eV)); // Calculate the Metropolis condition value
	
    // Compare the random value with the Metropolis criteria
    if (randomValue > metropolisValue) {
        ret = false; // Reject the jump
        #ifdef _TESTING
        jumpingNoMetropolis++;
        #endif
    }
	else {
        ret =  true; // Accept the jump
        #ifdef _TESTING
        jumpingMetropolis++;
        #endif
    }
    return ret;
}

// may be optimized	
// finding free place to jump cell with oxigen, if possible calculate penalty energy and fix current position [i][j]
bool findOxigen() {
    bool ret = false;
	int i = 0;
	int j = 0;
	
    while (true) {
		i = rand() % n;
        j = rand() % n;
	
        if (_OXIGEN(i,j)) {
            if (grid[i][j] == _IS_OXIGEN_PRESENTED) {
				jumpOldCell.cellTypeO.i = i;
                jumpOldCell.cellTypeO.j = j;
                ret = rombPenaltyCalculation(jumpOldCell);
				
				break;  
            }
        }
    }
    // true = is possible
    // false = not possible (no free spaces)
    return ret;
}

// may be optimized
void jumping() {
    bool isJump = false;
	
	// generate new cell to jump	
    int randomIndex = static_cast<int>(randomGenerator(0, jumpFreeSpaces.size() - 1)); // Generate a random index within the range of jumpFreeSpaces
	
    // Get the corresponding new cell to jump
    jumpNewCell.cellTypeO.i = jumpFreeSpaces[randomIndex].cellTypeO.i;
    jumpNewCell.cellTypeO.j = jumpFreeSpaces[randomIndex].cellTypeO.j;
    jumpNewCell.penalty = 0.0;

    rombPenaltyCalculation(jumpNewCell); // Calculate the penalty for the new cell

    if (jumpNewCell.penalty < jumpOldCell.penalty) {
        isJump = true;
        #ifdef _TESTING
        jumpingCounts++;
        #endif
    }
	else {
        isJump = metropolisCondition(jumpOldCell.penalty, jumpNewCell.penalty); // Implement metropolis condition and assign the result to isJump
    }

    if (isJump == true) {
		// Rewrite 2 to the new position in the grid vector
		grid[jumpNewCell.cellTypeO.i][jumpNewCell.cellTypeO.j] = _IS_OXIGEN_PRESENTED;
		grid[jumpOldCell.cellTypeO.i][jumpOldCell.cellTypeO.j] = _IS_JUMPED; // clear old place		
    }

    jumpFreeSpaces.clear(); // Clear all elements of jumpFreeSpaces
}

// main loop
void evolution(void) {
    
    while (iteration--) {
    	if (findOxigen() == true) {
    		jumping(); // checking local cell inside
		}
		#ifdef _TESTING
		else {
			jumpingRejected++;
		}
		#endif
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

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            outputFile << grid[i][j] << ' ';
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
		
		unsigned int value = 0;
		
		int countSi = 0;
		int countO = 0;
		int countGrid = 0;
		
		for (int i = 0; i < n; i++) {
		    for (int j = 0; j < n; j++) {
		    	
		        if (_SILICON(i, j)) { // [odd][odd] position calculation
		            value = isPresentedOxigen(i, j);
		            
		            if (colorMap.find(value) != colorMap.end()) {
		                outputFile << colorMap[value];
		                countSi++;
		            }
		        }
		        else {
		            if (_OXIGEN(i, j)) {
		                outputFile << colorMap[0];
		                countO++;
		            }
		            else {
		                // [even][even] position calculation - not permitted position
		                outputFile << colorMap[0];
		                countGrid++;
		            }
		        }
		    }
		    outputFile << "\n";
		}
        
        outputFile.close();
        std::cout << "Matrix has been saved " << fileName << std::endl;
        
        #if 1
        std::cout << "Count Si " << countSi << std::endl;
        std::cout << "Count Ox " << countO << std::endl;
        std::cout << "Count Grid " << countGrid << std::endl;
        #endif
    }
	else {
        std::cerr << "Unable to save the image " << fileName << std::endl;
    }
}

#ifdef _TESTING
// testing control output function
void countCells(void) {
	int countSi = 0;
	int countO = 0;
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (_SILICON(i,j)) {
				if (grid[i][j] == _IS_SILICON) {
					countSi++;
				}
			}
			else {
				if (_OXIGEN(i,j)) {
					if (grid[i][j] == _IS_OXIGEN_PRESENTED) {
						countO++;
					}
				}
				else {
					if (grid[i][j] != _IS_GRID_CORE) {
						std::cout << "Not permitted cell = " << grid[i][j] << std::endl;
					}
				}
			}
		}
	}
	std::cout << "Silicon stop = " << countSi << "; Silicon start = " << n*n / 4 << std::endl;
	std::cout << "Oxigen stop = " << countO << "; Oxigen start = " << numberO << std::endl;
	
	std::cout << "fast jumping = " << jumpingCounts << std::endl;
	std::cout << "no jumping full cells = " << jumpingRejected << std::endl;
	std::cout << "metrop jumping = " << jumpingMetropolis << std::endl;
	std::cout << "metrop rejected = " << jumpingNoMetropolis << std::endl;
	std::cout << "wrong index for Si = " << wrongSiIndex << std::endl;
}

void distributionSi(void) {
    unsigned int value[5] = {0, 0, 0, 0, 0};

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
        	if (_SILICON(i,j)) {
            	value[isPresentedOxigen(i, j) % 5]++;
        	}
        }
    }
    
    for (int i = 0; i < 5; ++i) {
    	std::cout << "Si[ " << i << " ] = " << value[i] << std::endl;
	}
}

#endif



// init function
void configurator(void) {

    std::cout << "Elements = " << n*n << std::endl;

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

#ifdef _TESTING
	distributionSi();
#endif
	
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

#ifdef _TESTING
    countCells();
    distributionSi();
#endif

    // Don't forget to free the memory when you're done
    for (int i = 0; i < n; ++i) {
        delete[] grid[i];
    }
    
	delete[] grid;
    

    std::cin.get();

    return 0;
}

