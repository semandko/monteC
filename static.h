#ifndef _OXIGEN_ARRAY_STATIC
#define _OXIGEN_ARRAY_STATIC

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <chrono>


int n = 500;
int numberO;
std::vector<std::vector<int>> grid;

void fillMatrix();
int countSurroundingOs(int i, int j);
void fillOxigen();

void fillMatrix() {	
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
			if ((i % 2 == 1) && (j % 2 == 1)) {
    			grid[i][j] = 0.0;
			}
        }
    } 
}

int countSurroundingOs(int i, int j) {
    int count = 0;

    if (i % 2 == 1 && j % 2 == 1) {
        int d[] = {-1, 1};
 
        for (int dir = 0; dir < 2; dir++) {
            int ni = (i + d[dir] + n) % n;
            int nj = (j + d[dir] + n) % n;
			
            if (grid[ni][nj] == -1) {
                count++;
            }
        }
    }
    else {
    	std::cout << "not [even][even] cell" << std::endl;
	}

    return count;
}

void fillOxigen() {
    while (numberO > 0) {
        int i = rand() % n;
        int j = rand() % n;
	
		if (((i % 2 == 0) && (j % 2 == 1)) || ((i % 2 == 1) && (j % 2 == 0))) {
	        if (grid[i][j] != -1) {
	            grid[i][j] = -1;
	            numberO--;
	        }
    	}
    }
}

#endif

