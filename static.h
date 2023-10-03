#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <chrono>


int n = 10; // Розмір масиву 10x10
int numberO;
std::vector<std::vector<int>> grid;

int countSurroundingOs(int i, int j);
void fillOxigen();


int countSurroundingOs(int i, int j) {
    int count = 0;

    // Перевірка, чи комірка "К" парна (парна, парна) - 0,1,2,3,4... парні{1,3,5}
    if (i % 2 == 1 && j % 2 == 1) {
        // Перевірка 4 можливих напрямків навколо комірки "К" [-i][j]; [+i][j]; [i][-j]; [i][+j]
        int d[] = {-1, 1};
 
        for (int dir = 0; dir < 2; dir++) {
            int ni = (i + d[dir] + n) % n;
            int nj = (j + d[dir] + n) % n;
			
            if (grid[ni][nj] == -1) {
                count++;
            }
        }
    }

    return count;
}

void fillOxigen() {
    // Заповнення комірок типу "О" випадково
    while (numberO > 0) {
        int i = rand() % n;
        int j = rand() % n;
	
        if (grid[i][j] == 0) {
            grid[i][j] = -1; // Заповнення комірки типу "О"
            numberO--;
        }
    }
}



int main() {

    grid.resize(n, std::vector<int>(n, 0)); // Ініціалізація масиву по комірках

    // Ввести відсоток комірок типу "О"
    int percentO;
    std::cout << "Input percentage of Oxygen: ";
    std::cin >> percentO;

    numberO = (n * n / 2) * percentO / 100;


	auto start = std::chrono::high_resolution_clock::now(); // Початок вимірювання часу
	
    fillOxigen();

    auto stop = std::chrono::high_resolution_clock::now(); // Кінець вимірювання часу
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Time taken by fillOxigen: " << duration.count() << " microseconds" << std::endl;

    int i = rand() % n;
    int j = rand() % n;

    std::cout << "i = " << i << std::endl;
    std::cout << "j = " << j << std::endl;
    std::cout << "Number of Oxygen atoms: " << countSurroundingOs(i, j) << std::endl;

    // Вивести масив у файл
    std::ofstream outputFile("grid.txt");
    
    std::cin.get();
    
    return 0;
}
