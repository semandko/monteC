#include <iostream>
#include <vector>
#include <fstream>
#include <mutex>
#include <thread>
#include <random>

int n = 50; // Розмір масиву
int numberO;
std::vector<std::vector<int>> grid; // оголосіть глобально або передавайте як параметр
const int num_threads = 8; // Виправте кількість потоків за необхідності
std::mutex result_mutex;

int countSurroundingOs(int i, int j) {
    int count = 0;

    // Перевірка, чи комірка "К" парна-непарна
    if (((i % 2 == 0) && (j % 2 == 1)) || ((i % 2 == 1) && (j % 2 == 0))) {
        // Перевірка 4 можливих напрямків навколо комірки "К"
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

void fillOxygen() {
    // Генерація комірок типу "О" випадково
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, n - 1);

    while (numberO > 0) {
        int i = dis(gen);
        int j = dis(gen);

		if (((i % 2 == 0) && (j % 2 == 1)) || ((i % 2 == 1) && (j % 2 == 0))) {
			
			// Захист результату від декількох потоків за допомогою mutex
			std::lock_guard<std::mutex> lock(result_mutex);

			if (grid[i][j] == 0) {
				grid[i][j] = -1; // Заповнення комірки типу "О"

				if (numberO > 0) {
					numberO--;
				} else {
					break;
				}
			}
		}
    }

    // Зупинка всіх потоків, якщо numberO == 0
    if (numberO == 0) {
        std::terminate();
    }
}

int main() {
    std::thread threads[num_threads];
    grid.resize(n, std::vector<int>(n, 0)); // Ініціалізація масиву по комірках

    // Ввести відсоток комірок типу "О"
    int percentO;
    std::cout << "Input percentage of Oxygen: ";
    std::cin >> percentO;

    numberO = (n * n / 2) * percentO / 100;

    auto start = std::chrono::high_resolution_clock::now(); // Початок вимірювання часу

    for (int i = 0; i < num_threads; ++i) {
        threads[i] = std::thread(fillOxygen);
    }

    // Wait for all threads to finish
    for (int i = 0; i < num_threads; ++i) {
        threads[i].join();
    }

    auto stop = std::chrono::high_resolution_clock::now(); // Кінець вимірювання часу
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Time taken by fillOxygen: " << duration.count() << " microseconds" << std::endl;

    int i = std::rand() % n;
    int j = std::rand() % n;

    std::cout << "i = " << i << std::endl;
    std::cout << "j = " << j << std::endl;
    std::cout << "Number of Oxygen atoms: " << countSurroundingOs(i, j) << std::endl;

    // Вивести масив у файл
    std::ofstream outputFile("grid.txt");

    std::cin.get();

    return 0;
}
