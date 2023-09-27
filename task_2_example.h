#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include <mutex>
#include <chrono>


const double p = 3.14159265359; // Pi
const int num_threads = 8; 		// Adjust the number of threads as needed
const int num_points = 1000000; // Adjust the number of points as needed
double total_area = 0.0; 		// common result
 
std::mutex result_mutex; // Mutex to protect the shared result

// Function to calculate the area under the curve for a given range
double calculate_area(double start, double end) {
    double h = (end - start) / num_points;
    double sum = 0.0;

    for (int i = 0; i < num_points; ++i) {
        double x0 = start + i * h;
        double x1 = x0 + h;
        double y0 = std::sin(x0);
        double y1 = std::sin(x1);
        sum += (y0 + y1) * h / 2.0;
    }

    return sum;
}

// Function to compute the area in a specific range and accumulate it
void compute_range(double start, double end) {
	
	auto start_time = std::chrono::high_resolution_clock::now(); // Start time

    double area = calculate_area(start, end);

    auto end_time = std::chrono::high_resolution_clock::now(); // End time
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	
	std::cout << "Execution Time: " << duration.count() << " microseconds" << std::endl;
	
    // Protect the shared result of total_area with a mutex
    std::lock_guard<std::mutex> lock(result_mutex);
    total_area += area;
}

//
int main() {
	
    std::vector<std::thread> threads;

    // Split the range into equal parts for parallel processing
    double step = p / num_threads;
    
    for (int i = 0; i < num_threads; ++i) {
        double start = i * step;
        double end = start + step;
        threads.emplace_back(compute_range, start, end);
    }

    // Wait for all threads to finish
    for (std::thread& thread : threads) {
        thread.join();
    }

    std::cout << "Total area under the curve: " << total_area << std::endl;

    return 0;
}
