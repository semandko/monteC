#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include <mutex>
#include <queue>

const double p = 3.14159265359; // Pi
const int mainThread = 1;
const int numThreads = 8 - mainThread; // Adjust the number of threads as needed
const int numPoints = 1000000; // Adjust the number of points as needed
double totalArea = 0.0; // Common result
double step = p / (numThreads); // Split the range into equal parts for parallel processing

std::mutex resultMutex; // Mutex to protect the shared result

typedef struct TSBlock {
    double start;
    double end;
} tsBlock;

std::queue<tsBlock> areaBlockQueue; // Queue to store sub-ranges

// Function to calculate the area under the curve for a given range
double calculateArea(double start, double end) {
    double h = (end - start) / numPoints;
    double sum = 0.0;

    for (int i = 0; i < numPoints; ++i) {
        double x0 = start + i * h;
        double x1 = x0 + h;
        double y0 = std::sin(x0);
        double y1 = std::sin(x1);
        sum += (y0 + y1) * h / 2.0;
    }

    return sum;
}

// Function to compute the area in a specific range and accumulate it
void computeRange(double start, double end) {
    double area = calculateArea(start, end);

    // Protect the shared result of total_area with a mutex
    std::lock_guard<std::mutex> lock(resultMutex);
    totalArea += area;
}

// Function to compute the sub-ranges and add them to the queue
void computeQueue() {
    for (int i = 0; i < numThreads; ++i) {
        tsBlock areaBlock;
        areaBlock.start = i * step;
        areaBlock.end = areaBlock.start + step;

        // Push the sub-range to the back of the queue
        std::lock_guard<std::mutex> lock(resultMutex);
        areaBlockQueue.push(areaBlock);
    }
}

// Worker function for each thread
void worker() {
    while (true) {
        tsBlock areaBlock;

        {
            // Lock the queue and pop an item from the front
            std::lock_guard<std::mutex> lock(resultMutex);
            if (areaBlockQueue.empty()) {
                break; // Exit the loop if the queue is empty
            }
            areaBlock = areaBlockQueue.front();
            areaBlockQueue.pop();
        }

        // Compute the area for the sub-range
        computeRange(areaBlock.start, areaBlock.end);
    }
}

int main() {
    std::thread mainThread;
    std::vector<std::thread> threads;

    mainThread = std::thread(computeQueue);

    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back(worker);
    }

    // Wait for the main thread and worker threads to finish
    mainThread.join();
    for (std::thread& thread : threads) {
        thread.join();
    }

    std::cout << "Total area under the curve: " << totalArea << std::endl;

    return 0;
}
