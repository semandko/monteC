#include <iostream>
#include <thread>
#include <mutex>
#include <vector>

const int num_philosophers = 5;
std::vector<std::mutex> forks(num_philosophers);

void philosopher(int id) {
    int left_fork = id;
    int right_fork = (id + 1) % num_philosophers; //0..4; example l = 4; r = 5%5 = 0

    for (auto i = 0; i < 10; ++i) {
        // Think
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        std::cout << "Philosopher " << id << " is thinking. " << std::endl;

        // Pick up forks
        forks[left_fork].lock();  // use mutex 1
        forks[right_fork].lock(); // use mutex 2

        // Eat
        std::cout << "Philosopher " << id << " is eating. " << std::endl;

        // Put down forks
        forks[left_fork].unlock();  // leave mutex 1
        forks[right_fork].unlock(); // leave mutex 2
    }
}

int main() {
    std::vector<std::thread> philosophers;

    for (int i = 0; i < num_philosophers; ++i) {
        philosophers.emplace_back(philosopher, i);
    }

    for (auto& philosopher_thread : philosophers) {
        philosopher_thread.join();
    }

	system("pause");
    return 0;
}
