#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

int main() {
    std::ofstream outfile("random_data.txt"); // 生成的文件名为 random_data.txt
    if (!outfile) {
        std::cerr << "Error: Could not open file for writing!" << std::endl;
        return 1;
    }

    std::srand(std::time(0));  // 设置随机数种子

    for (int i = 0; i < 16384; ++i) {
        unsigned long long random_value = static_cast<unsigned long long>(std::rand()) << 32 | std::rand(); // 生成一个随机的64位整数
        outfile << random_value << std::endl;
    }

    std::cout << "File random_data.txt has been generated with 16384 random 64-bit integers!" << std::endl;
    return 0;
}
