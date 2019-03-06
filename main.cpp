#include <iostream>
#include <string>

#include "matrices.hpp"

int main() {
    try {
        matrix m1, m2;

        std::cout << "Input first matrix (rows and columns on the first line, matrix on the following lines, all separated by spaces)" << std::endl;

        std::cin >> m1;

        std::cout << "Input second matrix (rows and columns on the first line, matrix on the following lines, all separated by spaces)" << std::endl;

        std::cin >> m2;

        std::cout << "First matrix:\n" << m1 << "Second matrix: \n" << m2 << "Matrix multiplication: \n" << m1 * m2 << std::endl;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
