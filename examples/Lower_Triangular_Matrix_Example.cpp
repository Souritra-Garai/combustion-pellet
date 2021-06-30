/**
 * @file Lower_Triangular_Matrix_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Example to test LowerTriangularMatrix class
 * @version 0.1
 * @date 2021-06-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>
#include <stdexcept>

#include "Lower_Triangular_Matrix.hpp"

int main(int argc, char const *argv[])
{
    LowerTriangularMatrix<float> A(5);
    float k = 1;

    try
    {
        for (int i = 0; i < 5; i++)

            for (int j = 0; j <= i+1; j++)

                A.setElement(i, j, k++);

        std::cout << "Lower Triangular Matrix" << std::endl;
        A.printMatrix();
    }
    catch(const std::out_of_range& e)
    {
        std::cerr << e.what() << '\n';
    }

    return 0;
}
