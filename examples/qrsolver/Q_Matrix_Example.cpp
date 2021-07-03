/**
 * @file Q_Matrix_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Example to test QMatrix class
 * @version 0.1
 * @date 2021-06-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>
#include <stdexcept>

#include "qrsolver/Q_Matrix.hpp"

int main(int argc, char const *argv[])
{
    QMatrix<float> A(5);
    float k = 1;

    float b[5] = {1, 2, 3, 4, 5};
    float x[5] = {0, 0, 0, 0, 0};

    try
    {
        for (int i = 0; i < 4; i++)

            for (int j = 0; j <= i+1; j++)

                A.setElement(i, j, k++);

        for (int j = 0; j < 5; j++)

            A.setElement(4, j, k++);

        std::cout << "Q Matrix" << std::endl;
        A.printMatrix();

        A.multiply(b, x);

        std::cout << "\nMultiplication Result" << std::endl;
        for (int i = 0; i < 5; i++) std::cout << x[i] << std::endl;


        A.identity();
        std::cout << "\nQ Matrix" << std::endl;
        A.printMatrix();
    }
    catch(const std::out_of_range& e)
    {
        std::cerr << e.what() << '\n';
    }

    return 0;
}
