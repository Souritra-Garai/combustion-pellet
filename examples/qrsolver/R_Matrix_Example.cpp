/**
 * @file R_Matrix_Example.cpp
 * @author Souritra Garai (you@domain.com)
 * @brief An example cpp program to test the RMatrix class
 * @version 0.1
 * @date 2021-06-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <stdio.h>
#include <iostream>
#include <exception>

#include "qrsolver/R_Matrix.hpp"

int main(int argc, char const *argv[])
{
    try
    {   
        RMatrix<float> A(5);

        int k = 1;
        A.setElement(0, 0, k++);
        A.setElement(0, 1, k++);

        for (int i = 1; i < 4; i++)
        {
            A.setElement(i, i-1, k++);
            A.setElement(i, i, k++);
            A.setElement(i, i+1, k++);
        }

        A.setElement(4, 3, k++);
        A.setElement(4, 4, k++);

        printf("R Matrix\n");
        A.printMatrix();

        printf("\n\nFlat Representation\n");
        A.print();

        RMatrix<long double> B(5);

        k = 1;

        for (int i = 0; i < 3; i++)
        {
            B.setElement(i, i, k++);
            B.setElement(i, i+1, k++);
            B.setElement(i, i+2, k++);
        }

        B.setElement(3, 3, k++);
        B.setElement(3, 4, k++);
        
        B.setElement(4, 4, k++);

        // B.setElement(1, 4, k++);

        printf("\n\n\nR Matrix\n");
        B.printMatrix();

        printf("\n\nFlat Representation\n");
        B.print();
    }

    catch (std::out_of_range error)
    {
        // printf(error.what());
        std::cerr << error.what() << std::endl;
    }

    return 0;
}
