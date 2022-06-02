/**
 * @file File_Utilities.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __FILE_GENERATOR__
#define __FILE_GENERATOR__

#include <stdio.h>
#include <string>

class FileGenerator
{
    private:
    
        std::string _folder_name;
        
    public:
        
        FileGenerator();

        FILE * getCSVFile(std::string name);
        FILE * getCSVFile(std::string name, std::string folder);

        FILE * getTXTFile(std::string name);        
        FILE * getTXTFile(std::string name, std::string folder);        
};

#endif