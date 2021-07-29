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

#include <fstream>

class FileGenerator
{
    private:
    
        char folder_name[100];
        
    public:
        
        FileGenerator();

        std::ofstream getCSVFile(const char * name);
        std::ofstream getTXTFile(const char * name);        
};

#endif