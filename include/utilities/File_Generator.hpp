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
    
        std::string _folder_name;
        
    public:
        
        FileGenerator();

        std::ofstream getCSVFile(std::string name);
        // std::ofstream getCSVFile(std::string name, std::string folder);

        std::ofstream getTXTFile(std::string name);        
        // std::ofstream getTXTFile(std::string name, std::string folder);        
};

#endif