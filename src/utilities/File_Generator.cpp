/**
 * @file File_Utilities.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <sys/stat.h>

#include <time.h>
#include <string.h>

#include "utilities/File_Generator.hpp"

FileGenerator::FileGenerator()
{
    time_t raw_time;
    struct tm * time_info;

    time(&raw_time);
    time_info = localtime(&raw_time);

    char * current_time = asctime(time_info);
    current_time[strlen(current_time) - 1] = '\0';

    strcpy(folder_name, "solutions/");
    strcat(folder_name, current_time);

    mkdir(folder_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

std::ofstream FileGenerator::getCSVFile(const char * file_name)
{
    char file[100];
    
    strcpy(file, folder_name);
    strcat(file, "/");

    strcat(file, file_name);
    strcat(file, ".csv");

    return std::ofstream(file);
}

std::ofstream FileGenerator::getTXTFile(const char * file_name)
{
    char file[100];
    
    strcpy(file, folder_name);
    strcat(file, "/");

    strcat(file, file_name);
    strcat(file, ".txt");

    return std::ofstream(file);
}