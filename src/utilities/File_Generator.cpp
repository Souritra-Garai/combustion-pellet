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

// #include <sys/stat.h>
// #include <filesystem>
// #include <Windows.h>

#include <time.h>
#include <string.h>

#include <iostream>

#include "utilities/File_Generator.hpp"

FileGenerator::FileGenerator()
{
    time_t raw_time;
    struct tm * time_info;

    time(&raw_time);
    time_info = localtime(&raw_time);

	std::cout << "Fetched time...\n";

    std::string current_time = asctime(time_info);
    current_time.erase(current_time.end() - 1);

	std::cout << "Formatted time string...\n";

	_folder_name.clear();
	_folder_name.append("solutions\\");
	// _folder_name.append(current_time);

	// std::cout << "Added string to folder...\n" << _folder_name <<  "\n";

	// std::filesystem::create_directory(_folder_name.c_str());

	// std::cout << CreateDirectory(_folder_name.c_str(), NULL);

	// std::cout << "Error :" << GetLastError() << "\n";

	// std::cout << "Created folder.\n";
    // mkdir(_folder_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

std::ofstream FileGenerator::getCSVFile(std::string file_name)
{
    std::string file;
    
	file.clear();
	file.append(_folder_name);
	
	file.append("/");
    file.append(file_name);
	file.append(".csv");

    return std::ofstream(file);
}

// std::ofstream FileGenerator::getCSVFile(std::string file_name, std::string folder)
// {
//     std::string file;
    
// 	file.clear();
// 	file.append(_folder_name);
	
// 	file.append("/");
// 	file.append(folder);

// 	// std::filesystem::create_directory(_folder_name.c_str());
// 	// mkdir(file.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

//    	file.append("/");
//     file.append(file_name);
// 	file.append(".csv");

//     return std::ofstream(file);
// }

std::ofstream FileGenerator::getTXTFile(std::string file_name)
{
    std::string file;
    
	file.clear();
	file.append(_folder_name);
	
	file.append("/");
    file.append(file_name);
	file.append(".txt");

    return std::ofstream(file);
}

// std::ofstream FileGenerator::getTXTFile(std::string file_name, std::string folder)
// {
//     std::string file;
    
//     file.clear();
// 	file.append(_folder_name);
	
// 	file.append("/");
// 	file.append(folder);

// 	std::filesystem::create_directory(_folder_name.c_str());
// 	// mkdir(file.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

//    	file.append("/");
//     file.append(file_name);
// 	file.append(".txt");

//     return std::ofstream(file);
// }