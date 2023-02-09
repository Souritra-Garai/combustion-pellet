add_library(ReadData Read-Data.cpp)
add_library(FileGenerator File-Generator.cpp)
add_library(KeyboardInterrupt Keyboard-Interrupt.cpp)

add_library(ProgramOptions Program-Options.cpp)
target_link_libraries(ProgramOptions PUBLIC ReadData)

add_library(MyUtilities INTERFACE)
target_link_libraries(MyUtilities INTERFACE ReadData FileGenerator ProgramOptions KeyboardInterrupt)