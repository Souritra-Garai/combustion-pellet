# Compiler
CC := g++

# Building this file is the main objective
TARGET := bin/simulate_pellet_combustion

# Directories
SRCDIR := src
INCDIR := include
BINDIR := bin
EXMDIR := examples
BUILDDIR := build

# Finding all source and object files
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) $(LIBDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# Flags required for compiler
CFLAGS := -fopenmp
LIB := -lm
INC := -I $(INCDIR)

$(TARGET) : $(OBJECTS)
	@echo "\nLinking all...";
	$(CC) $(OBJECTS) $(CFLAGS) $(LIB) -o $(TARGET)

$(BUILDDIR)/main.o : $(SRCDIR)/main.cpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -o $(BUILDDIR)/main.o

$(BUILDDIR)/Tridiagonal_Matrix.o : $(SRCDIR)/Tridiagonal_Matrix.cpp $(INCDIR)/Tridiagonal_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Tridiagonal_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Tridiagonal_Matrix.cpp -o $(BUILDDIR)/Tridiagonal_Matrix.o

$(BUILDDIR)/Lower_Triangular_Matrix.o : $(SRCDIR)/Lower_Triangular_Matrix.cpp $(INCDIR)/Lower_Triangular_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Lower_Triangular_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Lower_Triangular_Matrix.cpp -o $(BUILDDIR)/Lower_Triangular_Matrix.o

# Builds the example for Tridiagonal Matrix
Tridiagonal_Matrix_Example : $(BUILDDIR)/Tridiagonal_Matrix_Example.o $(BUILDDIR)/Tridiagonal_Matrix.o
	@echo "\nBuilding Tridiagonal_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/Tridiagonal_Matrix_Example.o $(BUILDDIR)/Tridiagonal_Matrix.o -o $(BINDIR)/Tridiagonal_Matrix_Example

$(BUILDDIR)/Tridiagonal_Matrix_Example.o : $(EXMDIR)/Tridiagonal_Matrix_Example.cpp $(INCDIR)/Tridiagonal_Matrix.hpp $(BUILDDIR)/Tridiagonal_Matrix.o
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Tridiagonal_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/Tridiagonal_Matrix_Example.cpp -o $(BUILDDIR)/Tridiagonal_Matrix_Example.o

# Builds the example for Lower Triangular Matrix
Lower_Triangular_Matrix_Example : $(BUILDDIR)/Lower_Triangular_Matrix_Example.o $(BUILDDIR)/Lower_Triangular_Matrix.o
	@echo "\nBuilding Lower_Triangular_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/Lower_Triangular_Matrix_Example.o $(BUILDDIR)/Lower_Triangular_Matrix.o -o $(BINDIR)/Lower_Triangular_Matrix_Example

$(BUILDDIR)/Lower_Triangular_Matrix_Example.o : $(EXMDIR)/Lower_Triangular_Matrix_Example.cpp $(INCDIR)/Lower_Triangular_Matrix.hpp $(BUILDDIR)/Lower_Triangular_Matrix.o
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Lower_Triangular_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/Lower_Triangular_Matrix_Example.cpp -o $(BUILDDIR)/Lower_Triangular_Matrix_Example.o

clean:
	@echo "\nCleaning..."; 
	@echo "$(RM) -r $(BUILDDIR)/*.o $(TARGET)"; $(RM) -r $(BUILDDIR)/* $(TARGET)