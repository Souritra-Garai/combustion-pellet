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
	@mkdir -p $(BINDIR);
	@echo "\nLinking all...";
	$(CC) $(OBJECTS) $(CFLAGS) $(LIB) -o $(TARGET)

$(BUILDDIR)/main.o : $(SRCDIR)/main.cpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -o $(BUILDDIR)/main.o

$(BUILDDIR)/R_Matrix.o : $(SRCDIR)/R_Matrix.cpp $(INCDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling R_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/R_Matrix.cpp -o $(BUILDDIR)/R_Matrix.o

$(BUILDDIR)/Q_Matrix.o : $(SRCDIR)/Q_Matrix.cpp $(INCDIR)/Q_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Q_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Q_Matrix.cpp -o $(BUILDDIR)/Q_Matrix.o

$(BUILDDIR)/G_Matrix.o : $(SRCDIR)/G_Matrix.cpp $(INCDIR)/G_Matrix.hpp $(INCDIR)/Q_Matrix.hpp $(INCDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling G_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/G_Matrix.cpp -o $(BUILDDIR)/G_Matrix.o 

$(BUILDDIR)/QR_Solver.o : $(SRCDIR)/QR_Solver.cpp $(INCDIR)/QR_Solver.hpp $(INCDIR)/G_Matrix.hpp $(INCDIR)/Q_Matrix.hpp $(INCDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling QR_Solver...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/QR_Solver.cpp -o $(BUILDDIR)/QR_Solver.o

# Builds the example for R Matrix
R_Matrix_Example : $(BUILDDIR)/R_Matrix_Example.o $(BUILDDIR)/R_Matrix.o
	@mkdir -p $(BINDIR);
	@echo "\nBuilding R_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/R_Matrix_Example.o $(BUILDDIR)/R_Matrix.o -o $(BINDIR)/R_Matrix_Example

$(BUILDDIR)/R_Matrix_Example.o : $(EXMDIR)/R_Matrix_Example.cpp $(INCDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling R_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/R_Matrix_Example.cpp -o $(BUILDDIR)/R_Matrix_Example.o

# Builds the example for Q Matrix
Q_Matrix_Example : $(BUILDDIR)/Q_Matrix_Example.o $(BUILDDIR)/Q_Matrix.o
	@mkdir -p $(BINDIR);
	@echo "\nBuilding Q_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/Q_Matrix_Example.o $(BUILDDIR)/Q_Matrix.o -o $(BINDIR)/Q_Matrix_Example

$(BUILDDIR)/Q_Matrix_Example.o : $(EXMDIR)/Q_Matrix_Example.cpp $(INCDIR)/Q_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Q_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/Q_Matrix_Example.cpp -o $(BUILDDIR)/Q_Matrix_Example.o

# Builds the example for Givens Rotation matrix
G_Matrix_Example : $(BUILDDIR)/G_Matrix_Example.o $(BUILDDIR)/G_Matrix.o $(BUILDDIR)/Q_Matrix.o $(BUILDDIR)/R_Matrix.o
	@mkdir -p $(BINDIR);
	@echo "\nBuilding G_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/G_Matrix_Example.o $(BUILDDIR)/G_Matrix.o $(BUILDDIR)/Q_Matrix.o $(BUILDDIR)/R_Matrix.o -o $(BINDIR)/G_Matrix_Example

$(BUILDDIR)/G_Matrix_Example.o : $(EXMDIR)/G_Matrix_Example.cpp $(INCDIR)/G_Matrix.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling G_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/G_Matrix_Example.cpp -o $(BUILDDIR)/G_Matrix_Example.o

# Builds the example for QR Solver
QR_Solver_Example : $(BUILDDIR)/QR_Solver_Example.o $(BUILDDIR)/QR_Solver.o $(BUILDDIR)/G_Matrix.o $(BUILDDIR)/Q_Matrix.o $(BUILDDIR)/R_Matrix.o
	@mkdir -p $(BINDIR);
	@echo "\nBuilding QR_Solver_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/QR_Solver_Example.o $(BUILDDIR)/QR_Solver.o $(BUILDDIR)/G_Matrix.o $(BUILDDIR)/Q_Matrix.o $(BUILDDIR)/R_Matrix.o -o $(BINDIR)/QR_Solver_Example

$(BUILDDIR)/QR_Solver_Example.o : $(EXMDIR)/QR_Solver_Example.cpp $(INCDIR)/QR_Solver.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling QR_Solver_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/QR_Solver_Example.cpp -o $(BUILDDIR)/QR_Solver_Example.o


clean:
	@echo "\nCleaning...";
	$(RM) -r $(BUILDDIR) $(BINDIR)