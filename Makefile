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
QRSDIR := qrsolver

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

# ----------------------------------------------------------------------------------------------------------
# Building QR Factorization solver

$(BUILDDIR)/$(QRSDIR)/R_Matrix.o : $(SRCDIR)/$(QRSDIR)/R_Matrix.cpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling R_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/R_Matrix.cpp -o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o

$(BUILDDIR)/$(QRSDIR)/Q_Matrix.o : $(SRCDIR)/$(QRSDIR)/Q_Matrix.cpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling Q_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/Q_Matrix.cpp -o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o

$(BUILDDIR)/$(QRSDIR)/G_Matrix.o : $(SRCDIR)/$(QRSDIR)/G_Matrix.cpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling G_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/G_Matrix.cpp -o $(BUILDDIR)/$(QRSDIR)/G_Matrix.o 

$(BUILDDIR)/$(QRSDIR)/QR_Solver.o : $(SRCDIR)/$(QRSDIR)/QR_Solver.cpp $(INCDIR)/$(QRSDIR)/QR_Solver.hpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling QR_Solver...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/QR_Solver.cpp -o $(BUILDDIR)/$(QRSDIR)/QR_Solver.o

# Builds the example for R Matrix
R_Matrix_Example : $(BUILDDIR)/$(QRSDIR)/R_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "\nBuilding R_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/R_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o -o $(BINDIR)/$(QRSDIR)/R_Matrix_Example

$(BUILDDIR)/$(QRSDIR)/R_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/R_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling R_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/R_Matrix_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/R_Matrix_Example.o

# Builds the example for Q Matrix
Q_Matrix_Example : $(BUILDDIR)/$(QRSDIR)/Q_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "\nBuilding Q_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/Q_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o -o $(BINDIR)/$(QRSDIR)/Q_Matrix_Example

$(BUILDDIR)/$(QRSDIR)/Q_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/Q_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling Q_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/Q_Matrix_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/Q_Matrix_Example.o

# Builds the example for Givens Rotation matrix
G_Matrix_Example : $(BUILDDIR)/$(QRSDIR)/G_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/G_Matrix.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "\nBuilding G_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/G_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/G_Matrix.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o -o $(BINDIR)/$(QRSDIR)/G_Matrix_Example

$(BUILDDIR)/$(QRSDIR)/G_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/G_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling G_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/G_Matrix_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/G_Matrix_Example.o

# Builds the example for QR Solver
QR_Solver_Example : $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o $(BUILDDIR)/$(QRSDIR)/QR_Solver.o $(BUILDDIR)/$(QRSDIR)/G_Matrix.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "\nBuilding QR_Solver_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o $(BUILDDIR)/$(QRSDIR)/QR_Solver.o $(BUILDDIR)/$(QRSDIR)/G_Matrix.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o -o $(BINDIR)/$(QRSDIR)/QR_Solver_Example

$(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o : $(EXMDIR)/$(QRSDIR)/QR_Solver_Example.cpp $(INCDIR)/$(QRSDIR)/QR_Solver.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling QR_Solver_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/QR_Solver_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o


clean:
	@echo "\nCleaning...";
	$(RM) -r $(BUILDDIR) $(BINDIR)