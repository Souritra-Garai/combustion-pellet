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
TRPDIR := thermo-physical-properties
PDEDIR := pde-problems
UTILDIR := utilities

# Finding all source and object files
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) $(LIBDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

QRSSRCS := $(shell find $(SRCDIR)/$(QRSDIR) -type f -name *.$(SRCEXT))
QRSOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(QRSSRCS:.$(SRCEXT)=.o))

TRPSRCS := $(shell find $(SRCDIR)/$(TRPDIR) -type f -name *.$(SRCEXT))
TRPOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(TRPSRCS:.$(SRCEXT)=.o))

PDESRCS := $(shell find $(SRCDIR)/$(PDEDIR) -type f -name *.$(SRCEXT))
PDEOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(PDESRCS:.$(SRCEXT)=.o))

UTILSRCS := $(shell find $(SRCDIR)/$(UTILDIR) -type f -name *.$(SRCEXT))
UTILOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(UTILSRCS:.$(SRCEXT)=.o))

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

$(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o : $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(TRPDIR)/Substance.hpp $(SRCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.cpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "\nCompiling Core-Shell-Combustion-Particle...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.cpp -o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o

$(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o : $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(QRSDIR)/QR_Solver.hpp $(SRCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "\nCompiling Core-Shell-Diffusion...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cpp -o $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o

# ----------------------------------------------------------------------------------------------------------
# Building utilities

$(BUILDDIR)/$(UTILDIR)/Keyboard_Interrupt.o : $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(SRCDIR)/$(UTILDIR)/Keyboard_Interrupt.cpp
	@mkdir -p $(BUILDDIR)/$(UTILDIR);
	@echo "\nCompiling Keyboard_Interrupt..."
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(UTILDIR)/Keyboard_Interrupt.cpp -o $(BUILDDIR)/$(UTILDIR)/Keyboard_Interrupt.o

$(BUILDDIR)/$(UTILDIR)/File_Generator.o : $(INCDIR)/$(UTILDIR)/File_Generator.hpp $(SRCDIR)/$(UTILDIR)/File_Generator.cpp
	@mkdir -p $(BUILDDIR)/$(UTILDIR);
	@echo "\nCompiling File_Generator..."
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(UTILDIR)/File_Generator.cpp -o $(BUILDDIR)/$(UTILDIR)/File_Generator.o

# ----------------------------------------------------------------------------------------------------------
# Building PDE Problems Examples

Core-Shell-Diffusion_Example : $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(QRSOBJS) $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "\nLinking Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(QRSOBJS) $(UTILOBJS) -o $(BINDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example

$(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(TRPDIR)/Substance.hpp $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "\nCompiling Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cpp -o $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o

# ----------------------------------------------------------------------------------------------------------
# Building Thermo-Physical Properties Examples

Substance_Example : $(BUILDDIR)/$(TRPDIR)/Substance_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "\nLinking Substance_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Substance_Example.o -o $(BINDIR)/$(TRPDIR)/Substance_Example

$(BUILDDIR)/$(TRPDIR)/Substance_Example.o : $(EXMDIR)/$(TRPDIR)/Substance_Example.cpp $(INCDIR)/$(TRPDIR)/Substance.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "\nCompiling Substance_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Substance_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Substance_Example.o

Core-Shell-Combustion-Particle_Example : $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "\nLinking Core-Shell-Combustion-Particle_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o -o $(BINDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example

$(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o : $(EXMDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.cpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "\nCompiling Core-Shell-Combustion-Particle_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o

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
QR_Solver_Example : $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o $(QRSOBJS)
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "\nBuilding QR_Solver_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o $(QRSOBJS) -o $(BINDIR)/$(QRSDIR)/QR_Solver_Example

$(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o : $(EXMDIR)/$(QRSDIR)/QR_Solver_Example.cpp $(INCDIR)/$(QRSDIR)/QR_Solver.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "\nCompiling QR_Solver_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/QR_Solver_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o

# ----------------------------------------------------------------------------------------------------------------------
# Clean everything

clean:
	@echo "\nCleaning...";
	$(RM) -r $(BUILDDIR) $(BINDIR)