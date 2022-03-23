# Compiler
CC := nvcc

# Building this file is the main objective
TARGET := bin/simulate_pellet_combustion

# Directories
SRCDIR := src
INCDIR := include
BINDIR := bin
EXMDIR := examples
BUILDDIR := build

LUSDIR := lusolver
QRSDIR := qrsolver
TRPDIR := thermo-physical-properties
PDEDIR := pde-problems
UTILDIR := utilities

# Finding all source and object files
SRCEXT := cu
SOURCES := $(shell find $(SRCDIR) $(LIBDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

LUSSRCS := $(shell find $(SRCDIR)/$(LUSDIR) -type f -name *.$(SRCEXT))
LUSOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(LUSSRCS:.$(SRCEXT)=.o))

QRSSRCS := $(shell find $(SRCDIR)/$(QRSDIR) -type f -name *.$(SRCEXT))
QRSOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(QRSSRCS:.$(SRCEXT)=.o))

TRPSRCS := $(shell find $(SRCDIR)/$(TRPDIR) -type f -name *.$(SRCEXT))
TRPOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(TRPSRCS:.$(SRCEXT)=.o))

PDESRCS := $(shell find $(SRCDIR)/$(PDEDIR) -type f -name *.$(SRCEXT))
PDEOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(PDESRCS:.$(SRCEXT)=.o))

UTILSRCS := $(shell find $(SRCDIR)/$(UTILDIR) -type f -name *.$(SRCEXT))
UTILOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(UTILSRCS:.$(SRCEXT)=.o))

SUBSTANCE_HPP := $(INCDIR)/$(TRPDIR)/Substance.hpp $(INCDIR)/$(TRPDIR)/IdealGas.hpp $(INCDIR)/$(TRPDIR)/Phase.hpp $(INCDIR)/$(TRPDIR)/Enthalpy.hpp $(INCDIR)/$(TRPDIR)/Thermal_Conductivity.hpp $(INCDIR)/substances/Argon.hpp $(INCDIR)/substances/Aluminium.hpp $(INCDIR)/substances/Nickel.hpp $(INCDIR)/substances/NickelAluminide.hpp

# Flags required for compiler
CFLAGS := # -fopenmp -O2
LIB := -lm
INC := -I $(INCDIR)

$(TARGET) : $(OBJECTS)
	@mkdir -p $(BINDIR);
	@echo "Linking all...";
	$(CC) $(OBJECTS) $(CFLAGS) $(LIB) -o $(TARGET)

main : $(BUILDDIR)/main.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking all files for main...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/main.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS) -lboost_program_options -o $(BINDIR)/Pellet-Flame-Propagation

$(BUILDDIR)/main.o : $(SRCDIR)/main.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SUBSTANCE_HPP) $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR);
	@echo "Compiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -I /usr/include/boost -o $(BUILDDIR)/main.o

# ----------------------------------------------------------------------------------------------------------
# Building thermo-physical properties source files

$(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o : $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(TRPDIR)/Substance.hpp $(SRCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.cpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Core-Shell-Combustion-Particle...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.cpp -o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o

$(BUILDDIR)/$(TRPDIR)/Packed-Pellet.o : $(INCDIR)/$(TRPDIR)/Packed-Pellet.hpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(TRPDIR)/Substance.hpp $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SRCDIR)/$(TRPDIR)/Packed-Pellet.cpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Packed-Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Packed-Pellet.cpp -o $(BUILDDIR)/$(TRPDIR)/Packed-Pellet.o

$(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o : $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SRCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.cpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.cpp -o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o

# ----------------------------------------------------------------------------------------------------------
# Building pde problems source files

$(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o : $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(LUSDIR)/LU_Solver.hpp $(SRCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Core-Shell-Diffusion...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cpp -o $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o

$(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation.o : $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Packed-Pellet.hpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(INCDIR)/$(LUSDIR)/LU_Solver.hpp $(SRCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.cpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.cpp -o $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation.o

# ----------------------------------------------------------------------------------------------------------
# Building utilities

$(BUILDDIR)/$(UTILDIR)/Keyboard_Interrupt.o : $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(SRCDIR)/$(UTILDIR)/Keyboard_Interrupt.cpp
	@mkdir -p $(BUILDDIR)/$(UTILDIR);
	@echo "Compiling Keyboard_Interrupt..."
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(UTILDIR)/Keyboard_Interrupt.cpp -o $(BUILDDIR)/$(UTILDIR)/Keyboard_Interrupt.o

$(BUILDDIR)/$(UTILDIR)/File_Generator.o : $(INCDIR)/$(UTILDIR)/File_Generator.hpp $(SRCDIR)/$(UTILDIR)/File_Generator.cpp
	@mkdir -p $(BUILDDIR)/$(UTILDIR);
	@echo "Compiling File_Generator..."
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(UTILDIR)/File_Generator.cpp -o $(BUILDDIR)/$(UTILDIR)/File_Generator.o

# ----------------------------------------------------------------------------------------------------------
# Building PDE Problem Examples

Core-Shell-Diffusion_Example : $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(LUSOBJS) $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(LUSOBJS) $(UTILOBJS) -o $(BINDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example

$(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(SUBSTANCE_HPP) $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cpp -o $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o

Pellet-Flame-Propagation_Example : $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS) -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example

$(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SUBSTANCE_HPP) $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cpp -o $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o

$(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SUBSTANCE_HPP) $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp -I /usr/include/boost -o $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o

Pellet-Flame-Propagation-Diffusion_Example : $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS) -lboost_program_options -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example

# ----------------------------------------------------------------------------------------------------------
# Building Thermo-Physical Properties Examples
Thermal_Conductivity_Pellet_Example : $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Thermal_Conductivity_Pellet_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(UTILOBJS) -o $(BINDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example

$(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.cpp $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SUBSTANCE_HPP) $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o

Adiabatic_Combustion_Temperature : $(BUILDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Adiabatic_Combustion_Temperature...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o $(UTILOBJS) -o $(BINDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature

$(BUILDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o : $(EXMDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.cpp $(SUBSTANCE_HPP) $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Adiabatic_Combustion_Temperature...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.cpp -o $(BUILDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o

Phase_Example : $(BUILDDIR)/$(TRPDIR)/Phase_Example.o $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Phase_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Phase_Example.o $(UTILOBJS) -o $(BINDIR)/$(TRPDIR)/Phase_Example

$(BUILDDIR)/$(TRPDIR)/Phase_Example.o : $(EXMDIR)/$(TRPDIR)/Phase_Example.cpp $(INCDIR)/$(TRPDIR)/Phase.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Phase_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Phase_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Phase_Example.o

Substance_Example : $(BUILDDIR)/$(TRPDIR)/Substance_Example.o $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Substance_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Substance_Example.o $(UTILOBJS) -o $(BINDIR)/$(TRPDIR)/Substance_Example

$(BUILDDIR)/$(TRPDIR)/Substance_Example.o : $(EXMDIR)/$(TRPDIR)/Substance_Example.cpp $(SUBSTANCE_HPP) $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Substance_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Substance_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Substance_Example.o

Core-Shell-Combustion-Particle_Example : $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Core-Shell-Combustion-Particle_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o -o $(BINDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example

$(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o : $(EXMDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.cpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(SUBSTANCE_HPP)
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Core-Shell-Combustion-Particle_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o

Packed-Pellet_Example : $(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o $(BUILDDIR)/$(TRPDIR)/Packed-Pellet.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Packed-Pellet_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o $(BUILDDIR)/$(TRPDIR)/Packed-Pellet.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o -o $(BINDIR)/$(TRPDIR)/Packed-Pellet_Example

$(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cpp $(INCDIR)/$(TRPDIR)/Packed-Pellet.hpp $(SUBSTANCE_HPP) $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Packed-Pellet_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o

# ----------------------------------------------------------------------------------------------------------
# Builds the example for LU Solver
LU_Solver_Example : $(BUILDDIR)/$(LUSDIR)/LU_Solver_Example.o
	@mkdir -p $(BINDIR)/$(LUSDIR);
	@echo "Building LU_Solver_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(LUSDIR)/LU_Solver_Example.o $(LUSOBJS) -o $(BINDIR)/$(LUSDIR)/LU_Solver_Example

$(BUILDDIR)/$(LUSDIR)/LU_Solver_Example.o : $(EXMDIR)/$(LUSDIR)/LU_Solver_Example.cu $(INCDIR)/$(LUSDIR)/LU_Solver.cuh
	@mkdir -p $(BUILDDIR)/$(LUSDIR);
	@echo "Compiling LU_Solver_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(LUSDIR)/LU_Solver_Example.cu -o $(BUILDDIR)/$(LUSDIR)/LU_Solver_Example.o

# Builds the example for Tridiagonal Matrix
Tridiagonal_Matrix_Example : $(BUILDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o 
	@mkdir -p $(BINDIR)/$(LUSDIR);
	@echo "Building Tridiagonal_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o -o $(BINDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example

$(BUILDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o : $(EXMDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.cu $(INCDIR)/$(LUSDIR)/Tridiagonal_Matrix.cuh
	@mkdir -p $(BUILDDIR)/$(LUSDIR);
	@echo "Compiling Tridiagonal_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.cu -o $(BUILDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o

# ----------------------------------------------------------------------------------------------------------
# Building QR Factorization solver

$(BUILDDIR)/$(QRSDIR)/R_Matrix.o : $(SRCDIR)/$(QRSDIR)/R_Matrix.cpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "Compiling R_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/R_Matrix.cpp -o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o

$(BUILDDIR)/$(QRSDIR)/Q_Matrix.o : $(SRCDIR)/$(QRSDIR)/Q_Matrix.cpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "Compiling Q_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/Q_Matrix.cpp -o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o

$(BUILDDIR)/$(QRSDIR)/G_Matrix.o : $(SRCDIR)/$(QRSDIR)/G_Matrix.cpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "Compiling G_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/G_Matrix.cpp -o $(BUILDDIR)/$(QRSDIR)/G_Matrix.o 

$(BUILDDIR)/$(QRSDIR)/QR_Solver.o : $(SRCDIR)/$(QRSDIR)/QR_Solver.cpp $(INCDIR)/$(QRSDIR)/QR_Solver.hpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "Compiling QR_Solver...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/QR_Solver.cpp -o $(BUILDDIR)/$(QRSDIR)/QR_Solver.o

# Builds the example for R Matrix
R_Matrix_Example : $(BUILDDIR)/$(QRSDIR)/R_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "Building R_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/R_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o -o $(BINDIR)/$(QRSDIR)/R_Matrix_Example

$(BUILDDIR)/$(QRSDIR)/R_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/R_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "Compiling R_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/R_Matrix_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/R_Matrix_Example.o

# Builds the example for Q Matrix
Q_Matrix_Example : $(BUILDDIR)/$(QRSDIR)/Q_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "Building Q_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/Q_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o -o $(BINDIR)/$(QRSDIR)/Q_Matrix_Example

$(BUILDDIR)/$(QRSDIR)/Q_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/Q_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "Compiling Q_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/Q_Matrix_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/Q_Matrix_Example.o

# Builds the example for Givens Rotation matrix
G_Matrix_Example : $(BUILDDIR)/$(QRSDIR)/G_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/G_Matrix.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "Building G_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/G_Matrix_Example.o $(BUILDDIR)/$(QRSDIR)/G_Matrix.o $(BUILDDIR)/$(QRSDIR)/Q_Matrix.o $(BUILDDIR)/$(QRSDIR)/R_Matrix.o -o $(BINDIR)/$(QRSDIR)/G_Matrix_Example

$(BUILDDIR)/$(QRSDIR)/G_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/G_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "Compiling G_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/G_Matrix_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/G_Matrix_Example.o

# Builds the example for QR Solver
QR_Solver_Example : $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o $(QRSOBJS)
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "Building QR_Solver_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o $(QRSOBJS) -o $(BINDIR)/$(QRSDIR)/QR_Solver_Example

$(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o : $(EXMDIR)/$(QRSDIR)/QR_Solver_Example.cpp $(INCDIR)/$(QRSDIR)/QR_Solver.hpp
	@mkdir -p $(BUILDDIR)/$(QRSDIR);
	@echo "Compiling QR_Solver_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/QR_Solver_Example.cpp -o $(BUILDDIR)/$(QRSDIR)/QR_Solver_Example.o

# ----------------------------------------------------------------------------------------------------------------------
# Clean everything

clean:
	@echo "Cleaning...";
	$(RM) -r $(BUILDDIR) $(BINDIR)