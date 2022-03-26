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
SPCDIR := species
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

SPECIES_HPP := $(INCDIR)/$(TRPDIR)/Species.cuh $(INCDIR)/$(TRPDIR)/Ideal_Gas.cuh $(INCDIR)/$(TRPDIR)/Phase.cuh $(INCDIR)/$(TRPDIR)/Enthalpy.cuh $(INCDIR)/$(TRPDIR)/Thermal_Conductivity.cuh $(INCDIR)/$(SPCDIR)/Argon.cuh $(INCDIR)/$(SPCDIR)/Aluminium.cuh $(INCDIR)/$(SPCDIR)/Nickel.cuh $(INCDIR)/$(SPCDIR)/NickelAluminide.cuh

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

$(BUILDDIR)/main.o : $(SRCDIR)/main.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPECIES_HPP) $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR);
	@echo "Compiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -I /usr/include/boost -o $(BUILDDIR)/main.o

Thermo-Physical_Properties : $(BUILDDIR)/Thermo-Physical_Properties.o $(BUILDDIR)/$(UTILDIR)/File_Generator.o
	@mkdir -p $(BINDIR);
	@echo "Linking all files for Thermo-Physical_Properties...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/Thermo-Physical_Properties.o $(BUILDDIR)/$(UTILDIR)/File_Generator.o -lboost_program_options -o $(BINDIR)/Thermo-Physical_Properties

$(BUILDDIR)/Thermo-Physical_Properties.o : $(SRCDIR)/Thermo-Physical_Properties.cu $(SPECIES_HPP)
	@mkdir -p $(BUILDDIR);
	@echo "Compiling Thermo-Physical_Properties...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Thermo-Physical_Properties.cu -I /usr/include/boost -o $(BUILDDIR)/Thermo-Physical_Properties.o

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

$(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(SPECIES_HPP) $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cpp -o $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o

Pellet-Flame-Propagation_Example : $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS) -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example

$(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPECIES_HPP) $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cpp -o $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o

$(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPECIES_HPP) $(INCDIR)/$(UTILDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp -I /usr/include/boost -o $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o

Pellet-Flame-Propagation-Diffusion_Example : $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTILOBJS) -lboost_program_options -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example

# ----------------------------------------------------------------------------------------------------------
# Building Thermo-Physical Properties Examples
Enthalpy_Example : $(BUILDDIR)/$(TRPDIR)/Enthalpy_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Enthalpy_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Enthalpy_Example.o -o $(BINDIR)/$(TRPDIR)/Enthalpy_Example

$(BUILDDIR)/$(TRPDIR)/Enthalpy_Example.o : $(EXMDIR)/$(TRPDIR)/Enthalpy_Example.cu $(INCDIR)/$(TRPDIR)/Enthalpy.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Enthalpy_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Enthalpy_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Enthalpy_Example.o

Thermal_Conductivity_Example : $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Thermal_Conductivity_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Example.o -o $(BINDIR)/$(TRPDIR)/Thermal_Conductivity_Example

$(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Example.o : $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Example.cu $(INCDIR)/$(TRPDIR)/Thermal_Conductivity.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Example.o

Ideal_Gas_Example : $(BUILDDIR)/$(TRPDIR)/Ideal_Gas_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Ideal_Gas_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Ideal_Gas_Example.o -o $(BINDIR)/$(TRPDIR)/Ideal_Gas_Example

$(BUILDDIR)/$(TRPDIR)/Ideal_Gas_Example.o : $(EXMDIR)/$(TRPDIR)/Ideal_Gas_Example.cu $(INCDIR)/$(TRPDIR)/Ideal_Gas.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Ideal_Gas_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Ideal_Gas_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Ideal_Gas_Example.o

Phase_Example : $(BUILDDIR)/$(TRPDIR)/Phase_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Phase_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Phase_Example.o -o $(BINDIR)/$(TRPDIR)/Phase_Example

$(BUILDDIR)/$(TRPDIR)/Phase_Example.o : $(EXMDIR)/$(TRPDIR)/Phase_Example.cu $(INCDIR)/$(TRPDIR)/Phase.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Phase_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Phase_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Phase_Example.o

Species_Example : $(BUILDDIR)/$(TRPDIR)/Species_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Species_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Species_Example.o -o $(BINDIR)/$(TRPDIR)/Species_Example

$(BUILDDIR)/$(TRPDIR)/Species_Example.o : $(EXMDIR)/$(TRPDIR)/Species_Example.cu $(INCDIR)/$(TRPDIR)/Species.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Species_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Species_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Species_Example.o

Thermal_Conductivity_Pellet_Example : $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Thermal_Conductivity_Pellet_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(UTILOBJS) -o $(BINDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example

$(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.cpp $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SPECIES_HPP) $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o

Adiabatic_Combustion_Temperature : $(BUILDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o $(UTILOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Adiabatic_Combustion_Temperature...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o $(UTILOBJS) -o $(BINDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature

$(BUILDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o : $(EXMDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.cpp $(SPECIES_HPP) $(INCDIR)/$(UTILDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Adiabatic_Combustion_Temperature...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.cpp -o $(BUILDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o

Core-Shell-Combustion-Particle_Example : $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Core-Shell-Combustion-Particle_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o -o $(BINDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example

$(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o : $(EXMDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.cpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(SPECIES_HPP)
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Core-Shell-Combustion-Particle_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.cpp -o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle_Example.o

Packed-Pellet_Example : $(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o $(BUILDDIR)/$(TRPDIR)/Packed-Pellet.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Packed-Pellet_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o $(BUILDDIR)/$(TRPDIR)/Packed-Pellet.o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o -o $(BINDIR)/$(TRPDIR)/Packed-Pellet_Example

$(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cpp $(INCDIR)/$(TRPDIR)/Packed-Pellet.hpp $(SPECIES_HPP) $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp
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