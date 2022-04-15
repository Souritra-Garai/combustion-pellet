# Compiler
CC := nvcc

# Directories
SRCDIR := src
INCDIR := include
BINDIR := bin
EXMDIR := examples
BLDDIR := build

LUSDIR := lu-solver
TRPDIR := thermo-physical-properties
PDEDIR := pde-problems
SPCDIR := species
UTLDIR := utilities

# UTLHPP	:= $(shell find $(SRCDIR)/$(UTLDIR) -type f -name *.hpp)
# UTLSRCS := $(shell find $(SRCDIR)/$(UTLDIR) -type f -name *.cpp)
# UTLOBJS := $(patsubst $(SRCDIR)/%, $(BLDDIR)/%, $(UTLSRCS:.cpp=.o))

# SPCHPP := $(shell find $(INCDIR)/$(SPCDIR) -type f -name *.cuh)
# TRPHPP := $(shell find $(INCDIR)/$(TRPDIR) -type f -name *.cuh)
# LUSHPP := $(shell find $(INCDIR)/$(LUSDIR) -type f -name *.cuh)
# PDEHPP := $(shell find $(INCDIR)/$(PDEDIR) -type f -name *.cuh)

EXEEXT :=.exe

# Flags required for compiler
CFLAGS := -fopenmp -O2
NFLAGS := -arch=sm_75
LBOOST := D:\Souritra\boost\lib\libboost_program_options-vc142-mt-s-x64-1_79.lib
IBOOST := -I D:\Souritra\boost\include
LIB := -lm
INC := -I $(INCDIR)

main : $(BLDDIR)/main.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking all files for main...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/main.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS) -lboost_program_options -o $(BINDIR)/Pellet-Flame-Propagation

$(BLDDIR)/main.o : $(SRCDIR)/main.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPCHPP) $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR);
	@echo "Compiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -I /usr/include/boost -o $(BLDDIR)/main.o

Thermo-Physical_Properties : $(SRCDIR)/Thermo-Physical_Properties.cu $(SRCDIR)/$(UTLDIR)/File_Generator.cpp
	$(CC) -o $(BINDIR)/$(TRPDIR)/Thermo-Physical_Properties$(EXEEXT) "-std=c++17" $(NFLAGS) $(INC) $(IBOOST) $(SRCDIR)/Thermo-Physical_Properties.cu $(SRCDIR)/$(UTLDIR)/File_Generator.cpp $(LBOOST)

# $(BLDDIR)/Thermo-Physical_Properties.o : $(SRCDIR)/Thermo-Physical_Properties.cu $(SPCHPP) $(TRPHPP) $(INCDIR)/$(UTLDIR)/File_Generator.hpp
# 	$(CC) -o $(BLDDIR)/Thermo-Physical_Properties.o $(INC) $(IBOOST) $(SRCDIR)/Thermo-Physical_Properties.cu D:\Souritra\boost\lib\libboost_program_options-vc142-mt-s-x64-1_79.lib

Core-Shell-Diffusion : $(BLDDIR)/Core-Shell-Diffusion.o $(UTLOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking all files for Core-Shell-Diffusion...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/Core-Shell-Diffusion.o $(UTLOBJS) -lboost_program_options -o $(BINDIR)/Core-Shell-Diffusion

$(BLDDIR)/Core-Shell-Diffusion.o : $(SRCDIR)/Core-Shell-Diffusion.cu $(SPCHPP) $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cuh $(INCDIR)/$(TRPDIR)/Core-Shell-Particle.cuh $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.cuh
	@mkdir -p $(BLDDIR);
	@echo "Compiling Core-Shell-Diffusion...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Core-Shell-Diffusion.cu -I /usr/include/boost -o $(BLDDIR)/Core-Shell-Diffusion.o

# ----------------------------------------------------------------------------------------------------------
# Building thermo-physical properties source files

$(BLDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o : $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(TRPDIR)/Substance.hpp $(SRCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.cpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Core-Shell-Combustion-Particle...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.cpp -o $(BLDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o

$(BLDDIR)/$(TRPDIR)/Packed-Pellet.o : $(INCDIR)/$(TRPDIR)/Packed-Pellet.hpp $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(TRPDIR)/Substance.hpp $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SRCDIR)/$(TRPDIR)/Packed-Pellet.cpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Packed-Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Packed-Pellet.cpp -o $(BLDDIR)/$(TRPDIR)/Packed-Pellet.o

$(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o : $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SRCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.cpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.cpp -o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o

# ----------------------------------------------------------------------------------------------------------
# Building pde problems source files

$(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation.o : $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Packed-Pellet.hpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(INCDIR)/$(LUSDIR)/LU_Solver.hpp $(SRCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.cpp
	@mkdir -p $(BLDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.cpp -o $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation.o

# ----------------------------------------------------------------------------------------------------------
# Building utilities

$(BLDDIR)/$(UTLDIR)/Keyboard_Interrupt.o : $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(SRCDIR)/$(UTLDIR)/Keyboard_Interrupt.cpp
	@mkdir -p $(BLDDIR)/$(UTLDIR);
	@echo "Compiling Keyboard_Interrupt..."
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(UTLDIR)/Keyboard_Interrupt.cpp -o $(BLDDIR)/$(UTLDIR)/Keyboard_Interrupt.o

$(BLDDIR)/$(UTLDIR)/File_Generator.o : $(INCDIR)/$(UTLDIR)/File_Generator.hpp $(SRCDIR)/$(UTLDIR)/File_Generator.cpp
	@mkdir -p $(BLDDIR)/$(UTLDIR);
	@echo "Compiling File_Generator..."
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(UTLDIR)/File_Generator.cpp -o $(BLDDIR)/$(UTLDIR)/File_Generator.o

# ----------------------------------------------------------------------------------------------------------
# Building PDE Problem Examples

Core-Shell-Diffusion_Example : $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o $(UTLOBJS) -o $(BINDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example

$(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cu $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cuh $(INCDIR)/$(TRPDIR)/Core-Shell-Particle.cuh $(SPCHPP) $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR)/$(PDEDIR);
	@echo "Compiling Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cu -o $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o

Pellet-Flame-Propagation_Example : $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS) -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example

$(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPCHPP) $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cpp -o $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o

$(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPCHPP) $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp -I /usr/include/boost -o $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o

Pellet-Flame-Propagation-Diffusion_Example : $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS) -lboost_program_options -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example

# ----------------------------------------------------------------------------------------------------------
# Building Thermo-Physical Properties Examples
Enthalpy_Example : $(EXMDIR)/$(TRPDIR)/Enthalpy_Example.cu $(INCDIR)/$(TRPDIR)/Enthalpy.cuh
	$(CC) -o $(BINDIR)/$(TRPDIR)/Enthalpy_Example$(EXEEXT) $(NFLAGS) $(INC) $(EXMDIR)/$(TRPDIR)/Enthalpy_Example.cu

Thermal_Conductivity_Example : $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Example.cu $(INCDIR)/$(TRPDIR)/Thermal_Conductivity.cuh
	$(CC) -o $(BINDIR)/$(TRPDIR)/Thermal_Conductivity_Example$(EXEEXT) $(NFLAGS) $(INC) $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Example.cu

Ideal_Gas_Example : $(EXMDIR)/$(TRPDIR)/Ideal_Gas_Example.cu $(INCDIR)/$(TRPDIR)/Ideal_Gas.cuh $(INCDIR)/$(TRPDIR)/Enthalpy.cuh $(INCDIR)/$(TRPDIR)/Thermal_Conductivity.cuh
	$(CC) -o $(BINDIR)/$(TRPDIR)/Ideal_Gas_Example$(EXEEXT) $(NFLAGS) $(INC) $(EXMDIR)/$(TRPDIR)/Ideal_Gas_Example.cu

Phase_Example : $(EXMDIR)/$(TRPDIR)/Phase_Example.cu $(INCDIR)/$(TRPDIR)/Phase.cuh $(INCDIR)/$(TRPDIR)/Enthalpy.cuh $(INCDIR)/$(TRPDIR)/Thermal_Conductivity.cuh
	$(CC) -o $(BINDIR)/$(TRPDIR)/Phase_Example$(EXEEXT) $(NFLAGS) $(INC) $(EXMDIR)/$(TRPDIR)/Phase_Example.cu

Species_Example : $(EXMDIR)/$(TRPDIR)/Species_Example.cu $(INCDIR)/$(TRPDIR)/Species.cuh $(INCDIR)/$(TRPDIR)/Phase.cuh $(INCDIR)/$(TRPDIR)/Enthalpy.cuh $(INCDIR)/$(TRPDIR)/Thermal_Conductivity.cuh
	$(CC) -o $(BINDIR)/$(TRPDIR)/Species_Example$(EXEEXT) $(NFLAGS) $(INC) $(EXMDIR)/$(TRPDIR)/Species_Example.cu

Core-Shell-Particle_Example : $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Core-Shell-Particle_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o -o $(BINDIR)/$(TRPDIR)/Core-Shell-Particle_Example

$(BLDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o : $(EXMDIR)/$(TRPDIR)/Core-Shell-Particle_Example.cu $(INCDIR)/$(TRPDIR)/Core-Shell-Particle.cuh $(SPCHPP)
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Core-Shell-Particle_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Core-Shell-Particle_Example.cu -o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o

Thermal_Conductivity_Pellet_Example : $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Thermal_Conductivity_Pellet_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(UTLOBJS) -o $(BINDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example

$(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.cpp $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SPCHPP) $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.cpp -o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o

Adiabatic_Combustion_Temperature : $(BLDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Adiabatic_Combustion_Temperature...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o $(UTLOBJS) -o $(BINDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature

$(BLDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o : $(EXMDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.cpp $(SPCHPP) $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Adiabatic_Combustion_Temperature...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.cpp -o $(BLDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o

Packed-Pellet_Example : $(BLDDIR)/$(TRPDIR)/Packed-Pellet_Example.o $(BLDDIR)/$(TRPDIR)/Packed-Pellet.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Packed-Pellet_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Packed-Pellet_Example.o $(BLDDIR)/$(TRPDIR)/Packed-Pellet.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o -o $(BINDIR)/$(TRPDIR)/Packed-Pellet_Example

$(BLDDIR)/$(TRPDIR)/Packed-Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cpp $(INCDIR)/$(TRPDIR)/Packed-Pellet.hpp $(SPCHPP) $(INCDIR)/$(TRPDIR)/Core-Shell-Combustion-Particle.hpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Packed-Pellet_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cpp -o $(BLDDIR)/$(TRPDIR)/Packed-Pellet_Example.o

# ----------------------------------------------------------------------------------------------------------
# Builds the example for LU Solver
LU_Solver_Example : $(EXMDIR)/$(LUSDIR)/LU_Solver_Example.cu $(INCDIR)/$(LUSDIR)/LU_Solver.cuh $(INCDIR)/$(LUSDIR)/Tridiagonal_Matrix.cuh
	$(CC) -o $(BINDIR)/$(LUSDIR)/LU_Solver_Example$(EXEEXT) $(NFLAGS) $(INC) $(EXMDIR)/$(LUSDIR)/LU_Solver_Example.cu

# Builds the example for Tridiagonal Matrix
Tridiagonal_Matrix_Example : $(EXMDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.cu $(INCDIR)/$(LUSDIR)/Tridiagonal_Matrix.cuh
	$(CC) -o $(BINDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example$(EXEEXT) $(NFLAGS) $(INC) $(EXMDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.cu

# ----------------------------------------------------------------------------------------------------------------------
# Clean everything

clean:
	rmdir -Recurse ./$(BLDDIR) ./$(BINDIR)