# Compiler
CC := nvcc

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
UTLDIR := utilities

HPPEXT := cuh

LUSHPP := $(shell find $(INCDIR)/$(LUSDIR) -type f -name *.$(HPPEXT))
TRPHPP := $(shell find $(INCDIR)/$(TRPDIR) -type f -name *.$(HPPEXT))
PDEHPP := $(shell find $(INCDIR)/$(PDEDIR) -type f -name *.$(HPPEXT))
SPCHPP := $(shell find $(INCDIR)/$(SPCDIR) -type f -name *.$(HPPEXT))

UTLHPP := $(shell find $(INCDIR)/$(TRPDIR) -type f -name *.hpp)
UTLSRCS := $(shell find $(SRCDIR)/$(UTLDIR) -type f -name *.cpp)
UTLOBJS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(UTLSRCS:.cpp=.o))

# Flags required for compiler
INC := -I $(INCDIR)
NFLAGS := # -std c++14
IBOOST := -I /usr/include/boost
LBOOST := -lboost_program_options

Thermo-Physical_Properties : $(BUILDDIR)/Thermo-Physical_Properties.o $(BUILDDIR)/$(UTLDIR)/File_Generator.o
	@mkdir -p $(BINDIR);
	@echo "Linking all files for Thermo-Physical_Properties...";
	$(CC) $(LBOOST) $(BUILDDIR)/Thermo-Physical_Properties.o $(BUILDDIR)/$(UTLDIR)/File_Generator.o -o $(BINDIR)/Thermo-Physical_Properties

$(BUILDDIR)/Thermo-Physical_Properties.o : $(SRCDIR)/Thermo-Physical_Properties.cu $(TRPHPP) $(SPCHPP) $(UTLHPP)
	@mkdir -p $(BUILDDIR);
	@echo "Compiling Thermo-Physical_Properties...";
	$(CC) $(NFLAGS) $(INC) $(IBOOST) -c $(SRCDIR)/Thermo-Physical_Properties.cu -o $(BUILDDIR)/Thermo-Physical_Properties.o

Core-Shell-Diffusion : $(BUILDDIR)/Core-Shell-Diffusion.o $(UTLOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking all files for Core-Shell-Diffusion...";
	$(CC) $(LBOOST) $(BUILDDIR)/Core-Shell-Diffusion.o $(UTLOBJS) -o $(BINDIR)/Core-Shell-Diffusion

$(BUILDDIR)/Core-Shell-Diffusion.o : $(SRCDIR)/Core-Shell-Diffusion.cu $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cuh $(TRPHPP) $(SPCHPP)
	@mkdir -p $(BUILDDIR);
	@echo "Compiling Core-Shell-Diffusion...";
	$(CC) $(NFLAGS) $(INC) $(IBOOST) -c $(SRCDIR)/Core-Shell-Diffusion.cu -o $(BUILDDIR)/Core-Shell-Diffusion.o

Hetergenous_Thermal_Conductivity_Models : $(BUILDDIR)/Hetergenous_Thermal_Conductivity_Models.o $(BUILDDIR)/$(UTLDIR)/File_Generator.o
	@mkdir -p $(BINDIR);
	@echo "Linking all files for Hetergenous_Thermal_Conductivity_Models...";
	$(CC) $(LBOOST) $(BUILDDIR)/Hetergenous_Thermal_Conductivity_Models.o $(BUILDDIR)/$(UTLDIR)/File_Generator.o -o $(BINDIR)/Hetergenous_Thermal_Conductivity_Models

$(BUILDDIR)/Hetergenous_Thermal_Conductivity_Models.o : $(SRCDIR)/Hetergenous_Thermal_Conductivity_Models.cu $(TRPHPP) $(SPCHPP)
	@mkdir -p $(BUILDDIR);
	@echo "Compiling Hetergenous_Thermal_Conductivity_Models...";
	$(CC) $(NFLAGS) $(INC) $(IBOOST) -c $(SRCDIR)/Hetergenous_Thermal_Conductivity_Models.cu -o $(BUILDDIR)/Hetergenous_Thermal_Conductivity_Models.o

# ----------------------------------------------------------------------------------------------------------
# Building utilities

$(BUILDDIR)/$(UTLDIR)/Keyboard_Interrupt.o : $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(SRCDIR)/$(UTLDIR)/Keyboard_Interrupt.cpp
	@mkdir -p $(BUILDDIR)/$(UTLDIR);
	@echo "Compiling Keyboard_Interrupt..."
	$(CC) $(NFLAGS) $(INC) -c $(SRCDIR)/$(UTLDIR)/Keyboard_Interrupt.cpp -o $(BUILDDIR)/$(UTLDIR)/Keyboard_Interrupt.o

$(BUILDDIR)/$(UTLDIR)/File_Generator.o : $(INCDIR)/$(UTLDIR)/File_Generator.hpp $(SRCDIR)/$(UTLDIR)/File_Generator.cpp
	@mkdir -p $(BUILDDIR)/$(UTLDIR);
	@echo "Compiling File_Generator..."
	$(CC) $(NFLAGS) $(INC) -c $(SRCDIR)/$(UTLDIR)/File_Generator.cpp -o $(BUILDDIR)/$(UTLDIR)/File_Generator.o

# ----------------------------------------------------------------------------------------------------------
# Building PDE Problem Examples

Core-Shell-Diffusion_Example : $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Core-Shell-Diffusion_Example...";
	$(CC) $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o -o $(BINDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example

$(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cu $(PDEHPP) $(TRPHPP) $(SPCHPP)
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Core-Shell-Diffusion_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cu -o $(BUILDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o

Pellet-Flame-Propagation_Example : $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation_Example...";
	$(CC) $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example

$(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cu $(PDEHPP) $(TRPHPP) $(SPCHPP)
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cu -o $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o

$(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPECIES_HPP) $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BUILDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp -I /usr/include/boost -o $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o

Pellet-Flame-Propagation-Diffusion_Example : $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(BUILDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS) -lboost_program_options -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example

# ----------------------------------------------------------------------------------------------------------
# Building Thermo-Physical Properties Examples
Enthalpy_Example : $(BUILDDIR)/$(TRPDIR)/Enthalpy_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Enthalpy_Example...";
	$(CC) $(BUILDDIR)/$(TRPDIR)/Enthalpy_Example.o -o $(BINDIR)/$(TRPDIR)/Enthalpy_Example

$(BUILDDIR)/$(TRPDIR)/Enthalpy_Example.o : $(EXMDIR)/$(TRPDIR)/Enthalpy_Example.cu $(INCDIR)/$(TRPDIR)/Enthalpy.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Enthalpy_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Enthalpy_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Enthalpy_Example.o

Thermal_Conductivity_Example : $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Thermal_Conductivity_Example...";
	$(CC) $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Example.o -o $(BINDIR)/$(TRPDIR)/Thermal_Conductivity_Example

$(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Example.o : $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Example.cu $(INCDIR)/$(TRPDIR)/Thermal_Conductivity.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Thermal_Conductivity_Example.o

Ideal_Gas_Example : $(BUILDDIR)/$(TRPDIR)/Ideal_Gas_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Ideal_Gas_Example...";
	$(CC) $(BUILDDIR)/$(TRPDIR)/Ideal_Gas_Example.o -o $(BINDIR)/$(TRPDIR)/Ideal_Gas_Example

$(BUILDDIR)/$(TRPDIR)/Ideal_Gas_Example.o : $(EXMDIR)/$(TRPDIR)/Ideal_Gas_Example.cu $(INCDIR)/$(TRPDIR)/Ideal_Gas.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Ideal_Gas_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Ideal_Gas_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Ideal_Gas_Example.o

Phase_Example : $(BUILDDIR)/$(TRPDIR)/Phase_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Phase_Example...";
	$(CC) $(BUILDDIR)/$(TRPDIR)/Phase_Example.o -o $(BINDIR)/$(TRPDIR)/Phase_Example

$(BUILDDIR)/$(TRPDIR)/Phase_Example.o : $(EXMDIR)/$(TRPDIR)/Phase_Example.cu $(INCDIR)/$(TRPDIR)/Phase.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Phase_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Phase_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Phase_Example.o

Species_Example : $(BUILDDIR)/$(TRPDIR)/Species_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Species_Example...";
	$(CC) $(BUILDDIR)/$(TRPDIR)/Species_Example.o -o $(BINDIR)/$(TRPDIR)/Species_Example

$(BUILDDIR)/$(TRPDIR)/Species_Example.o : $(EXMDIR)/$(TRPDIR)/Species_Example.cu $(INCDIR)/$(TRPDIR)/Species.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Species_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Species_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Species_Example.o

Core-Shell-Particle_Example : $(BUILDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Core-Shell-Particle_Example...";
	$(CC) $(BUILDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o -o $(BINDIR)/$(TRPDIR)/Core-Shell-Particle_Example

$(BUILDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o : $(EXMDIR)/$(TRPDIR)/Core-Shell-Particle_Example.cu $(INCDIR)/$(TRPDIR)/Core-Shell-Particle.cuh $(SPECIES_HPP)
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Core-Shell-Particle_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Core-Shell-Particle_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o

Packed-Pellet_Example : $(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Packed-Pellet_Example...";
	$(CC) $(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o -o $(BINDIR)/$(TRPDIR)/Packed-Pellet_Example

$(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cu $(INCDIR)/$(TRPDIR)/Packed-Pellet.cuh $(SPECIES_HPP) $(INCDIR)/$(TRPDIR)/Core-Shell-Particle.cuh
	@mkdir -p $(BUILDDIR)/$(TRPDIR);
	@echo "Compiling Packed-Pellet_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cu -o $(BUILDDIR)/$(TRPDIR)/Packed-Pellet_Example.o

# ----------------------------------------------------------------------------------------------------------
# Build LU Solver examples

LU_Solver_Example : $(BUILDDIR)/$(LUSDIR)/LU_Solver_Example.o
	@mkdir -p $(BINDIR)/$(LUSDIR);
	@echo "Building LU_Solver_Example...";
	$(CC) $(BUILDDIR)/$(LUSDIR)/LU_Solver_Example.o -o $(BINDIR)/$(LUSDIR)/LU_Solver_Example

$(BUILDDIR)/$(LUSDIR)/LU_Solver_Example.o : $(EXMDIR)/$(LUSDIR)/LU_Solver_Example.cu $(LUSHPP)
	@mkdir -p $(BUILDDIR)/$(LUSDIR);
	@echo "Compiling LU_Solver_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(LUSDIR)/LU_Solver_Example.cu -o $(BUILDDIR)/$(LUSDIR)/LU_Solver_Example.o

Tridiagonal_Matrix_Example : $(BUILDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o
	@mkdir -p $(BINDIR)/$(LUSDIR);
	@echo "Building Tridiagonal_Matrix_Example...";
	$(CC) $(BUILDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o -o $(BINDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example

$(BUILDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o : $(EXMDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.cu $(INCDIR)/$(LUSDIR)/Tridiagonal_Matrix.cuh
	@mkdir -p $(BUILDDIR)/$(LUSDIR);
	@echo "Compiling Tridiagonal_Matrix_Example...";
	$(CC) $(NFLAGS) $(INC) -c $(EXMDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.cu -o $(BUILDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o

# ----------------------------------------------------------------------------------------------------------------------
# Clean everything

clean:
	@echo "Cleaning...";
	$(RM) -r $(BUILDDIR) $(BINDIR)