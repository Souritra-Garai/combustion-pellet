# Compiler
CC := g++

# Directories
SRCDIR := src
INCDIR := include
BINDIR := bin
EXMDIR := examples
BLDDIR := build

LUSDIR := lusolver
QRSDIR := qrsolver
TRPDIR := thermo-physical-properties
PDEDIR := pde-problems
SPCDIR := species
UTLDIR := utilities

# Finding all header, source and object files
SRCEXT := cpp
HPPEXT := hpp

LUSHPPS := $(shell find $(INCDIR)/$(LUSDIR) -type f -name *.$(HPPEXT))
LUSSRCS := $(shell find $(SRCDIR)/$(LUSDIR) -type f -name *.$(SRCEXT))
LUSOBJS := $(patsubst $(SRCDIR)/%, $(BLDDIR)/%, $(LUSSRCS:.$(SRCEXT)=.o))

QRSHPPS := $(shell find $(INCDIR)/$(QRSDIR) -type f -name *.$(HPPEXT))
QRSSRCS := $(shell find $(SRCDIR)/$(QRSDIR) -type f -name *.$(SRCEXT))
QRSOBJS := $(patsubst $(SRCDIR)/%, $(BLDDIR)/%, $(QRSSRCS:.$(SRCEXT)=.o))

TRPHPPS := $(shell find $(INCDIR)/$(TRPDIR) -type f -name *.$(HPPEXT))
TRPSRCS := $(shell find $(SRCDIR)/$(TRPDIR) -type f -name *.$(SRCEXT))
TRPOBJS := $(patsubst $(SRCDIR)/%, $(BLDDIR)/%, $(TRPSRCS:.$(SRCEXT)=.o))

PDEHPPS := $(shell find $(INCDIR)/$(PDEDIR) -type f -name *.$(HPPEXT))
PDESRCS := $(shell find $(SRCDIR)/$(PDEDIR) -type f -name *.$(SRCEXT))
PDEOBJS := $(patsubst $(SRCDIR)/%, $(BLDDIR)/%, $(PDESRCS:.$(SRCEXT)=.o))

UTLHPPS := $(shell find $(INCDIR)/$(UTLDIR) -type f -name *.$(HPPEXT))
UTLSRCS := $(shell find $(SRCDIR)/$(UTLDIR) -type f -name *.$(SRCEXT))
UTLOBJS := $(patsubst $(SRCDIR)/%, $(BLDDIR)/%, $(UTLSRCS:.$(SRCEXT)=.o))

SPCHPPS := $(shell find $(INCDIR)/$(SPCDIR) -type f -name *.$(HPPEXT))

# Flags required for compiler
CFLAGS := -fopenmp -O2
LIB := -lm
INC := -I $(INCDIR)

$(TARGET) : $(OBJECTS)
	@mkdir -p $(BINDIR);
	@echo "Linking all...";
	$(CC) $(OBJECTS) $(CFLAGS) $(LIB) -o $(TARGET)

main : $(BLDDIR)/main.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking all files for main...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/main.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS) -lboost_program_options -o $(BINDIR)/Pellet-Flame-Propagation

$(BLDDIR)/main.o : $(SRCDIR)/main.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPCHPPS) $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR);
	@echo "Compiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -I /usr/include/boost -o $(BLDDIR)/main.o

Thermo-Physical_Properties : $(BLDDIR)/Thermo-Physical_Properties.o $(UTLOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking all files for Thermo-Physical_Properties...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/Thermo-Physical_Properties.o $(UTLOBJS) -o $(BINDIR)/Thermo-Physical_Properties

$(BLDDIR)/Thermo-Physical_Properties.o : $(SRCDIR)/Thermo-Physical_Properties.cpp $(TRPHPPS) $(SPCHPPS) $(UTLHPPS)
	@mkdir -p $(BLDDIR);
	@echo "Compiling Thermo-Physical_Properties...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Thermo-Physical_Properties.cpp -o $(BLDDIR)/Thermo-Physical_Properties.o

Core-Shell-Diffusion : $(BLDDIR)/Core-Shell-Diffusion.o $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking Core-Shell-Diffusion...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/Core-Shell-Diffusion.o $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(LUSOBJS) $(UTLOBJS) -o $(BINDIR)/Core-Shell-Diffusion

$(BLDDIR)/Core-Shell-Diffusion.o : $(SRCDIR)/Core-Shell-Diffusion.cpp $(PDEHPPS) $(TRPHPPS) $(SPCHPPS) $(UTLHPPS)
	@mkdir -p $(BLDDIR);
	@echo "Compiling Core-Shell-Diffusion...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Core-Shell-Diffusion.cpp -o $(BLDDIR)/Core-Shell-Diffusion.o

Heterogenous_Thermal_Conductivity_Models : $(BLDDIR)/Heterogenous_Thermal_Conductivity_Models.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(UTLOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking Heterogenous_Thermal_Conductivity_Models...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/Heterogenous_Thermal_Conductivity_Models.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(UTLOBJS) -o $(BINDIR)/Heterogenous_Thermal_Conductivity_Models

$(BLDDIR)/Heterogenous_Thermal_Conductivity_Models.o : $(SRCDIR)/Heterogenous_Thermal_Conductivity_Models.cpp $(TRPHPPS) $(SPCHPPS) $(UTLHPPS)
	@mkdir -p $(BLDDIR);
	@echo "Compiling Heterogenous_Thermal_Conductivity_Models...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Heterogenous_Thermal_Conductivity_Models.cpp -o $(BLDDIR)/Heterogenous_Thermal_Conductivity_Models.o

Pellet-Flame-Propagation : $(BLDDIR)/Pellet-Flame-Propagation.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR);
	@echo "Linking Pellet-Flame-Propagation...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/Pellet-Flame-Propagation.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS) -o $(BINDIR)/Pellet-Flame-Propagation

$(BLDDIR)/Pellet-Flame-Propagation.o : $(SRCDIR)/Pellet-Flame-Propagation.cpp $(PDEHPPS) $(TRPHPPS) $(LUSHPPS) $(SPCHPPS) $(UTLHPPS)
	@mkdir -p $(BLDDIR);
	@echo "Compiling Pellet-Flame-Propagation...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Pellet-Flame-Propagation.cpp -o $(BLDDIR)/Pellet-Flame-Propagation.o

# ----------------------------------------------------------------------------------------------------------
# Building thermo-physical properties source files

$(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o : $(SRCDIR)/$(TRPDIR)/Core-Shell-Particle.cpp $(TRPHPPS)
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Core-Shell-Particle...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Core-Shell-Particle.cpp -o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o

$(BLDDIR)/$(TRPDIR)/Packed-Pellet.o : $(SRCDIR)/$(TRPDIR)/Packed-Pellet.cpp $(TRPHPPS)
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Packed-Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Packed-Pellet.cpp -o $(BLDDIR)/$(TRPDIR)/Packed-Pellet.o

$(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o : $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SRCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.cpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.cpp -o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o

# ----------------------------------------------------------------------------------------------------------
# Building pde problems source files

$(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o : $(SRCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cpp $(PDEHPPS) $(LUSHPPS) $(TRPHPPS)
	@mkdir -p $(BLDDIR)/$(PDEDIR);
	@echo "Compiling Core-Shell-Diffusion...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(PDEDIR)/Core-Shell-Diffusion.cpp -o $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o

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

Core-Shell-Diffusion_Example : $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(LUSOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(LUSOBJS) -o $(BINDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example

$(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cpp $(PDEHPPS) $(TRPHPPS) $(SPCHPPS)
	@mkdir -p $(BLDDIR)/$(PDEDIR);
	@echo "Compiling Core-Shell-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.cpp -o $(BLDDIR)/$(PDEDIR)/Core-Shell-Diffusion_Example.o

Pellet-Flame-Propagation_Example : $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS) -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example

$(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cpp $(PDEHPPS) $(TRPHPPS) $(LUSHPPS) $(SPCHPPS) $(UTLHPPS)
	@mkdir -p $(BLDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.cpp -o $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation_Example.o

$(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o : $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp $(INCDIR)/$(PDEDIR)/Core-Shell-Diffusion.hpp $(INCDIR)/$(PDEDIR)/Pellet-Flame-Propagation.hpp $(INCDIR)/$(TRPDIR)/Arrhenius_Diffusivity_Model.hpp $(SPCHPPS) $(INCDIR)/$(UTLDIR)/Keyboard_Interrupt.hpp $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR)/$(PDEDIR);
	@echo "Compiling Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.cpp -I /usr/include/boost -o $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o

Pellet-Flame-Propagation-Diffusion_Example : $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(PDEDIR);
	@echo "Linking Pellet-Flame-Propagation-Diffusion_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example.o $(PDEOBJS) $(TRPOBJS) $(LUSOBJS) $(UTLOBJS) -lboost_program_options -o $(BINDIR)/$(PDEDIR)/Pellet-Flame-Propagation-Diffusion_Example

# ----------------------------------------------------------------------------------------------------------
# Building Thermo-Physical Properties Examples
Thermal_Conductivity_Pellet_Example : $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Thermal_Conductivity_Pellet_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(UTLOBJS) -o $(BINDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example

$(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.cpp $(INCDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.hpp $(SPCHPPS) $(INCDIR)/$(TRPDIR)/Core-Shell-Particle.hpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Thermal_Conductivity_Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.cpp -o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet_Example.o

Adiabatic_Combustion_Temperature : $(BLDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o $(UTLOBJS)
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Adiabatic_Combustion_Temperature...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o $(UTLOBJS) -o $(BINDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature

$(BLDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o : $(EXMDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.cpp $(SPCHPPS) $(INCDIR)/$(UTLDIR)/File_Generator.hpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Adiabatic_Combustion_Temperature...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.cpp -o $(BLDDIR)/$(TRPDIR)/Adiabatic_Combustion_Temperature.o

Phase_Example : $(BLDDIR)/$(TRPDIR)/Phase_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Phase_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Phase_Example.o -o $(BINDIR)/$(TRPDIR)/Phase_Example

$(BLDDIR)/$(TRPDIR)/Phase_Example.o : $(EXMDIR)/$(TRPDIR)/Phase_Example.cpp $(INCDIR)/$(TRPDIR)/Phase.hpp
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Phase_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Phase_Example.cpp -o $(BLDDIR)/$(TRPDIR)/Phase_Example.o

Condensed_Species_Example : $(BLDDIR)/$(TRPDIR)/Condensed_Species_Example.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Condensed_Species_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Condensed_Species_Example.o -o $(BINDIR)/$(TRPDIR)/Condensed_Species_Example

$(BLDDIR)/$(TRPDIR)/Condensed_Species_Example.o : $(EXMDIR)/$(TRPDIR)/Condensed_Species_Example.cpp $(SPCHPPS) $(TRPHPPS)
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Condensed_Species_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Condensed_Species_Example.cpp -o $(BLDDIR)/$(TRPDIR)/Condensed_Species_Example.o

Core-Shell-Particle_Example : $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Core-Shell-Particle_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o -o $(BINDIR)/$(TRPDIR)/Core-Shell-Particle_Example

$(BLDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o : $(EXMDIR)/$(TRPDIR)/Core-Shell-Particle_Example.cpp $(TRPHPPS) $(SPCHPPS)
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Core-Shell-Particle_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Core-Shell-Particle_Example.cpp -o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle_Example.o

Packed-Pellet_Example : $(BLDDIR)/$(TRPDIR)/Packed-Pellet_Example.o $(BLDDIR)/$(TRPDIR)/Packed-Pellet.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o
	@mkdir -p $(BINDIR)/$(TRPDIR);
	@echo "Linking Packed-Pellet_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(TRPDIR)/Packed-Pellet_Example.o $(BLDDIR)/$(TRPDIR)/Packed-Pellet.o $(BLDDIR)/$(TRPDIR)/Core-Shell-Particle.o $(BLDDIR)/$(TRPDIR)/Thermal_Conductivity_Pellet.o -o $(BINDIR)/$(TRPDIR)/Packed-Pellet_Example

$(BLDDIR)/$(TRPDIR)/Packed-Pellet_Example.o : $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cpp $(TRPHPPS)
	@mkdir -p $(BLDDIR)/$(TRPDIR);
	@echo "Compiling Packed-Pellet_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(TRPDIR)/Packed-Pellet_Example.cpp -o $(BLDDIR)/$(TRPDIR)/Packed-Pellet_Example.o

# ----------------------------------------------------------------------------------------------------------
# Building LU Decomposition solver
$(BLDDIR)/$(LUSDIR)/LU_Solver.o : $(SRCDIR)/$(LUSDIR)/LU_Solver.cpp $(INCDIR)/$(LUSDIR)/LU_Solver.hpp
	@mkdir -p $(BLDDIR)/$(LUSDIR);
	@echo "Compiling LU_Solver...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(LUSDIR)/LU_Solver.cpp -o $(BLDDIR)/$(LUSDIR)/LU_Solver.o

# Building Tridiagonal Matrix
$(BLDDIR)/$(LUSDIR)/Tridiagonal_Matrix.o : $(SRCDIR)/$(LUSDIR)/Tridiagonal_Matrix.cpp $(INCDIR)/$(LUSDIR)/Tridiagonal_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(LUSDIR);
	@echo "Compiling Tridiagonal_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(LUSDIR)/Tridiagonal_Matrix.cpp -o $(BLDDIR)/$(LUSDIR)/Tridiagonal_Matrix.o

# Builds the example for LU Solver
LU_Solver_Example : $(BLDDIR)/$(LUSDIR)/LU_Solver_Example.o $(LUSOBJS)
	@mkdir -p $(BINDIR)/$(LUSDIR);
	@echo "Building LU_Solver_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(LUSDIR)/LU_Solver_Example.o $(LUSOBJS) -o $(BINDIR)/$(LUSDIR)/LU_Solver_Example

$(BLDDIR)/$(LUSDIR)/LU_Solver_Example.o : $(EXMDIR)/$(LUSDIR)/LU_Solver_Example.cpp $(INCDIR)/$(LUSDIR)/LU_Solver.hpp
	@mkdir -p $(BLDDIR)/$(LUSDIR);
	@echo "Compiling LU_Solver_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(LUSDIR)/LU_Solver_Example.cpp -o $(BLDDIR)/$(LUSDIR)/LU_Solver_Example.o

# Builds the example for Tridiagonal Matrix
Tridiagonal_Matrix_Example : $(BLDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o $(BLDDIR)/$(LUSDIR)/Tridiagonal_Matrix.o
	@mkdir -p $(BINDIR)/$(LUSDIR);
	@echo "Building Tridiagonal_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o $(BLDDIR)/$(LUSDIR)/Tridiagonal_Matrix.o -o $(BINDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example

$(BLDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o : $(EXMDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.cpp $(INCDIR)/$(LUSDIR)/Tridiagonal_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(LUSDIR);
	@echo "Compiling Tridiagonal_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.cpp -o $(BLDDIR)/$(LUSDIR)/Tridiagonal_Matrix_Example.o

# ----------------------------------------------------------------------------------------------------------
# Building QR Factorization solver

$(BLDDIR)/$(QRSDIR)/R_Matrix.o : $(SRCDIR)/$(QRSDIR)/R_Matrix.cpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(QRSDIR);
	@echo "Compiling R_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/R_Matrix.cpp -o $(BLDDIR)/$(QRSDIR)/R_Matrix.o

$(BLDDIR)/$(QRSDIR)/Q_Matrix.o : $(SRCDIR)/$(QRSDIR)/Q_Matrix.cpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(QRSDIR);
	@echo "Compiling Q_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/Q_Matrix.cpp -o $(BLDDIR)/$(QRSDIR)/Q_Matrix.o

$(BLDDIR)/$(QRSDIR)/G_Matrix.o : $(SRCDIR)/$(QRSDIR)/G_Matrix.cpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(QRSDIR);
	@echo "Compiling G_Matrix...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/G_Matrix.cpp -o $(BLDDIR)/$(QRSDIR)/G_Matrix.o 

$(BLDDIR)/$(QRSDIR)/QR_Solver.o : $(SRCDIR)/$(QRSDIR)/QR_Solver.cpp $(INCDIR)/$(QRSDIR)/QR_Solver.hpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(QRSDIR);
	@echo "Compiling QR_Solver...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/$(QRSDIR)/QR_Solver.cpp -o $(BLDDIR)/$(QRSDIR)/QR_Solver.o

# Builds the example for R Matrix
R_Matrix_Example : $(BLDDIR)/$(QRSDIR)/R_Matrix_Example.o $(BLDDIR)/$(QRSDIR)/R_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "Building R_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(QRSDIR)/R_Matrix_Example.o $(BLDDIR)/$(QRSDIR)/R_Matrix.o -o $(BINDIR)/$(QRSDIR)/R_Matrix_Example

$(BLDDIR)/$(QRSDIR)/R_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/R_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/R_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(QRSDIR);
	@echo "Compiling R_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/R_Matrix_Example.cpp -o $(BLDDIR)/$(QRSDIR)/R_Matrix_Example.o

# Builds the example for Q Matrix
Q_Matrix_Example : $(BLDDIR)/$(QRSDIR)/Q_Matrix_Example.o $(BLDDIR)/$(QRSDIR)/Q_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "Building Q_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(QRSDIR)/Q_Matrix_Example.o $(BLDDIR)/$(QRSDIR)/Q_Matrix.o -o $(BINDIR)/$(QRSDIR)/Q_Matrix_Example

$(BLDDIR)/$(QRSDIR)/Q_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/Q_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/Q_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(QRSDIR);
	@echo "Compiling Q_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/Q_Matrix_Example.cpp -o $(BLDDIR)/$(QRSDIR)/Q_Matrix_Example.o

# Builds the example for Givens Rotation matrix
G_Matrix_Example : $(BLDDIR)/$(QRSDIR)/G_Matrix_Example.o $(BLDDIR)/$(QRSDIR)/G_Matrix.o $(BLDDIR)/$(QRSDIR)/Q_Matrix.o $(BLDDIR)/$(QRSDIR)/R_Matrix.o
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "Building G_Matrix_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(QRSDIR)/G_Matrix_Example.o $(BLDDIR)/$(QRSDIR)/G_Matrix.o $(BLDDIR)/$(QRSDIR)/Q_Matrix.o $(BLDDIR)/$(QRSDIR)/R_Matrix.o -o $(BINDIR)/$(QRSDIR)/G_Matrix_Example

$(BLDDIR)/$(QRSDIR)/G_Matrix_Example.o : $(EXMDIR)/$(QRSDIR)/G_Matrix_Example.cpp $(INCDIR)/$(QRSDIR)/G_Matrix.hpp
	@mkdir -p $(BLDDIR)/$(QRSDIR);
	@echo "Compiling G_Matrix_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/G_Matrix_Example.cpp -o $(BLDDIR)/$(QRSDIR)/G_Matrix_Example.o

# Builds the example for QR Solver
QR_Solver_Example : $(BLDDIR)/$(QRSDIR)/QR_Solver_Example.o $(QRSOBJS)
	@mkdir -p $(BINDIR)/$(QRSDIR);
	@echo "Building QR_Solver_Example...";
	$(CC) $(CFLAGS) $(LIB) $(BLDDIR)/$(QRSDIR)/QR_Solver_Example.o $(QRSOBJS) -o $(BINDIR)/$(QRSDIR)/QR_Solver_Example

$(BLDDIR)/$(QRSDIR)/QR_Solver_Example.o : $(EXMDIR)/$(QRSDIR)/QR_Solver_Example.cpp $(INCDIR)/$(QRSDIR)/QR_Solver.hpp
	@mkdir -p $(BLDDIR)/$(QRSDIR);
	@echo "Compiling QR_Solver_Example...";
	$(CC) $(CFLAGS) $(INC) -c $(EXMDIR)/$(QRSDIR)/QR_Solver_Example.cpp -o $(BLDDIR)/$(QRSDIR)/QR_Solver_Example.o

# ----------------------------------------------------------------------------------------------------------------------
# Clean everything

clean:
	@echo "Cleaning...";
	$(RM) -r $(BLDDIR) $(BINDIR)