# PDE Problems

This group of libraries defines classes to represent the discretized pellet and core-shell particles; and solve the discretized governing PDE of the respective object.

## Core-Shell Particle
The header `Core-Shell-Diffusion.hpp` defines the class `CoreShellDiffusion` that represents a discretized core-shell spherical particle and solves the discretized diffusion equation for two species A and B. Species A represents the elemental species initially present in the core and species B represents the same initially present in the shell.

The member function `updateMassFractions` sums the chemical concentrations of different species in the particle
```
n = number of grid points in radial coordinate
n-1 = number of intervals

m = molar_mass * integral( C(r) * 4 * pi * r^2, r )
  = 4 * pi * molar_mass *
    sum(
		0.5 * (C(r+Dr) * (r+Dr)^2   +   C(r) * r^2) * Dr,
	    r = 0 to (n-2)*Dr (= n-1 iterations)
	)

  = some_constant * molar_mass * (
	0.5 * C((n-1)*Dr) * ((n-1)*Dr)^2 +
	sum(
		C(r) * r^2,
		r = 1*Dr to (n-2)*Dr (= n-2 iterations)
	)
  )
```

The member function `setUpEquations` sets up the tridiagonal matrix equation to representing the discretized diffusion equation
```
g = 0.5 * D * (1/Dr)^2
e = g * ((r+Dr) / r)^2

s = e + g
f1 = 1/Dt + s
f2 = 1/Dt - s

b = e * x[i+1] + f2 * x[i] + g * x[i-1]

// Final equation for LUSolver
- e * x[i+1] + f1 * x[i] - g * x[i-1] = b
```
