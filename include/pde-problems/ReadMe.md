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

## Pellet Flame Propagation
The header `Pellet-Flame-Propagation.hpp` defines the class `PelletFlamePropagation` that represents a discretized cylindrical pellet and solves the discretized pellet enthalpy equation.

The member function `calcTransientTerm` calculates and returns the main heat source term in the discetized pellet, in form of a linear expression
```
a_0 + a_1 * (T - T0)
```
where `a_0` and `a_1` are constants, `T` is what the temperature is going to be in the next time step, `T0` is the known temperature in the current time step.

The values of `a`s are determined as follows
```
a_0 = sum(((1-g) * (DY_k / Dt)_exp + g * (DY_k / Dt)_imp) * h_k(T0))
	= (1-g) * sum(D(Y_k * h_k(T0) / Dt)) + g * sum(D(Y_k * h_k(T0)) / Dt)
	= (1-g) * (h_p1 - h_p2) / Dt + g * (h_p3 - h_p1) / Dt

a_1 = g * sum( (DY_k / DT) * h_k(T0) ) / Dt + c_p / Dt
	= g * sum(D(Y_k * h_k(T0)) / DT) / Dt + c_p / Dt
	= g * (h_p4 - h_p3) / (DT * Dt) + c_p / Dt

h_p1 = enthalpy of particle at previous time step composition and temperature T0
h_p2 = enthalpy of particle at current time step composition and temperature T0
h_p3 = enthalpy of particle at evolved time step (constant T) and temperature T0
h_p4 = enthalpy of particle at evolved time step (T + DT) and temperature T0
c_p = specific heat capacity of particle at temperature T0
```

The member function `setUpEquations` sets up the tridiagonal matrix equation representing the discretized pellet enthalpy equation
```
e * T[i+1] + f * T[i] + g * T[i-1] = b

e = - kappa * lambda[i+1] / Dx^2
g = - kappa * lambda[i-1] / Dx^2
f = phi_P * rho_P * a_1 + 
	phi_F * rho_F * c_F / Dt -
	b_1 - (e + g)
b = phi_P * rho_P * (a_1 * T[i] - a_0) +
	phi_F * rho_F * c_F * T[i] / Dt +
	(1-kappa) * (lambda[i+1] * T[i+1] - (lambda[i+1] + lambda[i]) * T[i] + lambda[i] * T[i-1]) / Dx^2 +
	b_0 - b_1 * T[i]
```
