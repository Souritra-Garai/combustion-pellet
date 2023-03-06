# Math Headers

## Data Type

`Data-Type.hpp` defines the data type, `real_t`, for all real valued physical quantities.
A large floating point data type like `long double` can improve accuracy but requires more computation time than `float` or `double`.

## Linear Expression

`Linear-Expression.hpp` defines the class `LinearExpression` for implementing linear expressions of the form
```
a_0 + a_1 * x
```
To directly evaluate the expression for a given `x`, the following member function can be used.
```
LinearExpression linear_expr(a_0, a_1);

linear_expr.evaluateExpression(x);
```
Moreover, there are provisions for performing quick arithmetic with the linear expressions, like adding two linear expressions or multiplying with a scaler, using operators like `+`, `*`, `+=` etc.

For certain scenarios, where the following expression is required
```
a_1 * x - a_0
```
one can use the same member function like this
```
- linear_expr.evaluateExpression(-x)
```
It is primarily used for representing source terms linearized in temperature
```
f(T) = f(T_0) + f'(T_0) * (T - T_0)
```

## Quadratic Expression

`Quadratic-Expression.hpp` defines the class `QuadraticExpression` for implementing quadratic expression/polynomial of the form
```
a_0 + a_1 * x + a_2 * x^2
```
Again, to directly evaluate the expression for a given `x`, the following member function can be used.
```
QuadraticExpression quad_expr(a_0, a_1, a_2);

quad_expr.evaluateExpression(x);
```
This class is primarily used for representing thermal conductivity as a function of temperature.
```
k(T) = a_0 + a_1 * T + a_2 * T^2
```
where `k` is thermal conductivity in W/m-K and `T` is temperature in K.

## Shomate Equation

`Shomate-Expression.hpp` defines the class `ShomateExpression` for implementing Shomate equation of the form
```
A + B * x + C * x^2 + D * x^3 + E / (x^2)
```
The class is used to model specific heat, `c`, of a species, as a function of temperature. To get the specific heat in SI units, i.e. J/mol.-K, at a temperature of `T` K, use the member function
```
ShomateEquation specific_heat_model(A, B, C, D, E, F);
c(T) = specific_heat_model.evaluateExpression(
	ShomateExpression::normalizeInput(T)
)
```
Similarly, to get the standard enthalpy. `h`, of the species in J/mol., at a temperature `T` K, use the member function
```
h(T) = ShomateEquation::normalizeIntegralOutput(
	specific_heat_model.evaluateExpressionIntegral(
		ShomateExpression::normalizeInput(T)
	)
)
```

## Sigmoid Function

`Sigmoid-Function.hpp` defines functions for implementing the sigmoid function of the form
```
sigma(x) = 1 / (1 + exp(-x))
```
The function is implemented using the `tanh` function from cmath library.
```
tanh(x) = (exp(2x) - 1) / (exp(2x) + 1)

sigma(x) =  0.5 * (1 + tanh(0.5 * x))
```
The derivative of the sigmoid function is 
```
(d/dx) sigma(x) = sigma(x) * (1 - sigma(x))
```

The two functions defined in this header are
```
real_t getSigmoid(real_t x, real_t origin, real_t scale);
// = sigma( scale * (x - origin) )

real_t getSigmoidDerivative(real_t x, real_t origin, real_t scale)
// = (d/dx) sigma(x) evaluated at x = scale * (x - origin)