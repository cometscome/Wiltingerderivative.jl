# Wiltingerderivative.jl

[![Build Status](https://github.com/cometscome/Wiltingerderivative.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cometscome/Wiltingerderivative.jl/actions/workflows/CI.yml?query=branch%3Amain)


# Definition
We consider the complex plane (in a sense of expressing a complex number $z = x+iy$ real numbers $x$ and $y$). The Wirtinger derivatives are defined as the following linear partial differential operators of first order:

```math
\begin{align}
\frac{\partial f}{\partial z} = \frac{1}{2} \left(\frac{\partial f}{\partial x}  - i \frac{\partial f}{\partial y}  \right) \\
\frac{\partial f}{\partial \bar{z}} = \frac{1}{2} \left(\frac{\partial f}{\partial x}  + i \frac{\partial f}{\partial y}  \right) 
\end{align}
```

# Basic properties
Linearity:
```math
\begin{align}
\frac{\partial}{\partial z} (\alpha f + \beta g) = \alpha \frac{\partial f}{\partial z} +  \beta \frac{\partial g}{\partial z} \\
\frac{\partial}{\partial \bar{z}} (\alpha f + \beta g) = \alpha \frac{\partial f}{\partial \bar{z}} +  \beta \frac{\partial g}{\partial \bar{z}} 
\end{align}
```
Product rule：
```math
\begin{align}
\frac{\partial}{\partial z} (f g) = \frac{\partial f}{\partial z} g + f \frac{\partial g}{\partial z} \\
\frac{\partial}{\partial \bar{z}} (f g) = \frac{\partial f}{\partial \bar{z}} g + f \frac{\partial g}{\partial \bar{z}}
\end{align}
```
Chain rule:
```math
\begin{align}
\frac{\partial}{\partial z} f(g(z)) = \frac{\partial f}{\partial g} \frac{\partial g}{\partial z} + \frac{\partial f}{\partial \bar{g}} \frac{\partial \bar{g}}{\partial z} \\
\frac{\partial}{\partial \bar{z}} f(g(z)) = \frac{\partial f}{\partial g} \frac{\partial g}{\partial \bar{z}} + \frac{\partial f}{\partial \bar{g}} \frac{\partial \bar{g}}{\partial \bar{z}}
\end{align}
```

# Examples
```math
f(z) = z^4 + 2 z \bar{z} + z
```

```math
\begin{align}
\frac{\partial f}{\partial z} = 4 z^3 + 2 \bar{z} + 1 \\
\frac{\partial f}{\partial \bar{z}} = 2 z 
\end{align}
```

# Code examples
```julia

a = ComplexField(2+3im)
b = ComplexField(4+im)
c = ComplexField(1+2im)
println(a*b)
println(a*b')
f(x) = real(x)
gnu = numerical_Wiltingerderivative(f,a)
println("Numerical grad: ", gnu)

g = gradient(f,a)[1]
println("Autograd: ", g)
@test abs(gnu-g.z) < 1e-4

f2(x) = real(x*x'+x+x'+x*x*x')
gnu = numerical_Wiltingerderivative(f2,a)
println("Numerical grad: ", gnu)
g = gradient(f2,a)[1]
println("Autograd: ", g)
@test abs(gnu-g.z) < 1e-4
```

