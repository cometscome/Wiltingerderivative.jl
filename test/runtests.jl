using Wiltingerderivative
using Test
using Zygote

function test1()
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

    #=
    z = ComplexField(2+3im)
    fz(z) = z^4+2*z*z' + z
    println(f(z))
    gnu = numerical_Wiltingerderivative(fz,z)
    println("Numerical grad: ", gnu)

    g = gradient(fz,z)[1]
    println("Autograd: ", g)

    =#
end

@testset "Wiltingerderivative.jl" begin
    test1()
    
    # Write your tests here.
end
