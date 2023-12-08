using Zygote
using ChainRulesCore
using LinearAlgebra

struct ComplexField{T} 
    z::T
end

function Base.adjoint(a::ComplexField{T}) where T
    return  ComplexField{T}(a.z')
end

function Base.:*(a::ComplexField,b::ComplexField) 
    return ComplexField(a.z*b.z)
end

function Base.:*(a::T,b::ComplexField) where T<:Number
    return ComplexField(a*b.z)
end

function Base.:*(b::ComplexField,a::T) where T<:Number
    return ComplexField(a*b.z)
end


function Base.:/(a::ComplexField,b::ComplexField) 
    return ComplexField(a.z/b.z)
end

function Base.:/(a::ComplexField,b::T)  where T<:Number
    return ComplexField(a.z/b)
end

function Base.:/(a::T,b::ComplexField)  where T<:Number
    return ComplexField(a/b.z)
end

function Base.:+(a::ComplexField,b::ComplexField) 
    return ComplexField(a.z+b.z)
end

function Base.:-(a::ComplexField,b::ComplexField) 
    return ComplexField(a.z-b.z)
end


function Base.:*(b::ComplexField,a::T) where T<:Number
    return ComplexField(a*b.z)
end

function Base.:^(a::ComplexField,n::Integer) 
    return ComplexField(a.z^n)
end

function Base.display(a::ComplexField)
    display(a.z)
end

function LinearAlgebra.tr(A::ComplexField)
    return ComplexField(tr(A.z))
end

function Base.real(a::ComplexField{T}) where T <: Number
    ar = (a.z+a.z')/2
    #println("real in")
    return real(ar)
end

function Base.one(a::ComplexField)
    return ComplexField(one(a.z))
end


function ChainRulesCore.rrule(::typeof(real),a::T1) where T1 <: ComplexField
    y = real(a)
    function pullback(ybar)
        sbar = NoTangent()
        #println("real ", ybar)
        fbar = ComplexField(ybar/2+0im)
        return sbar,fbar
    end
    return y, pullback
end
function ChainRulesCore.rrule(::typeof(*),a::T1,b::T1)  where T1 <: ComplexField
    y = a * b
    function pullback(ybar)
        sbar = NoTangent()
        fabar = ybar*b
        fbbar = a*ybar
        return sbar,fabar,fbbar
    end
    return y, pullback
end


function ChainRulesCore.rrule(::typeof(*),a::ComplexField{T1},b::ComplexField{T1})  where T1 <: AbstractMatrix
    y = a * b
    function pullback(ybar)
        sbar = NoTangent()
        fabar = ybar*b
        fbbar = a*ybar
        return sbar,fabar,fbbar
    end
    return y, pullback
end
function ChainRulesCore.rrule(::typeof(Base.adjoint),a::T1) where T1 <: ComplexField
    y =a'
    function pullback(ybar)
        sbar = NoTangent()
        
        fbar = ybar' #ZeroTangent()
       #println("fybar ad " , ybar)
        return sbar,fbar
    end
    return y, pullback
end
function ChainRulesCore.rrule(::typeof(+),a::T1,b::T1)  where T1 <: ComplexField
    y = a + b
    function pullback(ybar)
        sbar = NoTangent()
        fabar = ybar#*b
        fbbar = ybar
        return sbar,fabar,fbbar
    end
    return y, pullback
end

function ChainRulesCore.rrule(::typeof(-),a::T1,b::T1)  where T1 <: ComplexField
    y = a - b
    function pullback(ybar)
        sbar = NoTangent()
        fabar = ybar#*b
        fbbar = -1*ybar
        return sbar,fabar,fbbar
    end
    return y, pullback
end

function ChainRulesCore.rrule(::typeof(ComplexField),a) 
    y =ComplexField(a)
    function pullback(ybar)
        sbar = NoTangent()
        
        fbar = ybar.z #ZeroTangent()
       #println("fybar ad " , ybar)
        return sbar,fbar
    end
    return y, pullback
end

