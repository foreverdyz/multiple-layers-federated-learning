#new_method.jl

function new_estimate(data1::AbstractArray, data2::AbstractArray, k1::Int, k2::Int, lambdac::Real,lambdau::Real)
    let
        #initialize two vectors to store final models
        w1 = zeros(k1)
        w2 = zeros(k2)

        s1 = sum(data1)
        s2 = sum(data2)
        #global model
        w = (s1+s2)/(k1+k2)
        #class models
        c1 = (s1/k1+lambdac*w)/(1+lambdac)
        c2 = (s2/k2+lambdac*w)/(1+lambdac)
        #final models
        w1 = (data1.+(lambdau*c1))/(1+lambdau)
        w2 = (data2.+(lambdau*c2))/(1+lambdau)

        return w1,w2
    end
end
