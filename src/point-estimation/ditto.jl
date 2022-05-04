#ditto.jl

function ditto_estimate(data1::AbstractArray, data2::AbstractArray, k1::Int, k2::Int, lambda::Real)
    let
        #initialize two vectors to store final models
        w1 = zeros(k1)
        w2 = zeros(k2)
        #global model
        w = (sum(data1)+sum(data2))/(k1+k2)
        #final models
        w1 = (data1.+(lambda*w))/(1+lambda)
        w2 = (data2.+(lambda*w))/(1+lambda)

        return w1,w2
    end
end
