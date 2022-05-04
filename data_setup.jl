#data_setup.jl
using Random, Distributions
#c, center of distribution; k, number of devices; n, size of samples
function generate_data(c::Real, k::Int, n::Int, tau::Real)
    let
        #generate underlying model for each device
        v = rand(Normal(c, 0.1), k)

        #initialize data
        data = zeros(k)

        for i=1:k
            #generate data for device i
            data[i] = rand(Normal(v[i], tau/(n^0.5)), 1)[1]
        end

        return data,v
    end
end
