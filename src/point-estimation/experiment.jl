#experiment.jl
using Plots, BangBang, Random, Distributions
import PyPlot
pyplot()

#=
#balance data, error bound for different lambda
n=10
k1=10
k2=10
Random.seed!(1234)
let
    re_ditto=zeros(31)
    re_new=zeros(31)
    for t=1:100
        data1,v1 = generate_data(0, k1, n, 1);
        data2,v2 = generate_data(0.5, k2, n, 1);
        #ditto with different lambda
        ditto_re = zeros(31)
        for i in 1:31
            lambda = 0.4*(i-1)
            u1,u2 = ditto_estimate(data1, data2, k1, k2, lambda)
            ditto_re[i] = (sum((v1-u1).^2)+sum((v2-u2).^2)).^0.5/(k1+k2)
        end
        re_ditto = re_ditto + ditto_re
        #new method with different lambda
        new_re = zeros(31)
        for i in 1:31
            lambda = 0.4*(i-1)
            u1,u2 = new_estimate(data1, data2, k1, k2, 0.1*lambda, lambda)
            new_re[i] = (sum((v1-u1).^2)+sum((v2-u2).^2)).^0.5/(k1+k2)
        end
        re_new = re_new + new_re
    end
    plot(0.4*((1:31).-1), re_ditto/100, xticks=1:1:12, label="Ditto, \$λ\$", xlabel="λ", ylabel="test error")
    plot!(0.4*((1:31).-1), re_new/100, label="MLFL, \$(0.1λ,λ)\$", linestyle=:dashdot)
    savefig("two distribution.png")
end
=#
#unbalance data, error for different class
n=10
k1=10
#k2=10
let
    re_ditto1 = zeros(10)
    re_ditto2 = zeros(10)
    re_ditto3 = zeros(10)
    re_ditto4 = zeros(10)
    re_new1 = zeros(10)
    re_new2 = zeros(10)
    re_new3 = zeros(10)
    re_new4 = zeros(10)
    re_ditto5 = zeros(10)
    re_ditto6 = zeros(10)
    re_ditto7 = zeros(10)
    re_ditto8 = zeros(10)
    re_new5 = zeros(10)
    re_new6 = zeros(10)
    re_new7 = zeros(10)
    re_new8 = zeros(10)
    Random.seed!(1234)
    for i=1:10
        k2 = k1+(i-1)*10
        s1_ditto = 0
        s2_ditto = 0
        s1_new = 0
        s2_new = 0
        s3_ditto = 0
        s4_ditto = 0
        s3_new = 0
        s4_new = 0
        s5_ditto = 0
        s6_ditto = 0
        s5_new = 0
        s6_new = 0
        s7_ditto = 0
        s8_ditto = 0
        s7_new = 0
        s8_new = 0
        for t=1:100
            data1,v1 = generate_data(0, k1, n, 1);
            data2,v2 = generate_data(1, k2, n, 1);

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 2)
            s1_ditto += (sum((v1-u1).^2)^0.5)/k1
            s2_ditto += (sum((v2-u2).^2)^0.5)/k2

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 1)
            s3_ditto += (sum((v1-u1).^2)^0.5)/k1
            s4_ditto += (sum((v2-u2).^2)^0.5)/k2

            u1,u2 = new_estimate(data1, data2, k1, k2, 1, 4)
            s1_new += (sum((v1-u1).^2)^0.5)/k1
            s2_new += (sum((v2-u2).^2)^0.5)/k2

            u1,u2 = new_estimate(data1, data2, k1, k2, 0.5, 4)
            s3_new += (sum((v1-u1).^2)^0.5)/k1
            s4_new += (sum((v2-u2).^2)^0.5)/k2

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 0.5)
            s5_ditto += (sum((v1-u1).^2)^0.5)/k1
            s6_ditto += (sum((v2-u2).^2)^0.5)/k2

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 0.25)
            s7_ditto += (sum((v1-u1).^2)^0.5)/k1
            s8_ditto += (sum((v2-u2).^2)^0.5)/k2

            u1,u2 = new_estimate(data1, data2, k1, k2, 0.25, 4)
            s5_new += (sum((v1-u1).^2)^0.5)/k1
            s6_new += (sum((v2-u2).^2)^0.5)/k2

            u1,u2 = new_estimate(data1, data2, k1, k2, 0.1, 4)
            s7_new += (sum((v1-u1).^2)^0.5)/k1
            s8_new += (sum((v2-u2).^2)^0.5)/k2
        end
        re_ditto1[i] = s1_ditto/s2_ditto
        re_new1[i] = s1_new/s2_new
        re_ditto2[i] = s3_ditto/s4_ditto
        re_new2[i] = s3_new/s4_newœ
        re_ditto3[i] = s5_ditto/s6_ditto
        re_new3[i] = s5_new/s6_new
        re_ditto4[i] = s7_ditto/s8_ditto
        re_new4[i] = s7_new/s8_new
    end
    plot(1:10, re_ditto1, label="Ditto, \$λ=2\$", lw=2, xlabel="\$K_2/K_1\$", ylabel="\$error_1/error_2\$")
    plot!(1:10, re_ditto2, label="Ditto, \$λ=1\$", lw=2)
    plot!(1:10, re_ditto3, label="Ditto, \$λ=0.5\$", lw=2)
    plot!(1:10, re_ditto4, label="Ditto, \$λ=0.25\$", lw=2)
    plot!(1:10, re_new1, label="MLFL, \$λ=(1,4)\$", linestyle=:dashdot)
    plot!(1:10, re_new2, label="MLFL, \$λ=(0.5,4)\$", linestyle=:dashdot)
    plot!(1:10, re_new3, label="MLFL, \$λ=(0.25,4)\$", linestyle=:dashdot)
    plot!(1:10, re_new4, label="MLFL, \$λ=(0.1,4)\$", linestyle=:dashdot)
    savefig("fairness.png")
end
#=
n=10
k1=10
let
    s1_ditto = zeros(10)
    s2_ditto = zeros(10)
    s3_ditto = zeros(10)
    s4_ditto = zeros(10)
    s1_new = zeros(10)
    s2_new = zeros(10)
    s3_new = zeros(10)
    s4_new = zeros(10)
    Random.seed!(1234)
    for i=1:10
        k2 = k1+(i-1)*10
        for t=1:100
            data1,v1 = generate_data(0, k1, n, 2);
            data2,v2 = generate_data(1, k2, n, 2);

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 2)
            s1_ditto[i] += ((sum((v1-u1).^2)+sum((v2-u2).^2))^0.5)/(k1+k2)
            u1,u2 = ditto_estimate(data1, data2, k1, k2, 1)
            s2_ditto[i] += ((sum((v1-u1).^2)+sum((v2-u2).^2))^0.5)/(k1+k2)
            u1,u2 = ditto_estimate(data1, data2, k1, k2, 0.5)
            s3_ditto[i] += ((sum((v1-u1).^2)+sum((v2-u2).^2))^0.5)/(k1+k2)
            u1,u2 = ditto_estimate(data1, data2, k1, k2, 0.25)
            s4_ditto[i] += ((sum((v1-u1).^2)+sum((v2-u2).^2))^0.5)/(k1+k2)

            u1,u2 = new_estimate(data1, data2, k1, k2, 1, 4)
            s1_new[i] += ((sum((v1-u1).^2)+sum((v2-u2).^2))^0.5)/(k1+k2)
            u1,u2 = new_estimate(data1, data2, k1, k2, 0.5, 4)
            s2_new[i] += ((sum((v1-u1).^2)+sum((v2-u2).^2))^0.5)/(k1+k2)
            u1,u2 = new_estimate(data1, data2, k1, k2, 0.25, 4)
            s3_new[i] += ((sum((v1-u1).^2)+sum((v2-u2).^2))^0.5)/(k1+k2)
            u1,u2 = new_estimate(data1, data2, k1, k2, 0.1, 4)
            s4_new[i] += ((sum((v1-u1).^2)+sum((v2-u2).^2))^0.5)/(k1+k2)
        end
    end
    plot(1:10, s1_ditto/100, label="Ditto, \$λ=2\$", lw=2, xlabel="\$K_2/K_1\$", ylabel="test error")
    plot!(1:10, s2_ditto/100, label="Ditto, \$λ=1\$", lw=2)
    plot!(1:10, s3_ditto/100, label="Ditto, \$λ=0.5\$", lw=2)
    plot!(1:10, s4_ditto/100, label="Ditto, \$λ=0.25\$", lw=2)
    plot!(1:10, s1_new/100, label="MLFL, \$λ=(1,4)\$", linestyle=:dashdot)
    plot!(1:10, s2_new/100, label="MLFL, \$λ=(0.5,4)\$", linestyle=:dashdot)
    plot!(1:10, s3_new/100, label="MLFL, \$λ=(0.25,4)\$", linestyle=:dashdot)
    plot!(1:10, s4_new/100, label="MLFL, \$λ=(0.1,4)\$", linestyle=:dashdot)
    savefig("unbalance data.png")
end

#distance between distributions
n=10
k1=10
k2=10
let
    ditto_re1=zeros(21)
    ditto_re2=zeros(21)
    ditto_re3=zeros(21)
    ditto_re4=zeros(21)
    new_re1=zeros(21)
    new_re2=zeros(21)
    new_re3=zeros(21)
    new_re4=zeros(21)
    Random.seed!(1234)
    for i=1:21
        for t=1:100
            data1,v1 = generate_data(0, k1, n, 1);
            data2,v2 = generate_data((i-1)*0.5, k2, n, 1);

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 2)
            ditto_re1[i] += (sum((v1-u1).^2)+sum((v2-u2).^2))^0.5/20
            u1,u2 = ditto_estimate(data1, data2, k1, k2, 1)
            ditto_re2[i] += (sum((v1-u1).^2)+sum((v2-u2).^2))^0.5/20
            u1,u2 = ditto_estimate(data1, data2, k1, k2, 0.5)
            ditto_re3[i] += (sum((v1-u1).^2)+sum((v2-u2).^2))^0.5/20
            u1,u2 = ditto_estimate(data1, data2, k1, k2, 0.25)
            ditto_re4[i] += (sum((v1-u1).^2)+sum((v2-u2).^2))^0.5/20

            #new method with different lambda
            u1,u2 = new_estimate(data1, data2, k1, k2, 1, 4)
            new_re1[i] += (sum((v1-u1).^2)+sum((v2-u2).^2))^0.5/20
            u1,u2 = new_estimate(data1, data2, k1, k2, 0.5, 4)
            new_re2[i] += (sum((v1-u1).^2)+sum((v2-u2).^2))^0.5/20
            u1,u2 = new_estimate(data1, data2, k1, k2, 0.25, 4)
            new_re3[i] += (sum((v1-u1).^2)+sum((v2-u2).^2))^0.5/20
            u1,u2 = new_estimate(data1, data2, k1, k2, 0.1, 4)
            new_re4[i] += (sum((v1-u1).^2)+sum((v2-u2).^2))^0.5/20
        end
    end
    plot(((1:21).-1)*0.5, ditto_re1/100, label="Ditto, \$λ=2\$", lw=2, xlabel="\$θ\$", ylabel="test error", legend=:topleft)
    plot!(((1:21).-1)*0.5, ditto_re2/100, label="Ditto, \$λ=1\$", lw=2)
    plot!(((1:21).-1)*0.5, ditto_re3/100, label="Ditto, \$λ=0.5\$", lw=2)
    plot!(((1:21).-1)*0.5, ditto_re4/100, label="Ditto, \$λ=0.25\$", lw=2)
    plot!(((1:21).-1)*0.5, new_re1/100, label="MLFL, \$λ=(1,4)\$", linestyle=:dashdot)
    plot!(((1:21).-1)*0.5, new_re2/100, label="MLFL, \$λ=(0.5,4)\$", linestyle=:dashdot)
    plot!(((1:21).-1)*0.5, new_re3/100, label="MLFL, \$λ=(0.25,4)\$", linestyle=:dashdot)
    plot!(((1:21).-1)*0.5, new_re4/100, label="MLFL, \$λ=(0.1,4)\$", linestyle=:dashdot)
    savefig("distance checking.png")
end

#robustness checking
n=10
k1=10
k2=10
let
    re_ditto1 = zeros(21)
    re_ditto2 = zeros(21)
    re_new1 = zeros(21)
    re_new2 = zeros(21)
    re_ditto3 = zeros(21)
    re_ditto4 = zeros(21)
    re_new3 = zeros(21)
    re_new4 = zeros(21)
    Random.seed!(1234)
    for i=1:21
        tau = 1+0.2*(i-1)
        s1_ditto = 0
        s2_ditto = 0
        s1_new = 0
        s2_new = 0
        s3_ditto = 0
        s4_ditto = 0
        s3_new = 0
        s4_new = 0
        for t=1:100
            data1,v1 = generate_data(0, k1, n, tau);
            data2,v2 = generate_data(1, k2, n, tau);

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 1)
            u = append!!(v1-u1, v2-u2)
            s1_ditto += (sum(u.^2)).^0.5/(k1+k2)

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 2)
            u = append!!(v1-u1, v2-u2)
            s2_ditto += (sum(u.^2)).^0.5/(k1+k2)

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 0.5)
            u = append!!(v1-u1, v2-u2)
            s3_ditto += (sum(u.^2)).^0.5/(k1+k2)

            u1,u2 = ditto_estimate(data1, data2, k1, k2, 0.25)
            u = append!!(v1-u1, v2-u2)
            s4_ditto += (sum(u.^2)).^0.5/(k1+k2)

            u1,u2 = new_estimate(data1, data2, k1, k2, 0.1, 4)
            u = append!!(v1-u1, v2-u2)
            s1_new += (sum(u.^2)).^0.5/(k1+k2)

            u1,u2 = new_estimate(data1, data2, k1, k2, 0.25, 4)
            u = append!!(v1-u1, v2-u2)
            s2_new += (sum(u.^2)).^0.5/(k1+k2)

            u1,u2 = new_estimate(data1, data2, k1, k2, 0.5, 4)
            u = append!!(v1-u1, v2-u2)
            s3_new += (sum(u.^2)).^0.5/(k1+k2)

            u1,u2 = new_estimate(data1, data2, k1, k2, 1, 4)
            u = append!!(v1-u1, v2-u2)
            s4_new += (sum(u.^2)).^0.5/(k1+k2)
        end
        re_ditto1[i] = s1_ditto/100
        re_new1[i] = s1_new/100
        re_ditto2[i] = s2_ditto/100
        re_new2[i] = s2_new/100
        re_ditto3[i] = s3_ditto/100
        re_new3[i] = s3_new/100
        re_ditto4[i] = s4_ditto/100
        re_new4[i] = s4_new/100
    end
    plot(0.2*((1:21).-1).+1, re_ditto2, label="Ditto, \$λ=2\$", lw=2, legend=:topleft, xlabel="\$τ\$", ylabel="test error")
    plot!(0.2*((1:21).-1).+1, re_ditto1, label="Ditto, \$λ=1\$", lw=2)
    plot!(0.2*((1:21).-1).+1, re_ditto3, label="Ditto, \$λ=0.5\$", lw=2)
    plot!(0.2*((1:21).-1).+1, re_ditto4, label="Ditto, \$λ=0.25\$", lw=2)
    plot!(0.2*((1:21).-1).+1, re_new4, label="MLFL, \$λ=(1,4)\$", linestyle=:dashdot)
    plot!(0.2*((1:21).-1).+1, re_new3, label="MLFL, \$λ=(0.5,4)\$", linestyle=:dashdot)
    plot!(0.2*((1:21).-1).+1, re_new2, label="MLFL, \$λ=(0.25,4)\$", linestyle=:dashdot)
    plot!(0.2*((1:21).-1).+1, re_new1, label="MLFL, \$λ=(0.1,4)\$", linestyle=:dashdot)
    savefig("robustchecking.png")
end
=#
