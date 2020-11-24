using PyPlot
using DelimitedFiles

module DFT
import ..Constants
#Simulation Constants
const T = 298.0
const kT = T * Constants.k_B
const beta = 1.0 / (kT)
const ε_r = 78.3
const l_B = Constants.e^2 / (4.0 * pi * Constants.ε_0 * ε_r) * beta * 10^(10.0)
const h = 20.0
const sigma = -1.0 / 200.0
const mix = 0.01
const dz = 0.05
const lambda = (8.0*pi*0.01)^(1.0/3.0)

println("l_B is: $l_B")
println("surfdens is: $sigma")

mutable struct D
    n::Vector{Float64}
    q::Float64
    cp::Float64
    l::Int32
end


struct Energy{T<:Function}
    phi::T
    phiV::Vector{Float64}
    ext::T
    extV::Vector{Float64}
    hole::T
    holeV::Vector{Float64}
end


function find_root3p(a, b, c, guess)
    x0 = guess
    x1 = 0.0
    for i in 1:1000000
        x1 = x0 - (b*x0^(4.0) + a*x0 + c) / (4.0*b*x0^(3.0) + a)

        if(abs(x0 - x1) < 1e-10)
            break
        end
        x0 = x1
    end
    return x1
end

function find_root1p(a, b, c, guess)
    x0 = guess
    x1 = 0.0
    for i in 1:100
        x1 = x0 - (b*x0^(2.0) + a*x0 + c) / (2.0*b*x0 + a)
        if(abs(x0 - x1) < 1e-10)
            break
        end
        x0 = x1
    end
    return x1
end

function phi(r, q1 = 1.0, q2 = 1.0)
    #smphiz = 2.d0*pi*diffz*rrT*dz
    return -2.0 * pi * l_B * q1 * q2 * r
end


function phi_coulomb(r, q1 = 1.0, q2 = 1.0)
    return l_B * q1 * q2 / r
end

function external_diff(z, q = 1.0)
    return -2.0 * pi * l_B * q * (sigma*z + sigma*(h - z))
end


function external(z, q = 1.0)
    return -2.0 * pi * l_B * q * sigma * h
    #return 0.0
end



function get_charge(p::Vector{Float64}, n::Vector{Float64}, pq, nq)
    charge = (pq * sum(p) + nq * sum(n)) * dz
    return charge + 2.0 * sigma
end

function hole(r::Float64)
    return (1.0 - exp(-lambda * r))
end

function hole(r::Float64, lambda::Float64)
    return (1.0 - exp(-lambda * r))
end


function hole_coloumb(diffz, lambda, q1 = 1.0, q2 = 1.0)
    energy = -2.0*pi*l_B*q1*q2 * (diffz +  exp(-lambda*diffz) / lambda)
    return energy
end

function hole_coloumb2(z, lambda, q1 = 1.0, q2 = 1.0)
    energy = -2.0 * pi * l_B * q1 * q2 * z * (1.0 - exp(-lambda * z))
    return energy
end

function pre_cal(f, v::Vector{Float64})
    for i in 0:dz:h
        v[round(Int, i / dz) + 1] = f(i)
    end
end

function energy(n::Vector{Float64}, q1, q2, z, en::Energy)
    energy = 0.0
    
    #Integrate
    @simd for i in eachindex(n)
        diffz = abs(z - (i*dz - dz))
        @inbounds energy = energy + n[i] * phi(diffz, q1, q2) * dz
        #energy = energy + n[i] * q1 * q2 * en.phiV[round(Int, abs(z - i*dz) / dz) + 1] * dz
    end

    return energy
end

function energy_wh(n::Vector{Float64}, q1, q2, z, en::Energy)
    energy = 0.0
    
    lambda = sqrt(2.0*pi*abs(sigma) / abs(q1))
    #Integrate
    @simd for i in eachindex(n)
        diffz = abs(z - (i*dz - dz))
        #@inbounds energy = energy + n[i] * phi(diffz, q1, q2) * (1.0 - exp(-lambda * diffz)) * dz 
        @inbounds energy = energy + n[i] * hole_coloumb(diffz, lambda, q1, q2) * dz 
        #@inbounds energy = energy + n[i] * phi(abs(z - i*dz), q1, q2) * hole(diffz, lambda) * dz 
        #-2.d0*pi*rrT*(diffz+dexp(-tlambda*diffz)/tlambda)*dz
        #energy = energy - n[i]*2.0*pi*l_B*q1*q2 * (diffz -  exp(-lambda*diffz) / lambda) * dz 
        #en.holeV[round(Int, abs(z - i*dz) / dz) + 1]
        #energy = energy + n[i] * q1 * q2 * en.phiV[round(Int, abs(z - i*dz) / dz) + 1] * dz
    end
    
    return energy
end


function update_densities!(p::D, n::D, donn::Float64, en::Energy)
    global mix
    pn_temp = fill(0.0, p.l)
    nn_temp = fill(0.0, n.l)

    @simd for i in eachindex(p.n)
        z = i*dz - dz
        pe = energy(n.n, p.q, n.q, z, en) + energy_wh(p.n, p.q, p.q, z, en)
        ne = energy(p.n, n.q, p.q, z, en) + energy_wh(n.n, n.q, n.q, z, en)
        @inbounds pn_temp[i] = exp(p.cp + p.q * donn - external(z, p.q) - pe)
        @inbounds nn_temp[i] = exp(n.cp + n.q * donn - external(z, n.q) - ne)
    end

    charge = get_charge(pn_temp, nn_temp, p.q, n.q)

    #Neutralize
    if(abs(charge) > 1e-10)
        dDonn = log(find_root3p(2.0*sigma, p.q*sum(pn_temp)*dz, n.q*sum(nn_temp)*dz, exp(donn)))
        pn_temp = pn_temp .* exp(p.q * dDonn)
        nn_temp = nn_temp .* exp(n.q * dDonn)
        donn = donn + dDonn
    end
    p.n = (mix .* pn_temp .+ (1.0 - mix) .* p.n::Vector{Float64})
    n.n = (mix .* nn_temp .+ (1.0 - mix) .* n.n::Vector{Float64})  
    return donn
end


function run!(p::D, n::D, donn, en::Energy)
    for i in 0:2000
        charge = get_charge(p.n, n.n, p.q, n.q)
        if(i % 100 == 0)
            println()
            println("Iteration: $i")
            println("Charge: $charge")
            println("Donnan: $donn")
            println("sum pn: ", sum(p.n))
            println("sum nn: ", sum(n.n))
            if(isnan(charge))
                println("Charge is NaN exiting...")
                break
            end
        end

        if(abs(charge) > 1e-10)
            dDonn = log(find_root3p(2.0*sigma, p.q * sum(p.n)*dz, n.q * sum(n.n)*dz, exp(donn)))
            donn = donn + dDonn

            #p.n = p.n .* exp(p.q * dDonn)
            #n.n = n.n .* exp(n.q * dDonn)
        end

        donn = update_densities!(p, n, donn, en)

    end
    return donn
end
end

pq = 3.0
conc = 0.001 #Molars
donnan = 5.0;

dens = conc * Constants.N_A * 10^(-27.0)
pcp = log(dens)
ncp = log(pq * dens)

p_lambda = (2.0*Constants.pi*abs(DFT.sigma)/pq)^0.5
n_lambda = (2.0*Constants.pi*abs(DFT.sigma))^0.5 
println("p_lambda: $p_lambda")
println("n_lambda: $n_lambda")
pcp = pcp - 4.0*DFT.pi*pq*pq*DFT.l_B/p_lambda^2.0 * dens
ncp = ncp - 4.0*DFT.pi    *  DFT.l_B/n_lambda^2.0 * pq * dens
println("Cation chemical potential is: ", pcp)
println("Anion chemical potential is: ", ncp)
                                  #Density                  q    cp        l         
p = DFT.D(fill(conc * Constants.N_A * 10^(-27.0)*DFT.h, round(Int, DFT.h / DFT.dz) + 1), pq,   pcp, round(Int,DFT. h / DFT.dz) + 1)
n = DFT.D(fill(pq * conc * Constants.N_A * 10^(-27.0)*DFT.h, round(Int, DFT.h / DFT.dz) + 1), -1.0, ncp, round(Int, DFT.h / DFT.dz) + 1)
en = DFT.Energy{Function}(DFT.phi, fill(0.0, p.l), DFT.external, fill(0.0, p.l), DFT.hole, fill(0.0, p.l))
DFT.pre_cal(en.hole, en.holeV)
DFT.pre_cal(en.phi, en.phiV)
DFT.pre_cal(en.ext, en.extV)

#println(en.phiV)
#println(isabstracttype(en))
@time donnan = DFT.run!(p, n, donnan, en)
println("Final charge sum: ", p.q*sum(p.n)*DFT.dz + n.q*sum(n.n)*DFT.dz )


#PLOTTING
zs = [0.0:DFT.dz:DFT.h;];

data = reshape(readdlm("../../Downloads/jan_DFT/fcdfil"),:,3)

PyPlot.clf()

plot(data[:,1], data2[:,2], color="green", linewidth=2.0, linestyle="-")
plot(data[:,1], data2[:,3], color="blue", linewidth=2.0, linestyle="-")

plot(zs, p.n, color="red", linewidth=1.0, linestyle="-")
plot(zs, n.n, color="black", linewidth=1.0, linestyle="-")
xlim(0, 20)
plot(size=(200,200))
display(gcf())


#plot(zs, nn, color="blue", linewidth=1.0, linestyle="-")

#println(( 3.0*sum(data2[:,2]) - sum(data2[:,3]) ) * DFT.dz)
#println(( 3.0*sum(p.n)        - sum(n.n)        ) * DFT.dz)
#
#sigma = -1.0 / 200.0
#lambda = sqrt(2.0*pi*abs(sigma) / 1.0)
#
#DFT.phi(DFT.dz) * DFT.dz + DFT.hole(DFT.dz, lambda) * DFT.dz
#println("+-: ",-DFT.phi(DFT.dz * DFT.dz, 3.0, -1.0))
#
#println("--: ", DFT.hole_coloumb(DFT.dz, lambda, -1.0, -1.0) * DFT.dz )
#lambda = sqrt(2.0*pi*abs(sigma) / 3.0)
#println("++: ", DFT.hole_coloumb(DFT.dz, lambda, 3.0, 3.0) * DFT.dz)