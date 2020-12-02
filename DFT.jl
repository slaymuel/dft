include("Constants.jl")

module DFT
import ..D
import ..Constants

#Simulation Constants
const T = 298.0
const kT = T * Constants.k_B
const beta = 1.0 / (kT)
const ε_r = 78.3
const l_B = Constants.e^2 / (4.0 * pi * Constants.ε_0 * ε_r) * beta * 10.0^(10.0)
h = 40.0
const sigma = -1.0 / 200.0
const mix = 0.004
const dz = 0.05
const lambda = (8.0*pi*0.01)^(1.0/3.0)

println("l_B is: $l_B")
println("surfdens is: $sigma")

function set_h(hz)
    global h = hz
end

struct Energy{T<:Function}
    phi::T
    phiV::Vector{Float64}
    ext::T
    extV::Vector{Float64}
    hole::T
    holeVP::Vector{Float64}
    holeVN::Vector{Float64}
end


function find_root3p_v(a, b, c, guess)
    p0 = -20.0
    p1 = 20.0

    int = [p0, p1]
    mp = 0.0
    f(x) = (b*exp(4.0 * x) + a*exp(x) + c)

    if(f(p0) >= 0.0 || f(p1) <= 0.0)
        println()
        println("FAILED TO FIND ROOT!!!!!!!")
        println()
        return nothing
    end

    while true
        mp = (int[1] + int[2]) / 2.0
        if(f(mp) < 0.0)
            int[1] = mp
        else
            int[2] = mp
        end
        if(abs(int[1] - int[2]) < 1e-15)
            break
        end
        #if(abs(f(mp)) < 1e-10)
        #    break
        #end
    end

    return mp
    #x0 = guess
    #x1 = 0.0
    #for i in 1:1000000
    #    x1 = x0 - (b*exp(4.0 * x0) + a*exp(x0) + c) / (4.0*b*exp(4.0 * x0) + a*exp(x0))
#
    #    if((b*exp(4.0 * x1) + a*exp(x1) + c) < 1e-12)
    #        break
    #    end
    #    #if(abs(x0 - x1) < 1e-10)
    #    #    break
    #    #end
    #    x0 = x1
    #end
    #println(x1)
    #return x1
end


function find_root3p(a, b, c, guess)
    x0 = guess
    x1 = 0.0
    for i in 1:1000000
        x1 = x0 - (b*x0^(4.0) + a*x0 + c) / (4.0*b*x0^(3.0) + a)
        if(x1 < 0.0)
            guess = guess + 1.0
            x0 = guess
            continue

        elseif(abs(b*x1^(4.0) + a*x1 + c) < 1e-12)
            break
        end
        #if(abs(x0 - x1) < 1e-10)
        #    break
        #end
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

function hole(r::Float64, l::Float64)
    return (1.0 - exp(-l * r))
end


function hole_coloumb(diffz, lambda, q1 = 1.0)
    return (-2.0*pi*l_B*q1^2.0 * (diffz +  exp(-lambda*diffz) / lambda))
    #return energy
end

function hole_coloumb2(z, lambda, q1 = 1.0)
    energy = -2.0 * pi * l_B * q1^2.0 * z * (1.0 - exp(-lambda * z))
    return energy
end

function pre_cal(f, v::Vector{Float64})
    for i in 0:dz:h
        v[round(Int, i / dz) + 1] = f(i)
    end
end

function pre_cal_hole(f, v::Vector{Float64}, lambda)
    for i in 0:dz:h
        v[round(Int, i / dz) + 1] = f(i, lambda)
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

function energy_wh(s::D, z, en::Energy)
    energy = 0.0
    
    #lambda = sqrt(2.0*pi*abs(sigma) / abs(q1))
    #Integrate
    @simd for i in eachindex(s.n)
        diffz = abs(z - (i*dz - dz))
        #@inbounds energy = energy + n[i] * phi(diffz, q1, q2) * (1.0 - exp(-lambda * diffz)) * dz 
        @inbounds energy = energy + s.n[i] * hole_coloumb(diffz, s.lambda, s.q) * dz 
        #@inbounds energy = energy + n[i] * phi(abs(z - i*dz), q1, q2) * hole(diffz, lambda) * dz 
        #-2.d0*pi*rrT*(diffz+dexp(-tlambda*diffz)/tlambda)*dz
        #energy = energy - n[i]*2.0*pi*l_B*q1*q2 * (diffz -  exp(-lambda*diffz) / lambda) * dz 
        #en.holeV[round(Int, abs(z - i*dz) / dz) + 1]
        #energy = energy + n[i] * q1 * q2 * en.phiV[round(Int, abs(z - i*dz) / dz) + 1] * dz
    end
    
    return energy
end

function energy_whp(s::D, z, en::Energy)
    energy = 0.0
    @simd for i in eachindex(s.n)
        diffz = abs(z - (i*dz - dz))
        @inbounds energy = energy + s.n[i] * s.q^2.0* en.holeVP[round(Int, diffz / dz) + 1] * dz
    end
    return energy
end

function energy_whn(s::D, z, en::Energy)
    energy = 0.0
    @simd for i in eachindex(s.n)
        diffz = abs(z - (i*dz - dz))
        @inbounds energy = energy + s.n[i] * s.q^2.0 * en.holeVN[round(Int, diffz / dz) + 1] * dz
    end 
    return energy
end

function update_densities!(p::D, n::D, donn::Float64, en::Energy)
    global mix
    pn_temp = fill(0.0, length(p.n))
    nn_temp = fill(0.0, length(n.n))

    Threads.@threads for i in eachindex(p.n)
        z = i*dz - dz
        pe = energy(n.n, p.q, n.q, z, en) + energy_whp(p, z, en)
        ne = energy(p.n, n.q, p.q, z, en) + energy_whn(n, z, en)
        @inbounds pn_temp[i] = exp(p.cp + p.q * donn - external(z, p.q) - pe)
        @inbounds nn_temp[i] = exp(n.cp + n.q * donn - external(z, n.q) - ne)
    end

    charge = get_charge(pn_temp, nn_temp, p.q, n.q)
    #println(charge)
    #Neutralize
    if(abs(charge) > 1e-10)
        dDonn = log(find_root3p(2.0*sigma, p.q*sum(pn_temp)*dz, n.q*sum(nn_temp)*dz, exp(donn)))
        #dDonn = find_root3p_v(2.0*sigma, p.q*sum(pn_temp)*dz, n.q*sum(nn_temp)*dz, donn)

        pn_temp = pn_temp .* exp(p.q * dDonn)
        nn_temp = nn_temp .* exp(n.q * dDonn)
        #donn = donn + mix * dDonn
        donn = donn + dDonn
    end

    charge = get_charge(pn_temp, nn_temp, p.q, n.q)
    if(abs(charge) > 1e-10)
        println()
        println("Suggestion is not electroneutral!")
        println()
    end

    #for i in eachindex(p.n)
    #    trat = pn_temp[i]/p.n[i]
    #    if (trat > 1.0)
    #        ttfdm = 2.0*p.n[i]-(p.n[i]^(2.0)) / pn_temp[i]
    #        if (ttfdm < pn_temp[i])
    #            pn_temp[i] = ttfdm
    #        end
    #    #else
    #    #    ttfdm = p.n[i]-pn_temp[i]+(pn_temp[i]^(2.0))/p.n[i]
    #    #    if (ttfdm > pn_temp[i])
    #    #        pn_temp[i] = ttfdm
    #    #    end
    #    end
    #end
    #for i in eachindex(n.n)
    #    trat = nn_temp[i]/n.n[i]
    #    if (trat > 1.0)
    #        ttfdm = 2.0*n.n[i]-(n.n[i]^(2.0)) / nn_temp[i]
    #        if (ttfdm < nn_temp[i])
    #            nn_temp[i] = ttfdm
    #        end
    #    #else
    #    #    ttfdm = n.n[i]-nn_temp[i]+(nn_temp[i]^(2.0))/n.n[i]
    #    #    if (ttfdm > nn_temp[i])
    #    #        nn_temp[i] = ttfdm
    #    #    end
    #    end
    #end
    p.n = (mix .* pn_temp .+ (1.0 - mix) .* p.n::Vector{Float64})
    n.n = (mix .* nn_temp .+ (1.0 - mix) .* n.n::Vector{Float64})  

    #tfdm = new density
    #fdm = old

    #trat = tfdm/fdm
    #if (trat.gt.1.d0) then
    #   ttfdm = 2.d0*fdm-fdm**2/tfdm
    #   if (ttfdm.lt.tfdm) then
    #       rat = ttfdm/tfdm
    #       tfdm = ttfdm
    #       tfem = tfem*rat
    #   endif
    #endif
    return donn
end

function get_f(p::D, n::D, en::Energy)
    f = 0.0
    
    #for i in eachindex(p.n)
    for i in 1:1:round(Int, length(p.n) / 2.0)
        z  = i*dz - dz

        f_temp = energy(n.n, p.q, n.q, z, en)          #+-
        f_temp = f_temp + 0.5 * energy_wh(p, z, en) #++


        f = f + p.n[i]*(f_temp + external(z, p.q))
        f = f + p.n[i]*(log(p.n[i]) - 1.0 - p.cp)  #Ideal + legendre GC
    end

    #for i in eachindex(n.n)
    for i in 1:1:round(Int, length(n.n) / 2.0)
        z  = i*dz - dz

        f_temp = 0.5 * energy_wh(n, z, en) #--


        f = f + n.n[i]*(f_temp + external(z, n.q))
        f = f + n.n[i]*(log(n.n[i]) - 1.0 - n.cp)  #Ideal + legendre GC

        #fel = fds*(srnsval*(VL*z+VR*hz)+0.5d0*uss+Vexs(iz))
        #felns = fds*(0.5d0*uss+Vexs(iz))
        #aeta = pis*cdmonm(iz)
        #rxsib = 1.d0/(1.d0-aeta)
        #ccc = -fds+ccc
        #sumFsns = fds*(dlog(fds)-1.d0-chemps)+felns+sumFsns
        #sumFs = fds*(dlog(fds)-1.d0-chemps)+fel+sumFs
    end

    f = 2.0 * f * dz -2.0*pi*h * l_B * sigma^2.0
    return f
end

function get_pressure(p::D, n::D, p_B)

    fp1S = p.n[1]
    fn1S = p.n[1+2]
    c0Skv = p.n[1+1]
    c1Skv = (fp1S-fn1S)*0.5
    c2Skv = (fp1S+fn1S-2.0*c0Skv)*0.5
#      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
    fwcmq = (c0Skv+c1Skv*1.5+c2Skv*2.25)
#      ** EXTRAPOLATION BY USE OF LINEAR EXPRESSION **
    fwcml = fp1S+0.5*(fp1S-c0Skv)
    if (fwcml < 1e-11) fwcmq = 0.0 end
    #println("monomer contact density - quad. extr.: ", fwcmq)
    
    fp1S = n.n[1]
    fn1S = n.n[1 + 2]
    c0Skv = n.n[1 + 1]
    c1Skv = (fp1S - fn1S)*0.5
    c2Skv = (fp1S + fn1S - 2.0*c0Skv)*0.5
#      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
    fwcsq = (c0Skv + c1Skv*1.5 + c2Skv*2.25)
    fwcsl = fp1S+0.5*(fp1S-c0Skv)
    if (fwcsl < 1e-11) fwcsq = 0.0 end
    #println("solvent contact density - quad. extr.: ", fwcsq)
    VL = -2.0*pi*l_B*sigma

    Psubnet = fwcmq + fwcsq - p_B
    println("net solvation force (excl. csurf and surf.-surf.): ", Psubnet)

    Pnet = Psubnet - (p.q*0.5*sum(p.n) + n.q*0.5*sum(n.n))*dz * VL

    return Pnet
end

function run!(p::D, n::D, donn, en::Energy, iters)
    es22mm = -4.0*pi*p.q*p.q*l_B/(p.lambda^(2.0))
    es22ss = -4.0*pi*n.q*n.q*l_B/(n.lambda^(2.0))
    p_B = p.dens+n.dens+0.5*es22mm*p.dens^(2.0)+0.5*es22ss*n.dens^(2.0)

    println("+lambda: ", p.lambda)
    println("-lambda: ", n.lambda)
    println("+cp: ", p.cp)
    println("-cp", n.cp)
    println("Bulk pressure: ", p_B)
    println("sigma (wall charge): ", sigma)

    charge = get_charge(p.n, n.n, p.q, n.q)
    if(abs(charge) > 1e-10)
        dDonn = log(find_root3p(2.0*sigma, p.q * sum(p.n)*dz, n.q * sum(n.n)*dz, exp(donn)))
        #dDonn = find_root3p_v(2.0*sigma, p.q*sum(p.n)*dz, n.q*sum(n.n)*dz, donn)
        #donn = donn + dDonn
        donn = donn + dDonn
        #p.n = p.n .* exp(p.q * dDonn)
        #n.n = n.n .* exp(n.q * dDonn)
    end
    println("Donnan at start: ", donn)
    charge = get_charge(p.n .* exp(p.q*donn), n.n .* exp(n.q*donn), p.q, n.q)
    println("Initial charge: ", charge)

    for i in 0:iters
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
            #pressure = p.n[round(Int, 2.5 / dz)] + n.n[round(Int, 2.5 / dz)]
            #println("Pressure: ", pressure)

        end

        donn = update_densities!(p, n, donn, en)


        #Symmetrize
        for i in 1:1:floor(Int, length(p.n) / 2.0)
            p.n[i] = (p.n[i] + p.n[end-i+1]) / 2.0
            p.n[end-i+1] = p.n[i]
            n.n[i] = (n.n[i] + n.n[end-i+1]) / 2.0
            n.n[end-i+1] = n.n[i]
        end
        if(abs(charge) > 1e-10)
            #dDonn = log(find_root3p(2.0*sigma, p.q * sum(p.n)*dz, n.q * sum(n.n)*dz, exp(donn)))
            #dDonn = find_root3p_v(2.0*sigma, p.q*sum(p.n)*dz, n.q*sum(n.n)*dz, donn)
            #donn = donn + dDonn
            #donn = dDonn
            #p.n = p.n .* exp(p.q * dDonn)
            #n.n = n.n .* exp(n.q * dDonn)
        end
    end
    dp = get_pressure(p, n, p_B)
    gp = get_f(p, n, en)
    println()
    println("Grand Potential: ", gp)
    println("Disjoining pressure: ", dp)
    println("Final charge sum: ", p.q*sum(p.n)*DFT.dz + n.q*sum(n.n)*DFT.dz )
    return donn, dp, gp
end
end