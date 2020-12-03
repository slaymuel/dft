cd("/Users/samuel/Documents/dft_prog")

using PyPlot
using DelimitedFiles
include("Constants.jl")
using .Constants

mutable struct D
    n::Vector{Float64}
    q::Float64
    cp::Float64
    dens::Float64
    lambda::Float64
end

include("DFT.jl")
include("io.jl")

function main(sv)
    ps = zeros(0)
    gps = zeros(0)
    for i in 1:length(sim_vars["conc"])
        DFT.set_h(sv["h"][i])
        println("h: ", sv["h"][i])
        pq = 2.0
        conc = sv["conc"][i] #Molars
        donnan = sv["donnan"][i];
        dens = conc * Constants.N_A * 10^(-27.0)
        pcp = log(dens)
        ncp = log(pq * dens)

        n_lambda = 0.0
        p_lambda = 0.0
        pl1 = (8.0*Constants.pi*dens)^(1.0/3.0)
        pl2 = (2.0*Constants.pi*abs(DFT.sigma)/pq)^0.5
        if (pl1 > pl2)
            p_lambda = pl1
        else
            p_lambda = pl2
        end
        nl1 = (8.0*Constants.pi*pq*dens)^(1.0/3.0)
        nl2 = (2.0*Constants.pi*abs(DFT.sigma))^0.5 
        if (nl1 > nl2)
            n_lambda = nl1
        else
            n_lambda = nl2
        end

        pcp = pcp - 4.0*DFT.pi*pq*pq*DFT.l_B/(p_lambda^2.0) * dens
        println("pcp: ", pcp)
        ncp = ncp - 4.0*DFT.pi    *  DFT.l_B/(n_lambda^2.0) * pq * dens

        scd_corr = -2.0*DFT.sigma / (sv["h"][i] * pq)
        println("cbdm: ", scd_corr)
        println("bdm: ", dens)
        println("bds: ", pq * dens)
        pdens = dens + scd_corr
        ndens = pq * dens
        
                                        #Density                  q    cp        l         
        p = D(fill(pdens, round(Int, DFT.h / DFT.dz) + 1), 
            pq,   pcp, dens, p_lambda)
        n = D(fill(ndens, round(Int, DFT.h / DFT.dz) + 1), -1.0, 
            ncp, pq * dens, n_lambda)

        if(sv["p_file"][i] != nothing)
            println("Reading previous density")
            p_dens, n_dens = read_dens(sv["h"][i], sv["p_file"][i], sv["n_file"][i], pdens, ndens)
            p.n = p_dens
            n.n = n_dens
            println("Density loaded")
        end

        open(pwd() * "/pDens_initial", "w") do io
            writedlm(io, [0:DFT.dz:sv["h"][i] p.n])
        end
        open(pwd() * "/nDens_initial", "w") do io
            writedlm(io, [0:DFT.dz:sv["h"][i] n.n])
        end
        en = DFT.Energy{Function}(DFT.phi, fill(0.0, length(p.n)), DFT.external, fill(0.0, length(p.n)),
            DFT.hole_coloumb, fill(0.0, length(p.n)), fill(0.0, length(n.n)))
        DFT.pre_cal_hole(en.hole, en.holeVP, p.lambda)
        DFT.pre_cal_hole(en.hole, en.holeVN, n.lambda)

        DFT.pre_cal(en.phi, en.phiV)
        DFT.pre_cal(en.ext, en.extV)

        #println(isabstracttype(en))
        @time donnan, dp, gp = DFT.run!(p, n, donnan, en, 1000)
        append!(ps, dp)
        append!(gps, gp)
        if(isnan(gp))
            break
        end
        open(pwd() * "/pDens_" * string(sv["h"][i]) * "_" * string(DFT.sigma) * "_" * 
                                                        string(sim_vars["conc"][i]) * ".txt", "w") do io
            writedlm(io, [0:DFT.dz:sv["h"][i] p.n])
        end
        open(pwd() * "/nDens_" * string(sv["h"][i]) * "_" * string(DFT.sigma) * "_" * 
                                                        string(sim_vars["conc"][i]) * ".txt", "w") do io
            writedlm(io, [0:DFT.dz:sv["h"][i] n.n])
        end
        println()
        #PLOTTING
        zs = [0.0:DFT.dz:DFT.h;];

        data = reshape(readdlm("../../Downloads/jan_DFT/fcdfil"),:,3)
        
        PyPlot.clf()
        
        plot(data[:,1], data[:,2], color="green", linewidth=2.0, linestyle="-")
        #plot(data[:,1], data[:,3], color="blue", linewidth=2.0, linestyle="-")
        
        plot(zs, p.n, color="red", linewidth=1.0, linestyle="-")
        xlim(0, 40.0)
        #plot(zs, n.n, color="black", linewidth=1.0, linestyle="-")
        plot(size=(200,200))
        display(gcf())
    end 
    return ps, gps
end


concs = fill(0.001, 1)
hs = 50.0:10.0:50.0
#donnans = fill(0.7259712622913146, 1)
donnans = fill(2.437938353874042, 1)
               
#0.21473047053488292          
println(length(concs))
println(length(hs))
println(length(donnans))
sim_vars = Dict("conc"=>concs, 
                "h"=>hs,
                "donnan" => donnans,
                "p_file" => ["pDens_50.0_-0.005_0.001.txt"],
                "n_file"=> ["nDens_50.0_-0.005_0.001.txt"])


ps, gps = main(sim_vars)

open(pwd() * "/grand_potential.txt", isfile(pwd() * "/grand_potential.txt") ? "a" : "w") do io
    writedlm(io, [hs gps])
end
open(pwd() * "/Pnet.txt", isfile(pwd() * "/Pnet.txt") ? "a" : "w") do io
    writedlm(io, [hs ps])
end