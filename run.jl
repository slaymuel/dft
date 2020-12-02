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



function read_dens(h, pdens, ndens)
    h_vals = readdlm("../../Documents/dft_prog/pDens_150.0_-0.0025_0.1.txt", '\t', Float64, '\n')[:,1]
    p_dens = readdlm("../../Documents/dft_prog/pDens_150.0_-0.0025_0.1.txt", '\t', Float64, '\n')[:,2]
    n_dens = readdlm("../../Documents/dft_prog/nDens_150.0_-0.0025_0.1.txt", '\t', Float64, '\n')[:,2]

    vals = round(Int, h / DFT.dz) + 1

    if(vals > length(h_vals))
        println("Increasing size of density vector to: $vals")
        for i in 1:1:(vals - length(h_vals))
            insert!(p_dens, round(Int, length(h_vals) / 2.0), pdens)
            insert!(n_dens, round(Int, length(h_vals) / 2.0), ndens)
        end
    else
        for i in 1:1:(length(h_vals) - vals)
            deleteat!(p_dens, round(Int, length(p_dens) / 2.0))
            deleteat!(n_dens, round(Int, length(p_dens) / 2.0))
        end
    end

    @assert length(p_dens) == vals
    @assert length(n_dens) == vals
    return p_dens, n_dens
end

function main(sv)
    ps = zeros(0)
    gps = zeros(0)
    for i in 1:length(sim_vars["conc"])
        DFT.set_h(sv["h"][i])
        println("h: ", sv["h"][i])
        pq = 3.0
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

        #p_lambda = (8.0*Constants.pi*dens)^(1.0/3.0)
        #p_lambda = (2.0*Constants.pi*abs(DFT.sigma)/pq)^0.5
        #n_lambda = (2.0*Constants.pi*abs(DFT.sigma))^0.5 

        #es22mm = -4.0*pi*p.q*p.q*l_B/(p.lambda^(2.0))
        #es22ss = -4.0*pi*n.q*n.q*l_B/(n.lambda^(2.0))
        #chempp = dlog(bdm)+bdm*es22mm
        #chemps = dlog(bds)+es22ss*bds

        pcp = pcp - 4.0*DFT.pi*pq*pq*DFT.l_B/(p_lambda^2.0) * dens
        println("pcp: ", pcp)
        ncp = ncp - 4.0*DFT.pi    *  DFT.l_B/(n_lambda^2.0) * pq * dens

        scd_corr = -2.0*DFT.sigma / (sv["h"][i] * pq)
        println("cbdm: ", scd_corr)
        println("bdm: ", dens)
        println("bds: ", pq * dens)
        pdens = dens + scd_corr# * DFT.h
        ndens = pq * dens# * DFT.h
        
                                        #Density                  q    cp        l         
        p = D(fill(pdens, round(Int, DFT.h / DFT.dz) + 1), 
            pq,   pcp, dens, p_lambda)
        n = D(fill(ndens, round(Int, DFT.h / DFT.dz) + 1), -1.0, 
            ncp, pq * dens, n_lambda)
        println("Reading previous density")
        p_dens, n_dens = read_dens(sv["h"][i], pdens, ndens)
        p.n = p_dens
        n.n = n_dens
        println("Density loaded")
        open("../../Documents/dft_prog/pDens_initial", "w") do io
            writedlm(io, [0:DFT.dz:sv["h"][i] p.n])
        end
        open("../../Documents/dft_prog/nDens_initial", "w") do io
            writedlm(io, [0:DFT.dz:sv["h"][i] n.n])
        end
        en = DFT.Energy{Function}(DFT.phi, fill(0.0, length(p.n)), DFT.external, fill(0.0, length(p.n)),
            DFT.hole_coloumb, fill(0.0, length(p.n)), fill(0.0, length(n.n)))
        DFT.pre_cal_hole(en.hole, en.holeVP, p.lambda)
        DFT.pre_cal_hole(en.hole, en.holeVN, n.lambda)

        DFT.pre_cal(en.phi, en.phiV)
        DFT.pre_cal(en.ext, en.extV)

        #println(en.phiV)
        #println(isabstracttype(en))
        @time donnan, dp, gp = DFT.run!(p, n, donnan, en, 5000)
        append!(ps, dp)
        append!(gps, gp)
 
        open("../../Documents/dft_prog/pDens_" * string(sv["h"][i]) * "_" * string(DFT.sigma) * "_" * 
                                                        string(sim_vars["conc"][i]) * ".txt", "w") do io
            writedlm(io, [0:DFT.dz:sv["h"][i] p.n])
        end
        open("../../Documents/dft_prog/nDens_" * string(sv["h"][i]) * "_" * string(DFT.sigma) * "_" * 
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


concs = fill(0.1, 1)
hs = 200.0:10.0:200.0
#donnans = fill(0.7259712622913146, 1)
donnans = fill(0.5664351792981909, 1)
               
               
println(length(concs))
println(length(hs))
println(length(donnans))
sim_vars = Dict("conc"=>concs, 
                "h"=>hs,
                "donnan" => donnans)


ps, gps = main(sim_vars)

open("../../Documents/dft_prog/grand_potential.txt", "w") do io
    writedlm(io, [hs gps])
end
open("../../Documents/dft_prog/Pnet.txt", "w") do io
    writedlm(io, [hs ps])
end

using DelimitedFiles
using Interpolations

function pl()
    #hs = readdlm("../../Documents/dft_prog/0.005_0.0018/Pnet.txt", '\t', Float64, '\n')[:,1]
    #ps = readdlm("../../Documents/dft_prog/0.005_0.0018/Pnet.txt", '\t', Float64, '\n')[:,2]
    #println(hs)
    #gps = readdlm("../../Documents/dft_prog/0.005_0.0018/grand_potential.txt", '\t', Float64, '\n')[:,2]

    hs = readdlm("../../Documents/dft_prog/Pnet.txt", '\t', Float64, '\n')[:,1]
    ps = readdlm("../../Documents/dft_prog/Pnet.txt", '\t', Float64, '\n')[:,2]
    println(hs)
    gps = readdlm("../../Documents/dft_prog/grand_potential.txt", '\t', Float64, '\n')[:,2]
    der = zeros(0)
    xs = zeros(0)
    for i in 1:1:length(hs) - 1
        dyf = (gps[i + 1] + gps[i]) / 2.0 - gps[i + 1]
        dyb = (gps[i + 1] + gps[i]) / 2.0 - gps[i]
        #x = (hs[i + 1] + hs[i]) / 2.0
        x = (hs[i + 1] + hs[i]) / 2.0
        dxf = x - hs[i + 1]
        dxb = x - hs[i]
        dxyf = -dyf / dxf
        dxyb = -dyb / dxb
        append!(der, (dxyf + dxyb) / 2.0)
        append!(xs, x)
    end

    #A = gps
    #A_x = 10.0:2.0:150.0
    #nodes = (A_x,)
    #itp = interpolate(nodes, A, Gridded(Linear()))

    dgps = -diff(gps) ./ diff(hs)

    PyPlot.clf()
    plot(hs, ps, linewidth=2.0, color="red")

    #plot(hs[1:end-1], dgps)
    plot(xs, der, linewidth=2.0)
    #plot(xs[2:end], dgps2, color="green")
    #plot(hs[1:end], gps, color="blue")
    #plot(xs, itp(xs), color="red", linewidth=2.0, linestyle="--")
    xlim(20.0, 150.0)
    display(gcf())

    savefig("../../Documents/dft_prog/dis_press.png")
end


function gs()
    #hs = readdlm("../../Documents/dft_prog/0.005_0.0018/Pnet.txt", '\t', Float64, '\n')[:,1]
    #ps = readdlm("../../Documents/dft_prog/0.005_0.0018/Pnet.txt", '\t', Float64, '\n')[:,2]
    #println(hs)
    #gps = readdlm("../../Documents/dft_prog/0.005_0.0018/grand_potential.txt", '\t', Float64, '\n')[:,2]

    hs = readdlm("../../Documents/dft_prog/aurora/0.0018M_scd-0.005/grand_potential.txt", '\t', Float64, '\n')[:,1]
    ps = readdlm("../../Documents/dft_prog/aurora/0.0018M_scd-0.005/grand_potential.txt", '\t', Float64, '\n')[:,2]
    gps = readdlm("../../Documents/dft_prog/aurora/0.0018M_scd-0.005/grand_potential.txt", '\t', Float64, '\n')[:,2]
    println(hs)

    PyPlot.clf()
    plot(hs, (gps  .- gps[end]) ./ (hs), linewidth=2.0, color="red")

    xlim(20.0, 400.0)
    ylim(-0.00001, 0.000035)
    display(gcf())

    savefig("../../Documents/dft_prog/grand_0.0018.png")
    println((gps  .- gps[end]) ./ hs)
end
#gs()

function appsd()
    scd = -0.005
    hs = readdlm("../../Documents/dft_prog/aurora/0.0006M_scd-0.005/pDens_400.0_-0.005_0.0006.txt", 
                                                                                '\t', Float64, '\n')[:,1]
    pn = readdlm("../../Documents/dft_prog/aurora/0.0006M_scd-0.005/pDens_400.0_-0.005_0.0006.txt", 
                                                                                '\t', Float64, '\n')[:,2]
    nn = readdlm("../../Documents/dft_prog/aurora/0.0006M_scd-0.005/nDens_400.0_-0.005_0.0006.txt", 
                                                                                '\t', Float64, '\n')[:,2]

    p_sum = scd
    n_sum = 0.0
    app_scd = zeros(0)
    for i in eachindex(pn)
        p_sum = p_sum + 3.0*pn[i]*DFT.dz
        n_sum = n_sum - 1.0*nn[i]*DFT.dz
        append!(app_scd, p_sum + n_sum)
    end
    PyPlot.clf()
    plot(hs, app_scd, linewidth=2.0, color="red")

    ylim(-0.00008, 0.00008)
    #ylim(-0.00001, 0.000035)
    display(gcf())
    savefig("../../Documents/dft_prog/app_scd_0.0006.png")
end
#appsd()