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
        p_lambda = (2.0*Constants.pi*abs(DFT.sigma)/pq)^0.5
        n_lambda = (2.0*Constants.pi*abs(DFT.sigma))^0.5 
        pcp = pcp - 4.0*DFT.pi*pq*pq*DFT.l_B/p_lambda^2.0 * dens
        ncp = ncp - 4.0*DFT.pi    *  DFT.l_B/n_lambda^2.0 * pq * dens

                                        #Density                  q    cp        l         
        p = D(fill(dens *DFT.h, round(Int, DFT.h / DFT.dz) + 1), 
            pq,   pcp, dens, p_lambda)
        n = D(fill(pq * dens * DFT.h, round(Int, DFT.h / DFT.dz) + 1), -1.0, 
            ncp, pq * dens, n_lambda)

        en = DFT.Energy{Function}(DFT.phi, fill(0.0, length(p.n)), DFT.external, fill(0.0, length(p.n)),
            DFT.hole_coloumb, fill(0.0, length(p.n)), fill(0.0, length(n.n)))
        DFT.pre_cal_hole(en.hole, en.holeVP, p.lambda)
        DFT.pre_cal_hole(en.hole, en.holeVN, n.lambda)

        DFT.pre_cal(en.phi, en.phiV)
        DFT.pre_cal(en.ext, en.extV)
        
        #println(en.phiV)
        #println(isabstracttype(en))
        @time donnan, dp, gp = DFT.run!(p, n, donnan, en, 3000)
        append!(ps, dp)
        append!(gps, gp)
        #PLOTTING
        zs = [0.0:DFT.dz:DFT.h;];

        data = reshape(readdlm("../../Downloads/jan_DFT/fcdfil"),:,3)
        
        PyPlot.clf()
        
        plot(data[:,1], data[:,2], color="green", linewidth=2.0, linestyle="-")
        #plot(data[:,1], data[:,3], color="blue", linewidth=2.0, linestyle="-")
        
        plot(zs, p.n, color="red", linewidth=1.0, linestyle="-")
        xlim(0, 20.0)
        #plot(zs, n.n, color="black", linewidth=1.0, linestyle="-")
        plot(size=(200,200))
        display(gcf())
    end 
    return ps, gps
end

#mutable struct Result
#    h::Float64
#    conc::Float64
#    donnan::Float64
#    dp::Float64
#    gp::Float64
#end
#
#function main2(result::Result)
#    DFT.set_h(result.h)
#    println("h: ", result.h)
#    pq = 3.0
#    conc = result.conc #Molars
#    donnan = result.donnan;
#    dens = conc * Constants.N_A * 10^(-27.0)
#    pcp = log(dens)
#    ncp = log(pq * dens)
#    p_lambda = (2.0*Constants.pi*abs(DFT.sigma)/pq)^0.5
#    n_lambda = (2.0*Constants.pi*abs(DFT.sigma))^0.5 
#    pcp = pcp - 4.0*DFT.pi*pq*pq*DFT.l_B/p_lambda^2.0 * dens
#    ncp = ncp - 4.0*DFT.pi    *  DFT.l_B/n_lambda^2.0 * pq * dens
#
#                                    #Density                  q    cp        l         
#    p = D(fill(conc * dens *DFT.h, round(Int, DFT.h / DFT.dz) + 1), 
#        pq,   pcp, dens, p_lambda)
#    n = D(fill(pq * dens * DFT.h, round(Int, DFT.h / DFT.dz) + 1), -1.0, 
#        ncp, pq * dens, n_lambda)
#
#    en = DFT.Energy{Function}(DFT.phi, fill(0.0, length(p.n)), DFT.external, fill(0.0, length(p.n)),
#        DFT.hole, fill(0.0, length(p.n)))
#    DFT.pre_cal(en.hole, en.holeV)
#    DFT.pre_cal(en.phi, en.phiV)
#    DFT.pre_cal(en.ext, en.extV)
#
#    #println(en.phiV)
#    #println(isabstracttype(en))
#    @time donnan, dp, gp = DFT.run!(p, n, donnan, en, 10000)
#    result.dp = dp
#    result.gp = gp
#end


concs = fill(0.001, 11)
hs = 10.0:2.0:30.0
donnans = 5.0:1.0:15.0
println(length(concs))
println(length(hs))
println(length(donnans))
sim_vars = Dict("conc"=>concs, 
                "h"=>hs,
                "donnan" => donnans)

#resVec = Vector{Result}(undef, length(sim_vars["conc"]))
#
#for i in 1:length(sim_vars["conc"])
#    #res = Dict("conc"=>sim_vars["conc"][i], "h"=>sim_vars["h"][i], "donnan"=>sim_vars["donnan"][i])
#    resVec[i] = Result(sim_vars["h"][i], sim_vars["conc"][i], sim_vars["donnan"][i], 0.0, 0.0)
#    t = @task begin
#        main2(resVec[i])
#    end
#
#    schedule(t)
#
#end
#@sync println("Finished")
#for t in tasks
#    wait(t)
#end

ps, gps = main(sim_vars)

using DelimitedFiles
using Interpolations

function pl(ps, gps)
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

    A = gps
    A_x = 10.0:2.0:30.0
    nodes = (A_x,)
    itp = interpolate(nodes, A, Gridded(Linear()))

    dgps = -diff(gps) ./ diff(hs)

    PyPlot.clf()
    plot(hs, ps, linewidth=2.0, color="red")

    #plot(hs[1:end-1], dgps)
    plot(xs, der, linewidth=2.0)
    #plot(xs[2:end], dgps2, color="green")
    #plot(hs[1:end], gps, color="blue")
    #plot(xs, itp(xs), color="red", linewidth=2.0, linestyle="--")
    xlim(10.0, 30.0)
    display(gcf())

    savefig("../../Documents/dft_prog/dis_press.png")
end

open("../../Documents/dft_prog/grand_potential.txt", "w") do io
    writedlm(io, [hs gps])
end
open("../../Documents/dft_prog/Pnet.txt", "w") do io
    writedlm(io, [hs ps])
end








#plot(zs, nn, color="blue", linewidth=1.0, linestyle="-")

#println(( 3.0*sum(data2[:,2]) - sum(data2[:,3]) ) * DFT.dz)
#println(( 3.0*sum(p.n)        - sum(n.n)        ) * DFT.dz)
#
#sigma = -1.0 / 200.0
#lambda = sqrt(2.0*pi*abs(sigma) / 1.0)
#
#DFT.phi(DFT.dz) * DFT.dz + DFT.hole(DFT.dz, lambda) * DFT.dz
#println("+-: ",-DFT.phi(DFT.dz * DFT.dz, 3.0, -1.0))
##
#println("--: ", DFT.hole_coloumb(DFT.dz, lambda, -1.0, -1.0) * DFT.dz )
#lambda = sqrt(2.0*pi*abs(sigma) / 3.0)
#println("++: ", DFT.hole_coloumb(DFT.dz, lambda, 3.0, 3.0) * DFT.dz)
#
#data = reshape(readdlm("../../Downloads/jan_DFT/fcdfil"),:,3)
#data2 = reshape(readdlm("../../Downloads/jan_DFT/fcdfil_1"),:,3)
#    
#PyPlot.clf()
#
#plot(data[:,1], data2[:,2], color="green", linewidth=2.0, linestyle="-")
#plot(data[:,1], data2[:,3], color="blue", linewidth=2.0, linestyle="-")
#
#plot(data2[:,1], data2[:,2], linewidth=1.0, linestyle="-")
#plot(data2[:,1], data2[:,3], linewidth=1.0, linestyle="-")
#display(gcf())