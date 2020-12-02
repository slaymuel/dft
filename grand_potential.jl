
using PyPlot
include("../../Documents/dft_prog/DFT.jl")

function read_density(p_file, n_file)
    h_vals = readdlm(p_file, '\t', Float64, '\n')[:,1]
    p_dens = readdlm(p_file, '\t', Float64, '\n')[:,2]
    n_dens = readdlm(n_file, '\t', Float64, '\n')[:,2]

    return p_dens, n_dens
end


function get_gp(p, n, en)
    f = 0.0
    
    #for i in eachindex(p.n)
    for i in 1:1:round(Int, length(p.n) / 2.0)
        z  = i*DFT.dz - DFT.dz

        f_temp = DFT.energy(n.n, p.q, n.q, z, en)          #+-
        f_temp = f_temp + 0.5 * DFT.energy_whp(p, z, en) #++

        f = f + p.n[i]*(f_temp + DFT.external(z, p.q))
        f = f + p.n[i]*(log(p.n[i]) - 1.0 - p.cp)  #Ideal + legendre GC
    end

    #for i in eachindex(n.n)
    for i in 1:1:round(Int, length(n.n) / 2.0)
        z  = i*DFT.dz - DFT.dz

        f_temp = 0.5 * DFT.energy_whn(n, z, en) #--

        f = f + n.n[i]*(f_temp + DFT.external(z, n.q))
        f = f + n.n[i]*(log(n.n[i]) - 1.0 - n.cp)  #Ideal + legendre GC
    end

    f = 2.0 * f * DFT.dz -2.0*pi*DFT.h * DFT.l_B * DFT.sigma^2.0
    return f
end

function gp(p_file, n_file, conc)
    h_box = parse(Float64, match(r"pDens\_([0-9]+)", p_file)[1])
    DFT.set_h(h_box)
    println("h: ", h_box)
    pq = 3.0

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
    ncp = ncp - 4.0*DFT.pi    *  DFT.l_B/(n_lambda^2.0) * pq * dens
    scd_corr = -2.0*DFT.sigma / (h_box * pq)
    pdens = dens + scd_corr# * DFT.h
    ndens = pq * dens# * DFT.h
    
                                    #Density                  q    cp        l         
    p = D(fill(pdens, round(Int, DFT.h / DFT.dz) + 1), 
        pq,   pcp, dens, p_lambda)
    n = D(fill(ndens, round(Int, DFT.h / DFT.dz) + 1), -1.0, 
        ncp, pq * dens, n_lambda)
    p_dens, n_dens = read_density(p_file, n_file)
    p.n = p_dens
    n.n = n_dens


    en = DFT.Energy{Function}(DFT.phi, fill(0.0, length(p.n)), DFT.external, fill(0.0, length(p.n)),
        DFT.hole_coloumb, fill(0.0, length(p.n)), fill(0.0, length(n.n)))
    DFT.pre_cal_hole(en.hole, en.holeVP, p.lambda)
    DFT.pre_cal_hole(en.hole, en.holeVN, n.lambda)

    DFT.pre_cal(en.phi, en.phiV)
    DFT.pre_cal(en.ext, en.extV)
    println()
    return get_gp(p, n, en), h_box
end

using Glob
conc = 0.0002
p_files = glob("../../Documents/dft_prog/aurora_new/"*string(conc)*"M_scd-0.005/pDens*.txt")
n_files = glob("../../Documents/dft_prog/aurora_new/"*string(conc)*"M_scd-0.005/nDens*.txt")
p_files = sort(p_files, lt=(x,y)->isless(parse(Float64, match(r"pDens\_([0-9]+)", x)[1]), 
                               parse(Float64, match(r"pDens\_([0-9]+)", y)[1])))
n_files = sort(n_files, lt=(x,y)->isless(parse(Float64, match(r"nDens\_([0-9]+)", x)[1]), 
                               parse(Float64, match(r"nDens\_([0-9]+)", y)[1])))
println(parse(Float64, match(r"pDens\_([0-9]+)", p_files[9])[1]))
gps = zeros(0)
hs = zeros(0)
for i in eachindex(p_files)
    g, h = gp(p_files[i], n_files[i], conc)
    append!(gps, g)
    append!(hs, h)
end
grand = (gps .- gps[end]) ./ hs

PyPlot.clf() 
plot(hs, grand, color="green", linewidth=2.0, linestyle="-")
xlim(0, 400.0)
ylim(-1e-5, 3.5e-5)
plot(size=(200,200))
display(gcf())

savefig("../../Documents/dft_prog/aurora_new/g_"*string(conc)*".png")







function appsd()
    scd = DFT.sigma
    hs = readdlm("../../Documents/dft_prog/aurora_new/"*string(conc)*"M_scd-0.005/pDens_400.0_-0.005_"*string(conc)*".txt", 
                                                                                '\t', Float64, '\n')[:,1]
    pn = readdlm("../../Documents/dft_prog/aurora_new/"*string(conc)*"M_scd-0.005/pDens_400.0_-0.005_"*string(conc)*".txt", 
                                                                                '\t', Float64, '\n')[:,2]
    nn = readdlm("../../Documents/dft_prog/aurora_new/"*string(conc)*"M_scd-0.005/nDens_400.0_-0.005_"*string(conc)*".txt", 
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
    savefig("../../Documents/dft_prog/aurora_new/ascd_"*string(conc)*".png")
end
appsd()