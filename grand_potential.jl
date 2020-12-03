
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
    grand_pot = get_gp(p, n, en)
    println("Grand potential: ", grand_pot)
    return grand_pot, h_box, sum(p.n), sum(n.n)
end

function gp_jan()
    data = reshape(readdlm("../../Downloads/jan_DFT/fcdfil"),:,3)
    hs = data[:, 1]
    pn = data[:, 2]
    nn = data[:, 3]

    h_box = 110.0 - 2.0*DFT.dz
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
    p.n = pn
    n.n = nn

    en = DFT.Energy{Function}(DFT.phi, fill(0.0, length(p.n)), DFT.external, fill(0.0, length(p.n)),
        DFT.hole_coloumb, fill(0.0, length(p.n)), fill(0.0, length(n.n)))
    DFT.pre_cal_hole(en.hole, en.holeVP, p.lambda)
    DFT.pre_cal_hole(en.hole, en.holeVN, n.lambda)

    DFT.pre_cal(en.phi, en.phiV)
    DFT.pre_cal(en.ext, en.extV)
    grand_pot = get_gp(p, n, en)
    println("Grand potential: ", grand_pot)
end

#gp("../../Documents/dft_prog/aurora_new/0.0018M_scd-0.005/pDens_110.0_-0.005_0.0018.txt",
#   "../../Documents/dft_prog/aurora_new/0.0018M_scd-0.005/nDens_110.0_-0.005_0.0018.txt", 0.0018)

using Glob
conc = 0.0018
p_files = glob("../../Documents/dft_prog/aurora_new/"*string(conc)*"M_scd-0.005/pDens*.txt")
n_files = glob("../../Documents/dft_prog/aurora_new/"*string(conc)*"M_scd-0.005/nDens*.txt")
p_files = sort(p_files, lt=(x,y)->isless(parse(Float64, match(r"pDens\_([0-9]+)", x)[1]), 
                               parse(Float64, match(r"pDens\_([0-9]+)", y)[1])))
n_files = sort(n_files, lt=(x,y)->isless(parse(Float64, match(r"nDens\_([0-9]+)", x)[1]), 
                               parse(Float64, match(r"nDens\_([0-9]+)", y)[1])))

gps = zeros(0)
hs = zeros(0)
sps = zeros(0)
sns = zeros(0)
for i in eachindex(p_files)
    g, h, sp, sn = gp(p_files[i], n_files[i], conc)
    println()
    append!(gps, g)
    append!(hs, h)
    append!(sps, sp)
    append!(sns, sn)
end
grand = (gps .- gps[end])

PyPlot.clf() 
plot(hs, grand, color="green", linewidth=2.0, linestyle="-")
xlim(0, 400.0)
#ylim(-3e-5, 3.5e-5)
plot(size=(200,200))
display(gcf())

savefig("../../Documents/dft_prog/aurora_new/g_nh_"*string(conc)*".png")



hs_jan = [20, 30, 40, 50, 60, 70, 80, 90, 100, 110]
gps_jan = [6.5195970474412959E-003, 7.8596831642891024E-003, 8.3251057675676060E-003, 8.4643027182856184E-003,
           8.4849249396834747E-003, 8.4626598129401431E-003, 8.4251515875354982E-003, 8.3825062237983111E-003,
           8.3383822360606069E-003, 8.2940580177358753E-003]
PyPlot.clf() 
plot(hs_jan, gps_jan, color="orange", linewidth=2.0, linestyle="-")
plot(hs, gps, color="green", linewidth=1.0, linestyle="-")
xlim(0, 400.0)
#ylim(-3e-5, 3.5e-5)
plot(size=(200,200))
display(gcf())

savefig("../../Documents/dft_prog/aurora_new/jan_comp_0.0018_"*string(conc)*".png")

#6.5195970474412959E-003 20
#7.8596831642891024E-003 30
#8.3251057675676060E-003 40
#8.4643027182856184E-003 50
# 60
#
#
#
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
