function read_dens(h, p_file, n_file, pdens, ndens)
    h_vals = readdlm(p_file, '\t', Float64, '\n')[:,1]
    p_dens = readdlm(p_file, '\t', Float64, '\n')[:,2]
    n_dens = readdlm(n_file, '\t', Float64, '\n')[:,2]

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