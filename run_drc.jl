using YAML
using ArgParse
using DifferentialEquations
using DiffEqSensitivity
using ForwardDiff
using NPZ
using GRUtils
using Sobol
using Statistics


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "input_path"
            help = "path of input file"
            arg_type = String
            required = true
        "--output_path"
            help = "path of directory for output result"
            arg_type = String
            default = "output"
        "--delsa"
            help = "if using delsa mode for uncertain parameters"
            action = :store_true
        "--delsa_sample"
            help = "number of samples used in delsa"
            arg_type = Int
            default = 200
    end

    return parse_args(s)
end


function sep_coef_species(num_name)
    """
    Separate stoichiometry and species
    """
    coef = ""
    species = ""
    for (i, c) in enumerate(num_name)
        if (isdigit(c))
            coef = coef * c
        else
            species = num_name[i:end]
            break
        end
    end
    if coef == ""
        coef = '1'
    end
    return parse(Int64, coef), species
end


function get_species_coef(ads_species, gas_species, reactions)
    """
    get coefficient matrixs of lhs and rhs of reactions
    for adsorbed species and gas species
    """
    n_ads, n_gas, n_rxn = size(ads_species, 1), size(gas_species, 1), size(reactions, 1)
    ads_map = Dict([(k, v) for (v, k) in enumerate(ads_species)])
    gas_map = Dict([(k, v) for (v, k) in enumerate(gas_species)])
    lhs_ads_coef = zeros(Int64, n_rxn, n_ads)
    rhs_ads_coef = zeros(Int64, n_rxn, n_ads)
    lhs_gas_coef = zeros(Int64, n_rxn, n_gas)
    rhs_gas_coef = zeros(Int64, n_rxn, n_gas)
    for (ir, rxn) in enumerate(reactions)
        lhs, rhs = split(rxn, "<->")
        for lhs_s in split(lhs, '+')
            lhs_coef, lhs_species = sep_coef_species(strip(lhs_s))
            if (endswith(lhs_species, "gas"))
                lhs_gas_coef[ir, gas_map[lhs_species]] += lhs_coef
            else
                lhs_ads_coef[ir, ads_map[lhs_species]] += lhs_coef
            end
        end
                    
        for rhs_s in split(rhs, '+')
            rhs_coef, rhs_species = sep_coef_species(strip(rhs_s))
            if (endswith(rhs_species, "gas"))
                rhs_gas_coef[ir, gas_map[rhs_species]] += rhs_coef
            else
                rhs_ads_coef[ir, ads_map[rhs_species]] += rhs_coef
            end
        end
    end
    return lhs_ads_coef, rhs_ads_coef, lhs_gas_coef, rhs_gas_coef
end


function get_gas_factor(lhs_gas_coef, rhs_gas_coef, gas_pressures)
    """
    get gas factor
    """
    if (size(gas_species, 1) == 0)
        lhs_gas_factors = ones(size(lhs_gas_coef, 1))  # n reactions
        rhs_gas_factors = ones(size(lhs_gas_coef, 1))
    else
        gas_pressures = reshape(gas_pressures, 1, size(gas_pressures, 1))
        gas_pressures_rep = repeat(gas_pressures, inner=(size(lhs_gas_coef, 1), 1))
        lhs_gas_factors = prod(gas_pressures_rep .^ lhs_gas_coef, dims=2)
        rhs_gas_factors = prod(gas_pressures_rep .^ rhs_gas_coef, dims=2)
    end
    return lhs_gas_factors, rhs_gas_factors
end


function calc_rate(y, kf)
    global K, lhs_coef, rhs_coef, lhs_gas_factor, rhs_gas_factor
    kr = kf ./ K
    nr, ns = size(kf, 1), size(y, 1)
    y_tile = repeat(reshape(y, 1, ns), nr, 1)
    lhs_factor = prod(y_tile .^ lhs_coef, dims=2) .* lhs_gas_factor
    rhs_factor = prod(y_tile .^ rhs_coef, dims=2) .* rhs_gas_factor
    rates = lhs_factor .* kf - rhs_factor .* kr
    return rates
end


function calc_overall_rate(y, kf)
    r = calc_rate(y, kf)
    global overall_rate_mask
    out_r = sum(r .* overall_rate_mask)
end


function state_eqn(y, kf, t)
    global lhs_coef, rhs_coef
    r = calc_rate(y, kf)
    nr, ns = size(kf, 1), size(y, 1)
    r_tile = repeat(reshape(r, nr, 1), 1, ns)
    species_rates = sum(r_tile .* (-lhs_coef .+ rhs_coef), dims=1)
    return reshape(species_rates, ns)
end


function rate_wrapper(p)
    p = exp.(p)
    global sts, ste, snt, prob
    _prob = remake(prob, p=p)
    theta_sol = Array(solve(_prob, Kvaerno5(), saveat=LinRange(sts, ste, snt), sensealg=ForwardDiffSensitivity()))'
    p_matrix = reshape(p, 1, size(p, 1))
    p_repeat = repeat(p_matrix, snt, 1)
    r = Array{Real, 1}(undef, snt)
    for i in 1:snt
        r[i] = calc_overall_rate(theta_sol[i, :], p_repeat[i, :])
    end
    return log.(r)
end

function plot_delsa_sen(output, array, type)
    n_params = size(array, 3)  # n_sample, snt, nr
    for i_p in 1:n_params
        fig = Figure()
        avg = round(mean(array[:, end, i_p]), digits=3)
        histogram(array[:, end, i_p], nbins=15, color=0xfc4e03, xlabel=type, ylabel="Frequency", title="$(type)-$(i_p) mean: $(avg)")
        fname = type == "Sensitivity" ? "Sensit" : "DRC"
        savefig(joinpath(output, "delsa-$(fname)-$(i_p).png"))
    end
end


flags = parse_commandline()

if (~ isdir(flags["output_path"]))
    mkdir(flags["output_path"])
end

data = YAML.load_file(flags["input_path"])

y0_dict = data["C0"]
ts, te, nt, sts, ste, snt = Array(data["Time"])
nt, snt = Int64(nt), Int64(snt)
overall_rate_mask = Array(data["Rate_mask"])

ads_species = collect(keys(y0_dict))
y0 = Array{Float64}(undef, size(ads_species, 1))
for (ind, spe) in enumerate(ads_species)
    y0[ind] = y0_dict[spe]
end

press_dict = haskey(data, "Pressure") ? data["Pressure"] : Dict()
gas_species = collect(keys(press_dict))
lhs_coef, rhs_coef, lhs_gas_coef, rhs_gas_coef = get_species_coef(ads_species, gas_species, data["Reactions"])
gas_pressures = Array([press_dict[x] for x in gas_species])
lhs_gas_factor, rhs_gas_factor = get_gas_factor(lhs_gas_coef, rhs_gas_coef, gas_pressures)


if !flags["delsa"]
    kf0 = Array(data["kf"])
    K = Array(data["K"])

    tspan = (ts, te)
    prob = ODEProblem(state_eqn, y0, tspan, kf0)
    ode_sol = solve(prob, Kvaerno5(), saveat=LinRange(ts, te, nt), atol=1e-8, rtol=1e-8)

    drdp = ForwardDiff.jacobian(rate_wrapper, log.(kf0))

    ode_sol_path = joinpath(flags["output_path"], "ode-sol.npy")
    npzwrite(ode_sol_path, Array(ode_sol)')  # nt, ns

    drc_sol_path = joinpath(flags["output_path"], "drc-res.npy")
    npzwrite(drc_sol_path, drdp)  # snt, nr

    fig = Figure()
    plot(LinRange(ts, te, nt), Array(ode_sol)', xlabel="Time (s)", ylabel = "Coverage", linewidth=3, dpi=300)
    legend(ads_species..., location=12)
    savefig(joinpath(flags["output_path"], "Coverage-time.png"))

    fig = Figure()
    plot(LinRange(sts, ste, snt), drdp, xlabel="Time (s)", ylabel = "Degree of rate control", linewidth=3, dpi=300)
    drc_legend = ["DRC $i" for i in 1:size(drdp, 2)]
    legend(drc_legend..., location=12)
    savefig(joinpath(flags["output_path"], "DRC-time.png"))
else
    kf0_lb = Array(data["kf_lb"])
    kf0_ub = Array(data["kf_ub"])
    K_lb = Array(data["K_lb"])
    K_ub = Array(data["K_ub"])

    tspan = (ts, te)
    prob = ODEProblem(state_eqn, y0, tspan, kf0_lb)
    n_sample = flags["delsa_sample"]
    drdp_sample = Array{Float64, 3}(undef, n_sample, snt, size(kf0_lb, 1))

    s_id = 1
    sobol_lnK_eq = SobolSeq(log.(K_lb), log.(K_ub))
    for lnk in SobolSeq(log.(kf0_lb), log.(kf0_ub))
        global s_id, K
        K = exp.(next!(sobol_lnK_eq))        
        drdp_sample[s_id, :, :] = ForwardDiff.jacobian(rate_wrapper, lnk)
        if (s_id == n_sample)
            break
        end
        s_id += 1
    end

    delsa_res_path = joinpath(flags["output_path"], "delsa-drc.npy")
    npzwrite(delsa_res_path, drdp_sample)  # n_sample, snt, nr

    var_lnk = ((log.(kf0_ub) .- log.(kf0_lb)) .^ 2) ./ 12  # nr
    var_lnk_tile = repeat(reshape(var_lnk, 1, 1, size(kf0_lb, 1)), inner=(n_sample, snt, 1))
    var_single = (drdp_sample.^2) .* var_lnk_tile
    var_sum = sum(var_single, dims=3)
    var_sum = repeat(var_sum, inner=(1, 1, size(kf0_lb, 1)))
    sensit = var_single ./ var_sum
    
    sensit_path = joinpath(flags["output_path"], "delsa-sensit.npy")
    npzwrite(sensit_path, sensit)  # n_sample, snt, nr

    plot_delsa_sen(flags["output_path"], drdp_sample, "Degree of rate control")
    plot_delsa_sen(flags["output_path"], sensit, "Sensitivity")
end