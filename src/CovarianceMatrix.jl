module CovarianceMatrix

# External Packages
using TOML
using FITSIO
using OLUtils

# Internal Packages

# Exports
export create_covariance_matrix

struct Filter
    name::AbstractString
    CW::Float64
    σZP::Float64
    σFilter::Float64
end

struct Survey
    name::AbstractString
    filters::Vector{Filter}
end

function load_survey(survey_path::AbstractString)
    survey_name = splitdir(splitext(survey_path)[1])[end]
    @info "Found $survey_name"
    survey = open(survey_path, "r") do io
        lines = [line for line in readlines(io) if !occursin("#", line)][2:end] # Remove comments and ignore header
        filters = []
        for line in lines
            name, CW, σZP, σFilter = split(line, ",")
            filter = Filter("$(survey_name)_$(uppercase(strip(name)))", parse(Float64, strip(CW)), parse(Float64, strip(σZP)), parse(Float64, strip(σFilter)))
            push!(filters, filter)
        end
        return Survey(survey_name, filters)
    end
end

function create_covariance_matrix(toml::Dict, verbose::Bool)
    setup_global!(toml, verbose)
    config = toml["global"]
    
    @debug "Base path: $(config["base_path"])"
    @debug "Output path: $(config["output_path"])"
    
    # First load in filters, if any
    surveys = []
    if "survey" in keys(toml)
        survey_dir = toml["survey"]["path"]
        if !isabspath(survey_dir)
            survey_dir = joinpath(config["base_path"], survey_dir)
        end
        survey_dir = abspath(survey_dir)
        for file in readdir(survey_dir, join=true)
            if splitext(file)[end] == ".survey"
                survey = load_survey(file)
                push!(surveys, survey)
            end
        end
    end
    @info "Found $(length(surveys)) surveys"
    sort!(surveys, lt = (s1, s2) -> isless(s1.name, s2.name))
    @debug "Surveys: $([s.name for s in surveys])" 
    all_filters = []
    for survey in surveys
        all_filters = vcat(all_filters, survey.filters)
    end
    sort!(all_filters, lt = (f1, f2) -> isless(f1.name, f2.name))
    num_filters = length(all_filters)
    @info "Found $num_filters total filters"
    @debug "Filters: $([f.name for f in all_filters])"

    # Next generate an empty covariance matrix
    cov_size = 2 * num_filters
    @info "Generating a $(cov_size)x$(cov_size) covariance matrix"
    cov = zeros(cov_size, cov_size)

    # Calculate survey uncertainties
    @info "Calculating survey uncertainties"
    slope = 0.005
    waveStart = 300
    waveEnd = 1000
    central = 555.6
    value = slope / (waveEnd - waveStart)
    for (i, f1) in enumerate(all_filters)
        for (j, f2) in enumerate(all_filters)
            if i == j
                # ZP uncertainties
                cov[i, i] += f1.σZP * f1.σZP
                # Filter uncertainties
                cov[num_filters + i, num_filters + i] += f1.σFilter * f1.σFilter
            end
            # Covariant uncertainty
            e1 = value * (f1.CW - central)
            e2 = value * (f2.CW - central)
            cov[i, j] += e1 * e2
        end
    end

    # Save covariance matrix
    output = get(toml, "output", Dict("name" => "DES.fits"))
    save_file = joinpath(config["output_path"], output["name"])
    @info "Saving covariance matrix to $(output["name"])"
    FITS(save_file, "w") do io
        write(io, cov)
    end

    header = FITS(save_file, "r") do io
        return read_header(io[1])
    end

    for (i, filter) in enumerate(all_filters)
        header["$(filter.name)"] = i-1 
    end
    header["dZP"] = string([0, num_filters - 1])
    header["dFilter"] = string([num_filters, 2 * num_filters])
    FITS(save_file, "w") do io
        write(io, cov; header=header)
    end

end

function create_covariance_matrix(toml_path::AbstractString, verbose::Bool)
    toml = TOML.parsefile(toml_path)
    if !("global" in keys(toml))
        toml["global"] = Dict()
    end
    toml["global"]["toml_path"] = dirname(abspath(toml_path))
    return create_covariance_matrix(toml, verbose)
end

end # module
