module CovarianceMatrix

# External Packages
using TOML
using LoggingExtras
using FITSIO

# Internal Packages

# Exports
export create_covariance_matrix

function setup_global_config!(toml::Dict)
    config = get(toml, "global", Dict())
    # Base path is where everything relative will be relative to
    # Defaults to the directory containing the toml path
    # Can be relative (to the toml path) or absolute
    base_path = get(config, "base_path", nothing)
    if isnothing(base_path)
        base_path = dirname(toml["toml_path"])
    elseif !isabspath(base_path)
        base_path = joinpath(dirname(toml["toml_path"]), base_path)
    end
    base_path = abspath(base_path)
    config["base_path"] = base_path
    # Output path is where all output (figures, files, etc...) will be placed
    # Defaults to base_path / Output
    # Can be relative (to base_path) or absolute
    output_path = get(config, "output_path", nothing)
    if isnothing(output_path)
        output_path = joinpath(base_path, "Output")
    elseif !isabspath(output_path)
        output_path = joinpath(base_path, output_path)
    end
    config["output_path"] = abspath(output_path)
    # Logging sets whether or not to setup and use Supernovae's logging
    logging = get(config, "logging", true)
    config["logging"] = logging
    # Log file is the name of the log file. This will only work if logging is true
    # Can only be relative to output_path
    # Defaults to log.txt
    log_file = get(config, "log_file", nothing)
    if logging
        if isnothing(log_file)
            log_file = "log.txt"
        end
        log_file = abspath(joinpath(output_path, log_file))
    end
    if !logging & !isnothing(log_file)
        @warn "Logging set to false, so log file $log_file will not be written. Please add `logger=true` to your [ global ] config"
    end
    config["log_file"] = log_file
    toml["global"] = config
end

function setup_logger(log_file::AbstractString, verbose::Bool)
    if verbose
        level = Logging.Debug
    else
        level = Logging.Info
    end
    function fmt(io, args)
        if args.level == Logging.Error
            color = :red
            bold = true
        elseif args.level == Logging.Warn
            color = :yellow
            bold = true
        elseif args.level == Logging.Info
            color = :cyan
            bold = false
        else
            color = :white
            bold = false
        end
        printstyled(io, args._module, " | ", "[", args.level, "] ", args.message, "\n"; color = color, bold = bold)
    end
    logger = TeeLogger(
        MinLevelLogger(FormatLogger(fmt, open(log_file, "w")), level),
        MinLevelLogger(FormatLogger(fmt, stdout), level)
    )
    global_logger(logger)
    @info "Logging to $log_file"
end

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
    setup_global_config!(toml)
    config = toml["global"]
    
    # Ensure all path's exist
    if !isdir(config["base_path"])
        mkpath(config["base_path"])
    end
    if !isdir(config["output_path"])
        mkpath(config["output_path"])
    end
    # Optionally set up logging
    if config["logging"]
        setup_logger(config["log_file"], verbose)
    end
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
    toml["toml_path"] = abspath(toml_path)
    return create_covariance_matrix(toml, verbose)
end

end # module
