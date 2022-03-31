module CovarianceMatrix

# External Packages
using TOML
using LoggingExtras

# Internal Packages

# Exports

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

end

function create_covariance_matrix(toml_path::AbstractString, verbose::Bool)
    toml = TOML.parsefile(toml_path)
    toml["toml_path"] = abspath(toml_path)
    return create_covariance_matrix(toml, verbose)
end

end # module
