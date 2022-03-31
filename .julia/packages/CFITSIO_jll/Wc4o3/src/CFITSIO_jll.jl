# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule CFITSIO_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("CFITSIO")
JLLWrappers.@generate_main_file("CFITSIO", UUID("b3e40c51-02ae-5482-8a39-3ace5868dcf4"))
end  # module CFITSIO_jll
