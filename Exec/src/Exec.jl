
module Exec
using ArgParse
include("interpreter.jl")
#fooifier_path() = joinpath(artifact"fooifier", "bin", "fooifier" * (Sys.iswindows() ? ".exe" : ""))

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "filename"
            help = "filename argument"
            required = true
    end
    return parse_args(s)
end

function real_main()
    parsed_args = parse_commandline()
    filename = parsed_args["filename"]
    open(filename, "r") do openedFile
        input = read(openedFile, String)
        lines = split(input, '\n')
        interpreter = create_model(input)
        println(interpreter.output)
    end
end

function julia_main()
    try
        real_main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    real_main()
end

end # module