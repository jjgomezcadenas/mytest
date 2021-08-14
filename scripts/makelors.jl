using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using DataFrames
using CSV
using ArgParse
using Logging
using Printf
using ATools
include("../src/ANema.jl")

import ANema: getlorpath
   

logger = SimpleLogger(stdout, Logging.Warn)
old_logger = global_logger(logger)


function makelors(args)

	rdir   = args["dir"]
	conf   = args["conf"]
	pathf = ANema.pathfinder(rdir, conf)
	println("pathfinder tuple =", pathf)

	n3p = DataFrame(CSV.File(pathf.phipath))
	n3z = DataFrame(CSV.File(pathf.zpath))

	n3mp = ANema.setunits(n3p)
	n3mz = ANema.setunits(n3z)

	println("phi data frame has length of ", nrow(n3mp))
	println("z data frame has length of ", nrow(n3mz))
	println("phi data frame names ", names(n3mp))
	println("z data frame names ", names(n3mz))

	lors_first_true = ANema.dftolor(n3mz, ANema.dtfirst, ANema.postrue)
	ANema.write_lors_hdf5(getlorpath(pathf.lordir,pathf.zlp, "first-true"), lors_first_true)
	lors_minimum_reco = ANema.dftolor(n3mz, ANema.dtminimum, ANema.posreco)
	ANema.write_lors_hdf5(getlorpath(pathf.lordir,pathf.zlp, "minimum-reco"), lors_minimum_reco)
	lors_minimum_recall = ANema.dftolor(n3mz, ANema.dtminimum, ANema.posrecall)
	ANema.write_lors_hdf5(getlorpath(pathf.lordir,pathf.zlp, "minimum-recall"), lors_minimum_recall)
	lors_average_recall = ANema.dftolor(n3mz, ANema.dtaverage, ANema.posrecall)
	ANema.write_lors_hdf5(getlorpath(pathf.lordir,pathf.zlp, "average-recall"), lors_average_recall)

	lors_first_true = ANema.dftolor(n3mp, ANema.dtfirst, ANema.postrue)
	ANema.write_lors_hdf5(getlorpath(pathf.lordir,pathf.philp, "first-true"), lors_first_true)
	lors_minimum_reco = ANema.dftolor(n3mp, ANema.dtminimum, ANema.posreco)
	ANema.write_lors_hdf5(getlorpath(pathf.lordir,pathf.philp, "minimum-reco"), lors_minimum_reco)
	lors_minimum_recall = ANema.dftolor(n3mp, ANema.dtminimum, ANema.posrecall)
	ANema.write_lors_hdf5(getlorpath(pathf.lordir,pathf.philp, "minimum-recall"), lors_minimum_recall)
	lors_average_recall = ANema.dftolor(n3mp, ANema.dtaverage, ANema.posrecall)
	ANema.write_lors_hdf5(getlorpath(pathf.lordir,pathf.philp, "average-recall"), lors_average_recall)

end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dir", "-d"
            help = "directory with nema data"
            arg_type = String
            default = "../data"
		"--conf", "-c"
            help = "configuration"
            arg_type = String
    end

    return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Running makelors with arguments", parsed_args)
	makelors(parsed_args)
end

@time main()
