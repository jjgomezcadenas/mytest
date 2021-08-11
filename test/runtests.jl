using ANema
using ATools 
using DataFrames
using Test

dir = "../data/"
conf = "n3-w-20mm-phot"
path = joinpath(dir, conf)

@testset "ANema.jl" begin
    ndf = readdf(path) 
    @test nrow(ndf) == 180079
    @test length(names(ndf)) ==  46

    ndfq =selectinterval(ndf, "q1", "q2", 1700.0, 2100.0);
    @test nrow(ndfq) == 42875

    ndfz =selectinterval(ndfq, "zstd1", "zstd2", 1.0, 16.0)
    @test nrow(ndfz) == 42806

    fz, pz = fit_profile(ndfz.zstd1, ndfz.r1, "σz", "r", "pol2")
    fpars = fz.fitpar
    fstds = fz.fitstd
    @test isapprox(fpars[1], 372.96, atol=0.1) 

    r1z = fz.g.(ndfz.zstd1)
    fg,p = fitg1(r1z - ndfz.r1, "r1 - r1z", 200, -5.0, 5.0, xgmin=-1.5, xgmax=1.5)
    @test isapprox(fg.std[1], 0.4, atol=0.1) 
    
 
end
