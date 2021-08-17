using DataFrames
using Unitful
using ATools

import Unitful:
    nm, μm, mm, cm, ns, μs, ms, ps, s

"""
    distance_lor_to_point(x1::Real,y1::Real,x2::Real,y2::Real,
                          x0::Real=0.0,y0::Real=0.0)

Compute the distance between lor defined by (x1,y1) and (x2, y2) and point (x0,y0)
"""
function distance_lor_to_point(x1::Real,y1::Real,x2::Real,y2::Real,
                               x0::Real=0.0,y0::Real=0.0)

    num = abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1))
    den = sqrt((x2 -x1)^2 + (y2 - y1)^2)
    num/den
end


"""
    function sinogramdf(lors)

Return a DataFrame from which sinograms can be computed
"""
function sinogramdf(lors)

    ldf = DataFrame(x1=[lors[i].x1 for i in 1:length(lors)],
                    y1=[lors[i].y1 for i in 1:length(lors)],
                    z1=[lors[i].z1 for i in 1:length(lors)],
                    x2=[lors[i].x2 for i in 1:length(lors)],
                    y2=[lors[i].y2 for i in 1:length(lors)],
                    z2=[lors[i].z2 for i in 1:length(lors)])

    # add z, θ and r columns 
    ldf[!, "zl"] = ldf.z1 - ldf.z2;
    ldf[!, "rl"] = distance_lor_to_point.(ldf.x1, ldf.y1, ldf.x2, ldf.y2, 
                                      zeros(nrow(ldf)), zeros(nrow(ldf)) )
    ldf[!, "tl"] = atan.(ldf.x1, ldf.y1)  
    return ldf
end   


"""
    ztsinogram(sdf; nbinz=3, zmin=-350.0, zmax=350.0,
                    nbint=3, tmin=-Float64(π), tmax=Float64(π))

Return a dictionary binning in z and θ for sinogram calculations 
"""
function ztsinogram(sdf::DataFrame; 
                    nbinz::Integer=3, zmin::Float64=-350.0, zmax::Float64=350.0,
                    nbint::Integer=3, tmin::Float64=-Float64(π), tmax::Float64=Float64(π))

    if nbinz == 1
        RL = [sdf]
    else
        # bins in z 
        hz  = hist1d(sdf.zl, nbinz, zmin, zmax)
        hze = hz.edges[1]
        RL = [select_by_column_value_interval(sdf, "zl", 
        hze[i], hze[i+1]) for i in 1:length(hze)-1]
    end

    # select in theta bins
    dRL = Dict()
    for (i, rl) in enumerate(RL)
        if nbint == 1
            RT =[rl]
        else
            ht  = hist1d(rl.zl, nbint, tmin, tmax)
            hte = ht.edges[1]
            RT = [select_by_column_value_interval(sdf, 
                  "tl", hte[i], hte[i+1]) for i in 1:length(hte)-1]
        end
        dRL[i] = RT
    end
    return dRL
end


"""
    thetasgrm(drl; nbint, nbin=20, rmin=0.0, rmax=75.0)

Take the dictionary drl and return projected θ  sinograms 
that is histograms in θ for each theta in the dictionary 
"""
function thetasgrm(drl; nbint, nbin=20, rmin=0.0, rmax=75.0)
    nbinz = 1
    PR = []
    HR = []
    for i in 1:nbint
        sng = drl[1][i]
        h,p = hist1d(sng.rl, "rl", nbin, rmin, rmax)
        push!(PR, p)
        push!(HR, h)
    end
    return HR, PR
end


"""
    zsgrm(drl; nbinz, nbin=20, rmin=0.0, rmax=75.0)
Take the dictionary drl and return projected z sinograms 
that is histograms in r for each z in the dictionary 
"""
function zsgrm(drl; nbinz, nbin=20, rmin=0.0, rmax=75.0)
    nbint = 1
    PR = []
    HR = []
    for i in 1:nbinz
        sng = drl[i][1]
        h,p = hist1d(sng.rl, "rl", nbin, rmin, rmax)
        push!(PR, p)
        push!(HR, h)
    end
    return HR, PR
end

