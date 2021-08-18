using DataFrames
using Unitful
using StatsBase
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


function zxdf(sdf::DataFrame;  nbinz::Integer=3, zmin::Float64=-350.0, zmax::Float64=350.0)

    if nbinz == 1
        ZDF = [sdf]
    else
        # bins in z 
        hz  = hist1d(sdf.zl, nbinz, zmin, zmax)
        hze = hz.edges[1]
        ZDF = [select_by_column_value_interval(sdf, "zl", 
        hze[i], hze[i+1]) for i in 1:length(hze)-1]
    end
    return ZDF
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
    thetasgrm(df::DataFrame; ntproj::Integer, nbinr=20, rmin=0.0, rmax=75.0)

"""
function thetasgrm(df::DataFrame; ntproj::Integer, nbinr=20, rmin=0.0, rmax=75.0)
    
    ht  = hist1d(df.tl, ntproj, -Float64(π), Float64(π))
    hte = ht.edges[1]
    RT = [select_by_column_value_interval(df, 
                                         "tl", 
                                         hte[i], hte[i+1]) for i in 1:length(hte)-1]
    rH = [hist1d(dfx.rl, nbinr, rmin, rmax) for dfx in RT]
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

"""
    stp(hrl::Histogram; i0::Integer = 5, i1::Integer = 15)

Return Scatters, Trues and Prompts from an Histogram of nema4 radial distance
to source
    i0 :  defines the bin of background (left)
    i1 :  defines the bin of background (right)
    
Background is computed as linear interpolation between i0 and i1
"""
function stp(hrl::Histogram; i0::Integer = 5, i1::Integer = 15)
    w = Float64.(hrl.weights)
    e = hrl.edges[1]
    w0 = w[i0]                        # contents of bin
    w1 = w[i1]
    wm, im, em = find_max_xy(w, e)
    spb = w[i0:i1]
    fxy = gline2p(i0,w0,i1, w1)
    bkg = [fxy(i) for i in i0:i1]
    sgn = spb - bkg
    sgn = [s >0 ? ceil(s) : 0.0 for s in sgn]
    T = sum(sgn)
    P = sum(w)
    S = P - T

    @info "i0 = %d, im= %d, i1 = %d" i0 im i1
    @info "w0 = %5.1f, wm= %5.1f, w1 =%5.1f" w0 wm w1

    return S,T,P
end

"""
    pok1(λ::Float64)
    Poisson probability for k=0, 1 (k! = 1)

"""
pok1(λ::Float64) = exp(-λ) * (1.0 + λ)


"""
    pdt(rkcps::Float64, wmus::Float64)

Given an activity rkcps (in kcps) and a dead time wmus (in mus), the average number
of events in that window is kcps/(1/wmus) = kcps * wmus. 
"""
function pdt(rkcps::Float64, wmus::Float64)
    λ = rkcps * wmus * 1e-3  # λ = rate/(1/w) = rate * w where rate in kcps, w in mus
    return pok1(λ)
end

"""
    cnec(n_kcps::Float64, 
         eff_gp::Float64, 
         eff_sp::Float64, 
         SOT::Float64       =1.6, 
         eff_cwns::Float64  = 1e-9,
        deadt_mus::Float64 = 0.6)

Computes NEC curve 
"""
function cnec(n_kcps::Float64, 
    eff_gp::Float64, 
    eff_sp::Float64, 
    SOT::Float64       =1.6, 
    eff_cwns::Float64  = 1e-9,
    deadt_mus::Float64 = 0.6)

# n_kcps is the number of counts in kcps 
# eff_gp is the efficiency for good prompts
# eff_sp is the efficiency for single prompts
# SOT = S/T, the ratio of plots to trues 
# eff_cwns is the coincidence window in ns. Since rates come in kcps multiply by 10^3

eff_w =  eff_cwns *1e+3
nT = n_kcps * eff_gp               # number of good coincidences in scanner
nS = nT * SOT                      # number of scatters in scanner
n1 = n_kcps * eff_sp                # number of singles in scanner
nR = 2 * n1^2 * eff_w              # number of randoms
ntot = nT + nS + nR

prob = pdt(ntot, deadt_mus)
nT   *= prob
nS   *= prob
#nR   *= prob
ntot *= prob

nec = nT / (1.0 + nS/nT + nR/nT)   # nec


return nT, nS, nR, ntot, nec, prob
end