using ATools 
using DataFrames
using Glob
using CSV
using Plots



"""
    readdf(dir)

Reads all csv files found in dir and returns a single DataFrame.
"""
function readdf(dir)
    files = glob("*.csv",dir)
    dfs =[DataFrame(CSV.File(file)) for file in files]
    evtdf=vcat(dfs...)
end


function writemdf(dir, file, df)
    path = string(dir,"/", file)
    CSV.write(path, df)
end


function pathfinder(rootdir, conf)
    mdf = string(conf, "-mdf")
    lor = string(conf, "-lor")
    phif = string("mdf-phistd-",conf,".csv")
    zf = string("mdf-zstd-",conf,".csv")

    mdfdir = joinpath(rootdir, mdf)
    phipath = joinpath(mdfdir, phif)
    zpath = joinpath(mdfdir, zf)
    
    lordir = joinpath(rootdir, lor)

    if isdir(lordir) == false
        mkdir(lordir)
    end

    philp = string("lor-phistd-",conf)
    zlp = string("lor-zstd-",conf)
    return (phipath= phipath, zpath = zpath, lordir=lordir, philp = philp, zlp = zlp)

end


function getlorpath(lordir, prefix, posfix="first-true")
    fullpath = string(prefix,"-",posfix,".h5")
    joinpath(lordir, fullpath)
end

"""
    selectinterval(df, column1, column2, xmin, xmax)

Select interval over 2 hemisphere columns (e.g, q1, q2) for df
"""
function selectinterval(df, column1, column2, xmin, xmax)
    e1 = select_by_column_value_interval(df, column1, xmin,xmax)
    select_by_column_value_interval(e1, column2, xmin,xmax)
end

function pout(p,filename)
    if filename != ""
        png(p,filename)
    end
end


function plot_and_save(p1::Plots.Plot, tit::String="", filename::String="", legend=true)
    if tit != ""
        p = plot(p1, layout= (1, 1), title=tit, legend=legend, fmt = :png, size = (1000, 400),
            left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    else
        p = plot(p1, layout= (1, 1), legend=legend, fmt = :png, size = (1000, 400),
            left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    end
    pout(p,filename)
    return p
end


function plot_and_save(p1::Plots.Plot,p2::Plots.Plot, 
                       tit::String="", filename::String="", legend=true)
    if tit != ""
        p = plot(p1, p2, layout= (1, 2), title=tit, legend=legend, fmt = :png, size = (1000, 400),
            left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    else
        p = plot(p1, p2, layout= (1, 2), legend=legend, fmt = :png, size = (1000, 400),
            left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    end
    pout(p,filename)
    return p
end


function plot_and_save(p1::Plots.Plot,p2::Plots.Plot, p3::Plots.Plot, tit::String="", 
    filename::String="", legend=true)
    if tit != ""
        p = plot(p1, p2, p3, layout= (1, 3), title=tit, legend=legend, fmt = :png, size = (1000, 400),
            left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    else
        p = plot(p1, p2, p3, layout= (1, 3), legend=legend, fmt = :png, size = (1000, 400),
            left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    end
    pout(p,filename)
    return p
end


# function plot_and_save(p1,p2,p3,p4, tit="", filename="")
#     if tit != ""
#         p = plot(p1, p2, p3, p4, layout= (2, 2), title=tit, legend=false, fmt = :png, size = (1000, 400),
#             left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
#     else
#         p = plot(p1, p2, p3, p4, layout= (2, 2), legend=false, fmt = :png, size = (1000, 400),
#             left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
#     end
#     pout(p,filename)
#     return p
# end


# function condplot(p1, tit, filename, save)
#     if save
#         plot_and_save(p1, tit, filename)
#     else
#         plot(p1, layout= (1, 1), legend=false, fmt = :png, size = (1000, 400),
#             left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
#     end
# end


function condplot(p1, p2, tit, filename, save)
    if save
        plot_and_save(p1,p2, tit, filename)
    else
        plot(p1, p2, layout= (1, 2), legend=false, fmt = :png, size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    end
end


# function condplot(p1, p2, p3, tit, filename, save)
#     if save
#         plot_and_save(p1,p2, p3, tit, filename)
#     else
#         plot(p1, p2, p3, layout= (1, 3), legend=false, fmt = :png, size = (1000, 400),
#         left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
#     end
# end


function condplot(p1, p2, p3, p4, tit, filename, save)
    if save
        plot_and_save(p1,p2, p3, p4, tit, filename)
    else
        plot(p1, p2, p3, p4, layout= (2, 2), legend=false, fmt = :png, size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    end
end


function q1vsq2(df, nbins=100, qmin=50.0, qmax=5000.0; tit="", filename="", save=false)
    h1,p1 = hist2d(df.q1, df.q2, nbins, "q1 (pes)","q2 (pes)", qmin, qmax, qmin, qmax)
    h2,p2 = hist1d(df.q1, "q1 (pes)", nbins, qmin, qmax);
    condplot(p1, p2, tit, filename, save)

end


function r1q1(df, nbins=100, qmin=50.0, qmax=5000.0, rmin=300.0, rmax=400.0; tit="", filename="", save=false)
    h1,p1 = hist2d(df.q1, df.r1, nbins, "q1 (pes) ","r1 (mm)",qmin, qmax, rmin,rmax)
    h2,p2 = hist1d(df.r1, "r1", nbins, rmin,rmax)

    condplot(p1, p2, tit, filename, save)
end


function zstd(df, nbins=100, zmin=0.0, zmax=40.0,rmin=300.0, rmax=400.0; tit="", filename="", save=false)
    h1,p1 = hist2d(df.zstd1, df.r1, nbins, "σz (mm) ","1 (mm)",zmin, zmax, rmin, rmax)
    h2,p2 = hist1d(df.zstd1, "σz (mm)", nbins, zmin, zmax)
    condplot(p1, p2, tit, filename, save)
    
end


function phistd(df, nbins=100,phimin=0.0, phimax=0.1,rmin=300.0, rmax=400.0; tit="", filename="", save=false)
    h1,p1 = hist2d(df.phistd1, df.r1, nbins, "σϕ (mm) ","r1 (mm)",phimin, phimax,rmin, rmax)
    h2,p2 = hist1d(df.phistd1, "σϕ (mm)", 100, phimin, phimax)
    condplot(p1, p2, tit, filename, save)
end


function plotreso(r1t, r1x, tx1, ty1, xmin, xmax, nbins=150; tit="", filename="", save=false)
    h1,p1 = hist2d(r1x, r1t, nbins, tx1, ty1)
    h2,p2 = hist1d(r1t - r1x, tit, nbins, xmin, xmax)
    condplot(p1, p2, tit, filename, save)
end



