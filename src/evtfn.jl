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


"""
    selectinterval(df, column1, column2, xmin, xmax)

Select interval over 2 hemisphere columns (e.g, q1, q2) for df
"""
function selectinterval(df, column1, column2, xmin, xmax)
    e1 = select_by_column_value_interval(df, column1, xmin,xmax)
    select_by_column_value_interval(e1, column2, xmin,xmax)
end



