## --- XRF data file cleaned up
using DelimitedFiles
data = readdlm("xrf_data.csv", ',')

headers = findall(data[:,1] .== "File #")
alldata = findall(isa.(data[:,1], Number))

allcolumns = unique(data[headers,:])
data_parsed = Array{Any}(undef, length(alldata)+1, length(allcolumns))
data_parsed .= NaN
data_parsed[1,:] = allcolumns
parsedrow = 2
for row in alldata
    headerrow = headers[findlast(x->x<row, headers)]
    for col in 1:size(data, 2)
        parsedcolumn = findfirst(x->x==data[headerrow, col], allcolumns)
        value = data[row,col]
        if value == "< LOD"
            value = 0.0
        end
        data_parsed[parsedrow, parsedcolumn] = value
    end
    parsedrow += 1
end

# ds = elementify(data_parsed, importas=:Tuple)
open("xrf_data_cleaned.csv", "w") do io
    writedlm(io, data_parsed, ',')
end

## --- Individual element scatter plots to investigate relationships
using StatGeochem, Plots
ds = importdataset("xrf_data_cleaned.csv", importas=:Tuple)

scatter(ds.Cl, ds.S,
    xerror = ds.Cl_Err,
    yerror = ds.S_Err,
    xlabel = "Cl",
    ylabel = "S",
    framestyle = :box,
    label = "",
 )

 
h = scatter(ds.SiO2, ds.Cl,
    xerror = ds.SiO2_Err,
    yerror = ds.Cl_Err,
    xlabel = "SiO2",
    ylabel = "Cl",
    framestyle = :box,
    label = "",
)
savefig(h, "scatter_sio2_cl.pdf")

h = scatter(ds.SiO2, ds.S,
    xerror = ds.SiO2_Err,
    yerror = ds.S_Err,
    xlabel = "SiO2",
    ylabel = "S",
    framestyle = :box,
    label = "",
)
savefig(h, "scatter_sio2_s.pdf")

scatter(ds.SiO2, ds.Al2O3,
    xerror = ds.SiO2_Err,
    yerror = ds.Al2O3_Err,
    xlabel = "SiO2",
    ylabel = "Al2O3",
    framestyle = :box,
    label = "",
)

h = scatter(ds.SiO2, ds.Ca,
    xerror = ds.SiO2_Err,
    yerror = ds.Ca_Err,
    xlabel = "SiO2",
    ylabel = "Ca",
    framestyle = :box,
    label = "",
)
savefig(h, "scatter_sio2_ca.pdf")

scatter(ds.SiO2, ds.K2O,
    xerror = ds.SiO2_Err,
    yerror = ds.K2O_Err,
    xlabel = "SiO2",
    ylabel = "K2O",
    framestyle = :box,
    label = "",
)

h = scatter(ds.Fe, ds.Cl,
    xerror = ds.Fe_Err,
    yerror = ds.Cl_Err,
    xlabel = "Fe",
    ylabel = "Cl",
    framestyle = :box,
    label = "",
)
savefig(h, "scatter_fe_cl.pdf")

h = scatter(ds.Fe, ds.S,
    xerror = ds.Fe_Err,
    yerror = ds.S_Err,
    xlabel = "Fe",
    ylabel = "S",
    framestyle = :box,
    label = "",
)
savefig(h, "scatter_fe_s.pdf")

scatter(ds.K2O, ds.Cl,
    xerror = ds.K2O_Err,
    yerror = ds.Cl_Err,
    xlabel = "K2O",
    ylabel = "Cl",
    framestyle = :box,
    label = "",
)

scatter(ds.K2O, ds.S,
    xerror = ds.K2O_Err,
    yerror = ds.S_Err,
    xlabel = "K2O",
    ylabel = "S",
    framestyle = :box,
    label = "",
)

scatter(ds.K2O, ds.Ca,
    xerror = ds.K2O_Err,
    yerror = ds.Ca_Err,
    xlabel = "K2O",
    ylabel = "Ca",
    framestyle = :box,
    label = "",
)

# h = plot(layout=(isqrt(length(elements))+1, isqrt(length(elements))+1))
# for i in eachindex(elements)
#     elem = ds[Symbol(elements[i])]
#     if !isempty(elem) && !all(isnan, elem)
#         @info elem
#         histogram!(h[i], elem)
#     end
# end
# display(h)

## --- Plot of all XRF elements across both salars

elements = filter(x->(!contains(x, "Err") && length(x)<3), string.(keys(ds)))
elementsandoxides = ["SiO2", "Al2O3", "MgO", "K2O", elements...,]
h = plot(xticks=(1:length(elementsandoxides), elementsandoxides),
    yscale=:log10,
    ylims = (1e-4, 100),
    xrotation=-45,
    xlabel="Elements",
    ylabel="Abundance (Weight Percent)",
    title="XRF Elemental Abundance Across Atacama Salars",
    framestyle=:box,
    size = (1000, 400),
)
for i in eachindex(elementsandoxides)
    elem = ds[Symbol(elementsandoxides[i])]
    scatter!(fill(i, length(elem)), elem, alpha=0.5, label="", )
end
display(h)
savefig(h, "Fig.1.pdf")

## --- Plot of all XRF elements separated by those in Salar Elvira vs Quisquiro (color coded)

elements = filter(x->(!contains(x, "Err") && length(x)<3), string.(keys(ds)))
elementsandoxides = ["SiO2", "Al2O3", "MgO", "K2O", elements...,]
#elementsandoxides = [:SiO2, :Al2O3, :Fe, :K2O, :Rb, :Y, :Ti,  :P, :Ca, :S, :Sr, :As, :Cl, :Ba, :MgO, :Ni, :Mn, :Cu, :Zn,] # Additional, but less data :Pb, :Zr, :Nb,

h = plot(xticks=(1:length(elementsandoxides), elementsandoxides),
    yscale=:log10,
    ylims = (1e-4, 100),
    xrotation=-45,
    xlabel="Elements",
    ylabel="Abundance (Weight Percent)",
    title="XRF Elemental Abundance Across Atacama Salars",
    framestyle=:box,
    size = (1000, 400),
)

te = ds.Salar_Name .== "Elvira"
tq = ds.Salar_Name .== "Quisquiro"
for i in eachindex(elementsandoxides)
    elem = ds[Symbol(elementsandoxides[i])]
    elvira = filter(x->x>0, elem[te])
    if !isempty(elvira)
        scatter!(h, fill(i, length(elvira)), elvira, 
            alpha = 0.15,
            label = "",
            color = :red,
        )
    end
    quisquiro = filter(x->x>0, elem[tq])
    if !isempty(quisquiro)
        scatter!(h, fill(i, length(quisquiro)), quisquiro, 
            alpha = 0.15,
            label = "",
            color = :blue,
        )
    end
end

display(h)

## --- Plot of all XRF elements in Salar Elvira ONLY 
using Measures
elements = filter(x->(!contains(x, "Err") && length(x)<3), string.(keys(ds)))
#elementsandoxides = ["SiO2", "Al2O3", "MgO", "K2O", elements...,]

# Specific elements with XRF data were picked
elementsandoxides = [:SiO2, :Al2O3, :MgO, :K2O, :K, :Fe, :Mn, :Ni, :Cu, :Zn, :Rb, :Y, :Ti, :Zr, :Pb, :P, :Ca, :S, :Cl, :Sr, :As, :Ba, :V, :Cr, :Ga, :Se, :Nb, :Mo, :Rh, :Pd, :Sn,]

h = plot(xticks=(1:length(elementsandoxides), elementsandoxides),
    yscale=:log10,
    ylims = (1e-4, 100),
    xrotation=-45,
    xlabel="Elements",
    ylabel="Abundance (Weight Percent)",
    title="XRF Elemental Abundance: Salar Elvira",
    framestyle=:box,
    size = (1000, 400), 
    left_margin=15pt
)

te = ds.Salar_Name .== "Elvira"
for i in eachindex(elementsandoxides)
    elem = ds[Symbol(elementsandoxides[i])]
    elvira = filter(x->x>0, elem[te])
    if !isempty(elvira)
        scatter!(h, fill(i, length(elvira)), elvira, 
            alpha=0.5, 
            label = "",
            color = :red,
        )
    end
end

display(h)
savefig(h, "allxrf_Elvira.pdf")

## --- Plot of all XRF elements in Salar Quisquiro ONLY 

elements = filter(x->(!contains(x, "Err") && length(x)<3), string.(keys(ds)))
#elementsandoxides = ["SiO2", "Al2O3", "MgO", "K2O", elements...,]
elementsandoxides = [:SiO2, :Al2O3, :MgO, :K2O, :K, :Fe, :Mn, :Ni, :Cu, :Zn, :Rb, :Y, :Ti, :Zr, :Pb, :P, :Ca, :S, :Cl, :Sr, :As, :Ba, :V, :Cr, :Ga, :Se, :Nb, :Mo, :Rh, :Pd, :Sn,]

h = plot(xticks=(1:length(elementsandoxides), elementsandoxides),
    yscale=:log10,
    ylims = (1e-4, 100),
    xrotation=-45,
    xlabel="Elements",
    ylabel="Abundance (Weight Percent)",
    title="XRF Elemental Abundance: Salar Quisquiro",
    framestyle=:box,
    size = (1000, 400),
    left_margin=15pt
)

te = ds.Salar_Name .== "Quisquiro"
for i in eachindex(elementsandoxides)
    elem = ds[Symbol(elementsandoxides[i])]
    quisquiro = filter(x->x>0, elem[te])
    if !isempty(quisquiro)
        scatter!(h, fill(i, length(quisquiro)), quisquiro, 
            alpha=0.5, 
            label = "",
            color = :blue,
        )
    end
end
display(h)
savefig(h, "allxrf_Quisquiro.pdf")

## --- Plot histograms comparing both salars, absolute counts
elements = filter(x->(!contains(x, "Err") && length(x)<3), string.(keys(ds)))
elementsandoxides = ["SiO2", "Al2O3", "MgO", "K2O", elements...,]
logbins = -4:0.2:2
for i in eachindex(elementsandoxides)
    te = ds.Salar_Name .== "Elvira"
    tq = ds.Salar_Name .== "Quisquiro"
    elem = ds[Symbol(elementsandoxides[i])]
    elvira = filter(x->x>0, elem[te])
    if !isempty(elvira)
        h = histogram(log10.(elvira), bins=logbins,
            xlabel="log10 $(elementsandoxides[i])",
            ylabel="N",
            xlims = (-4, 2),
            framestyle=:box,
            alpha = 0.5,
            label = "Elvira",
            color = :darkred,
        )
        display(h)
        savefig(h, "Histogram $(elementsandoxides[i]).pdf")
    end
    quisquiro = filter(x->x>0, elem[tq])
    if !isempty(quisquiro)
        histogram!(log10.(quisquiro), bins=logbins,
            xlabel="log10 $(elementsandoxides[i])",
            ylabel="N",
            xlims = (-4, 2),
            framestyle=:box,
            alpha = 0.5,
            label = "Quisquiro",
            color = :darkblue,
        )
        display(h)
        savefig(h, "Histogram $(elementsandoxides[i]).pdf")
    end
end

## --- Plot histograms comparing both salars, normalized data
elements = filter(x->(!contains(x, "Err") && length(x)<3), string.(keys(ds)))
elementsandoxides = ["SiO2", "Al2O3", "MgO", "K2O", elements...,]
logbins = -4:0.2:2
for i in eachindex(elementsandoxides)
    te = ds.Salar_Name .== "Elvira"
    tq = ds.Salar_Name .== "Quisquiro"
    elem = ds[Symbol(elementsandoxides[i])]
    elvira = filter(x->x>0, elem[te])
    if !isempty(elvira)
        h = histogram(log10.(elvira), bins=logbins, normalized=true,
            xlabel="log10 $(elementsandoxides[i])",
            ylabel="Probability density",
            xlims = (-4, 2),
            framestyle=:box,
            alpha = 0.5,
            label = "Elvira",
            color = :red,
        )
        display(h)
       savefig(h, "Histogram $(elementsandoxides[i]).pdf")
    end
    quisquiro = filter(x->x>0, elem[tq])
    if !isempty(quisquiro)
        histogram!(log10.(quisquiro), bins=logbins, normalized=true,
            xlabel="log10 $(elementsandoxides[i])",
            ylabel="Probability density",
            xlims = (-4, 2),
            framestyle=:box,
            alpha = 0.5,
            label = "Quisquiro",
            color = :blue,
        )
        display(h)
        savefig(h, "Histogram $(elementsandoxides[i]) normalized.pdf")
    end
end

## --- Plots of XRF surface vs subsurface measurements, normalized
elements = filter(x->(!contains(x, "Err") && length(x)<3), string.(keys(ds)))
elementsandoxides = ["SiO2", "Al2O3", "MgO", "K2O", elements...,]
logbins = -4:0.2:2
for i in eachindex(elementsandoxides)
    te = ds.Location .== "Surface"
    tq = ds.Location .== "Subsurface"
    elem = ds[Symbol(elementsandoxides[i])]
    surface = filter(x->x>0, elem[te])
    h = plot(
        xlabel="log10 $(elementsandoxides[i])",
        ylabel="Probability density",
        xlims = (-4, 2),
        framestyle=:box,
    )
    if !isempty(surface)
        histogram!(h, log10.(surface), bins=logbins, normalized=true,
            alpha = 0.5,
            label = "Surface",
            color = :cyan,
        )

    end
    subsurface = filter(x->x>0, elem[tq])
    if !isempty(subsurface)
        histogram!(log10.(subsurface), bins=logbins, normalized=true,
            alpha = 0.5,
            label = "Subsurface",
            color = :purple,
        )
    end
    display(h)
    savefig(h, "Histogram $(elementsandoxides[i]) Location.pdf")
end

## --- PCA statistical test on all XRF data
using MultivariateStats

datamatrix = unelementify(ds, elementsandoxides, floatout=true)
datamatrix[datamatrix .== 0] .= NaN
nanstandardize!(datamatrix, dims=1) # Standardize colums to zero mean and unit standard deviation
datamatrix[isnan.(datamatrix)] .= 0

M = fit(PCA, datamatrix'; pratio=1, maxoutdim=2)
proj = projection(M)
data_transformed = transform(M, datamatrix')

contribution = vec(sum(abs.(proj), dims=2))
I = sortperm(contribution, rev=true)

h = plot(xlabel="PC1", ylabel="PC2", framestyle=:box) # A few formatting options
te = ds.Salar_Name .== "Elvira"
plot!(h, data_transformed[1,te], data_transformed[2,te], seriestype=:scatter, label="Elvira",color=:red)
tq = ds.Salar_Name .== "Quisquiro"
plot!(h, data_transformed[1,tq], data_transformed[2,tq], seriestype=:scatter, label="Quisquiro",color=:blue)
for i in I[1:6] # Four most important elements
    plot!([0,proj[i,1]], [0,proj[i,2]], arrow=true, label=elementsandoxides[i], legend=:topleft)
end
display(h)
savefig(h, "PCA.pdf")

## --- PCA test as above, but separated by surface and subsurface measurements

h = plot(xlabel="PC1", ylabel="PC2", framestyle=:box) # A few formatting options
te =ds.Location .== "Surface"
plot!(h, data_transformed[1,te], data_transformed[2,te], seriestype=:scatter, label="Surface",color=:cyan)
tq = ds.Location .== "Subsurface"
plot!(h, data_transformed[1,tq], data_transformed[2,tq], seriestype=:scatter, label="Subsurface",color=:purple)
for i in I[1:6] # Four most important elements
    plot!([0,proj[i,1]], [0,proj[i,2]], arrow=true, label=elementsandoxides[i], legend=:topleft)
end
display(h)
savefig(h, "PCA_surface_andsub.pdf")

## --- Plot loadings of PC1 and PC2 from first PCA test

contribution = vec(sum(abs2.(proj), dims=2))
tl = contribution .> 0.0001
h = plot(1:count(tl), proj[tl,1],
    xticks=(1:count(tl), elementsandoxides[tl]),
    xrotation=-45,
    xlabel="Elements",
    ylabel="Loading",
    title="PCA loadings",
    framestyle=:box,
    size = (800, 400),
    label = "PC1 loadings",
    left_margin=10pt,
    bottom_margin=15pt,
    linewidth=2.0,
    linecolor=:orange
)
plot!(1:count(tl), proj[tl,2],
    label = "PC2 loadings",
    linewidth=2.0,
    linecolor=:green
)
display(h)
savefig(h, "PCA loadings.pdf")

## --- find outlier in first PCA plot
t = findall((0.6 .< data_transformed[1,:] .< 1.5) .& (2 .< data_transformed[2,:] .< 2.7))
# Ans is 23, row 23 is gypsum crystal

## --- Pearson correlation matrix

# elementswithdata = filter(e->count(x->x>0, ds[e]) > 5, Symbol.(elementsandoxides))
elementswithdata = [:SiO2, :Al2O3, :Fe, :K2O, :Rb, :Y, :Ti,  :P, :Ca, :S, :Sr, :As, :Cl, :Ba, :MgO, :Ni, :Mn, :Cu, :Zn,] # Additional, but less data :Pb, :Zr, :Nb,

datamatrix = unelementify(ds, elementswithdata, floatout=true)
datamatrix[datamatrix .== 0] .= NaN
C = nancor(datamatrix, dims=1)
h = heatmap(C, 
    yflip=true,
    framestyle=:box,
    xticks = (1:length(elementswithdata), elementswithdata),
    yticks = (1:length(elementswithdata), elementswithdata),
    xrotation = 45,
    colorbartitle = "Correlation",
)
savefig(h, "correlation_matrix.pdf")