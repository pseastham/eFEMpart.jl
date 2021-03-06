using JLD, Plots, LaTeXStrings
font = Plots.font("DejaVu Sans", 24)
pyplot(size=(1200,800),linewidth=4,xscale=:log10,yscale=:log10,
       xtickfont=font,ytickfont=font,markersize=8,
       legendfont=font,guidefont=font,titlefont=font,
       xlabel="Δx",ylabel="||error||",
       title="Poisson Equation Convergence",grid=false)

function main()
    # import jld files
    jldfile = "data/TEMPpoisson_validation.jld"

    o1L1arr   = load(jldfile,"o1L1arr")
    o1L2arr   = load(jldfile,"o1L2arr")
    o1LInfarr = load(jldfile,"o1LInfarr")
    o2L1arr   = load(jldfile,"o2L1arr")
    o2L2arr   = load(jldfile,"o2L2arr")
    o2LInfarr = load(jldfile,"o2LInfarr")
    harr      = load(jldfile,"harr")
    κarr      = load(jldfile,"κarr")
    o1timearr = load(jldfile,"o1timearr")
    o2timearr = load(jldfile,"o2timearr")

    # plot files
    MS = 10         # marker size
    LW = 3          # line width
    FS = 15         # font size
    FStitle = FS+5

    # create plots
    # error plots
    p1 = plot(harr,o1L1arr,color=:red,linestyle=:dash,label="L1 (linear ϕ)")
    plot!(p1,harr,o1L2arr,color=:blue,linestyle=:dash,label="L2")
    plot!(p1,harr,o1LInfarr,color=:green,linestyle=:dash,label="Linf")
    scatter!(p1,harr,o1L1arr,label="",color=:black)
    scatter!(p1,harr,o1L2arr,label="",color=:black)
    scatter!(p1,harr,o1LInfarr,label="",color=:black)
    plot!(p1,harr,harr.^2,color=:black,linestyle=:dash,label="O(Δx^2)")

    plot!(p1,harr,o2L1arr,color=:red,label="L1 (quadratic ϕ)")
    plot!(p1,harr,o2L2arr,color=:blue,label="L2")
    plot!(p1,harr,o2LInfarr,color=:green,label="Linf")
    scatter!(p1,harr,o2L1arr,color=:black,label="")
    scatter!(p1,harr,o2L2arr,color=:black,label="")
    scatter!(p1,harr,o2LInfarr,color=:black,label="") 
    plot!(p1,harr,0.01*harr.^4,color=:black,label="O(Δx^4)")   

    # time plots
    p2 = plot(title="Poisson Time-Accuracy",xlabel="||error||",ylabel="elapsed time (s)")
    colorArr = [:red,:green,:blue,:orange,:yellow,:purple]
    for i=1:length(o1timearr)
        scatter!([o1LInfarr[i]],[o1timearr[i]],
                 color=colorArr[i],label="",markershape=:star,
                 markersize=16)
        scatter!([o2LInfarr[i]],[o2timearr[i]],
                 color=colorArr[i],label="",markershape=:circle,
                 markersize=12)
    end

    # condition number plot
    p3 = plot(harr,κarr,label="",linewidth=LW,color=:blue)
    scatter!(p3,harr,κarr,label="",ylabel=L"\kappa",
            title="Poisson condition number",markersize=MS,
            color=:blue)

    # export plots
    run(`mkdir -p figures`)
    plot(p1)
    savefig(p1,"figures/TEMPpoisson_validation.png")
    plot(p2)
    savefig(p2,"figures/TEMPpoisson_time_accuracy.png")
    plot(p3)
    savefig(p3,"figures/TEMP_poisson_cond.png")

    nothing
end

main()
