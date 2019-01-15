using JLD, Plots, LaTeXStrings
font = Plots.font("DejaVu Sans", 24)
pyplot(size=(1200,800),linewidth=4,xscale=:log10,yscale=:log10,
       xtickfont=font,ytickfont=font,markersize=8,
       legendfont=font,guidefont=font,titlefont=font,
       xlabel="Δx",ylabel="||error||",
       title="MPB Equation Convergence",grid=false)

function main()
    # import jld file
    jldfile = "data/TEMP_mpb2D_validation.jld"

    VL1arr   = load(jldfile,"VL1arr")
    VL2arr   = load(jldfile,"VL2arr")
    VLInfarr = load(jldfile,"VLInfarr")
    PL1arr   = load(jldfile,"PL1arr")
    PL2arr   = load(jldfile,"PL2arr")
    PLInfarr = load(jldfile,"PLInfarr")
    harr     = load(jldfile,"harr")
    κarr     = load(jldfile,"κarr")
    timearr  = load(jldfile,"timearr")

    # create plots
    # error plots
    p1 = plot(harr,PL1arr,color=:red,linestyle=:dash,label="L1 (pressure)")
    plot!(p1,harr,PL2arr,color=:blue,linestyle=:dash,label="L2")
    plot!(p1,harr,PLInfarr,color=:green,linestyle=:dash,label="Linf")
    scatter!(p1,harr,PL1arr,color=:red,label="")
    scatter!(p1,harr,PL2arr,color=:blue,label="")
    scatter!(p1,harr,PLInfarr,color=:green,label="")
    plot!(p1,harr,harr.^2,color=:black,linestyle=:dash,label="O(Δx^2)")

    plot!(p1,harr,VL1arr,color=:red,label="L1 (velocity)")
    plot!(p1,harr,VL2arr,color=:blue,label="L2")
    plot!(p1,harr,VLInfarr,color=:green,label="Linf")
    scatter!(p1,harr,VL1arr,color=:red,label="")
    scatter!(p1,harr,VL2arr,color=:blue,label="")
    scatter!(p1,harr,VLInfarr,color=:green,label="") 
    plot!(p1,harr,0.01*harr.^4,color=:black,label="O(Δx^4)")   

    # time plots
    p2 = plot(title="MPB Time-Accuracy",xlabel="||error||",ylabel="elapsed time (s)")
    colorArr = [:red,:green,:blue,:orange,:yellow,:purple]
    for i=1:length(timearr)
        scatter!([PLInfarr[i]],[timearr[i]],
                 color=colorArr[i],label="",markershape=:star,
                 markersize=16)
        scatter!([VLInfarr[i]],[timearr[i]],
                 color=colorArr[i],label="",markershape=:circle,
                 markersize=12)
    end

    # condition number plot
    p3 = plot(harr,κarr,label="",linewidth=8,color=:blue)
        scatter!(p3,harr,κarr,label="",ylabel=L"\kappa",
            title="MPB condition number",markersize=10,
            color=:blue)
        plot!(p3,harr,10^6*harr.^(-3.0),label="Δx^-3",linewidth=2,color=:black)

    # export plots
    run(`mkdir -p figures`)
    plot(p1)
    savefig(p1,"figures/TEMPMPB2D_validation.png")
    plot(p2)
    savefig(p2,"figures/TEMPMPB2D_time_accuracy.png")
    plot(p3)
    savefig(p3,"figures/TEMPMPB2D_cond.png")

    nothing
end

main()