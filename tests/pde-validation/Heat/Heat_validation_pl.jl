using JLD, Plots, LaTeXStrings
font = Plots.font("DejaVu Sans", 24)
pyplot(size=(1200,800),linewidth=3,xscale=:log10,yscale=:log10,
       xtickfont=font,ytickfont=font,markersize=8,
       legendfont=font,guidefont=font,titlefont=font,
       ylim = (5e-6,1e-3),grid=false)

function main()
  jldfile = "data/TEMPHeat_validation.jld"
  L1arr = load(jldfile, "L1arr")
  L2arr = load(jldfile, "L2arr")
  LInfarr = load(jldfile, "LInfarr")
  timearr = load(jldfile, "timearr")
  dtarr = load(jldfile, "dtarr")

  # create plots
  p1 = Plots.plot(dtarr,L1arr,label="L1",
            title="Heat Equation Convergence N=16")
  Plots.plot!(p1,dtarr,L2arr,label="L2")
  Plots.plot!(p1,dtarr,LInfarr,label="Linf")
  scatter!(p1,dtarr,L1arr,color=:black,label="")
  scatter!(p1,dtarr,L2arr,color=:black,label="")
  scatter!(p1,dtarr,LInfarr,color=:black,label="")
  Plots.plot!(p1,dtarr,dtarr,color=:black,label="O(Δt)")
  Plots.plot(p1)

  # save plots
  run(`mkdir -p figures`)
  savefig("figures/HEATvalidation.png")
end

main()
