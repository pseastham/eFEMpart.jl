# install script for Linux & MacOS machines.

using Pkg

function install()
    checkIneFEMpart()
    intro()
    skipline()
    checkOK()
    skipline()
    checkJuliaVersion()
    skipline()
    addStartupFile()
    skipline()
    install_fgt()
    skipline()
    addModules()
    skipline()
    finished()
end

function checkIneFEMpart()
    if !(basename(pwd()) == "eFEMpart")
        error("current directory is not eFEMpart")
    end
end

function intro()
    println("You have chosen to install eFEMpart,")
    println("a package for implementing finite element")
    println("code with particle simulators in Julia")
end

"""
    input(prompt::String="")::String

Read a string from STDIN. The trailing newline is stripped.

The prompt string, if given, is printed to standard output without a
trailing newline before reading input.
"""
function input(prompt::String)
    print(prompt)
    return chomp(readline())
end

function checkOK()
    println("This install script does the following:")
    println("   - Adds eFEMpart packageto your Julia LOAD_PATH")
    println("   - Downloads and installs figtree package")
    println("   - Adds the following dependencies:")
    println("      + JLD")
    println("      + Plots")
    println("      + PyPlot")
    println("      + PyCall")
    println("      + LaTeXStrings")
    println("      + Parameters")

    keepasking = true
    while keepasking
        OK = input("Is this OK? (y/n) ")
        
        if OK == "y"
            keepasking = false
        elseif OK == "n"
            error("User did not want to install")
        else
            println("Sorry, that isn't an option.")
        end
    end

    println()
    println("Great! Installation will begin now.")
    println("This could take several minutes...")
end

function skipline()
    println()
end

function checkJuliaVersion()
    print("Checking correct Julia version...")
    if VERSION>v"1.0.0"
        println("OK")
    else
        println("FAILED")
        error("Julia version ",VERSION," needs to be at least 1.0")
    end
end

function addStartupFile()
    print("Adding eFEMpart to your Julia PATH...")
    eFEMpartPATH = string(pwd(),"/src/")
    newLoadPath = string("push!(LOAD_PATH,\"",eFEMpartPATH,"\")")

    # create startup.jl file
    mkpath("$(homedir())/.julia/config")
    touch("$(homedir())/.julia/config/startup.jl")

    fr = open("$(homedir())/.julia/config/startup.jl","r")

    # check that LOAD_PATH doesn't already contain what you're about to add
    alreadyThere = false
    lines = readlines(fr)
    for i=1:length(lines)
        if occursin(newLoadPath, lines[i])
            alreadyThere = true
        end
    end
    close(fr)
    
    # add startup file
    if !(alreadyThere)
        fa = open("$(homedir())/.julia/config/startup.jl","a")
        write(fa,"\n",newLoadPath)
        close(fa)
    end

    # load startup file -- so that you don't have to restart julia
    include("$(homedir())/.julia/config/startup.jl")

    println("OK")
end

function addModules()
    # JLD
    println("Adding JLD package...")
    Pkg.add("JLD");
   
    # Plots
    println("Adding Plots package...")
    Pkg.add("Plots");
   
    # PyPlot
    println("Adding PyPlot package...")
    Pkg.add("PyPlot");
   
    # PyCall
    println("Adding PyCall package...")
    Pkg.add("PyCall");
   
    # LaTeXStrings
    println("Adding LaTeXStrings package...")
    Pkg.add("LaTeXStrings");

    # Parameters
    println("Adding Parameters package...")
    Pkg.add("Parameters");
end

function install_fgt()
    println("installing figtree...requires system owner permission")
    run(`sudo mkdir -p src/part`)
    cd("src/part")
    # download zip
    downloadurl = "https://sourceforge.net/projects/figtree/files/latest/download"
    run(`curl -OL $(downloadurl)`);
 
    # unzip
    run(`unzip download`);
    #run(`unzip download && mv figtree download`)
 
    # compile
    run(`make -C figtree-0.9.3`);
 
    # move relevant *.so file
    run(`mv figtree-0.9.3/lib/libfigtree.so libfigtree.so`)
    run(`mv figtree-0.9.3/lib/libann_figtree_version.so libann_figtree_version.so`)

    # remove unnecessary package
    run(`rm -r figtree-0.9.3/ download`)
 
    # add *.so folder to library path
    #pwdir = pwd()
    #bsh = "$(homedir())/.bash_profile"
    #f=open(bsh,"a")
    #println(f,"export LD_LIBRARY_PATH=$(pwdir)") 
    #close(f)
 
    # create file so systme recognizes *.so 
    run(`sudo touch /etc/ld.so.conf.d/myLibs.conf`)
    # give permission to edit
    run(`sudo chmod 777 /etc/ld.so.conf.d/myLibs.conf`)
    f=open("/etc/ld.so.conf.d/myLibs.conf","a")
    println(f,pwd())
    close(f)
    # close permission to edit
    run(`sudo chmod 400 /etc/ld.so.conf.d/myLibs.conf`)

    # reload shared object paths (or something like that ?)
    run(`sudo ldconfig`)

    cd("../..")

    println("succesfully installed figtree!")
end
function finished()
    skipline()
    println("All components installed correctly!")
    println("You can now use eFEMpart!")
end


install()

