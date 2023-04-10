if get(ENV, "HF_BUILD_FROM_SOURCE", "false") == "true"
    extension = Sys.isapple() ? "dylib" : Sys.islinux() ? "so" : Sys.iswindows() ? "dll" : ""
    make = Sys.iswindows() ? "mingw32-make" : "make"
    script = """
        set -e
        set -x
        if [ -d "HypergeometricFunctions" ]; then
            cd HypergeometricFunctions
            #git fetch
            #git checkout master
            #git pull
            cd ..
        #else
            #git clone https://github.com/MikaelSlevinsky/FastTransforms.git FastTransforms
        fi
        cd HypergeometricFunctions
        export CC=gcc-8
        export CFLAGS="-O3"
        $make lib
        cd ..
        mv -f HypergeometricFunctions/libhypergeometricfunctions.$extension libhypergeometricfunctions.$extension
    """
    try
        run(`bash -c $(script)`)
    catch
        error(
            "HypergeometricFunctions could not be properly installed.\n Please check that you have all dependencies installed. " *
            "Sample installation of dependencies:\n" *
            (Sys.isapple() ? "On MacOS\n\tbrew install mpfr\n" :
             Sys.islinux() ? "On Linux\n\tsudo apt-get install libmpfr-dev\n" :
             Sys.iswindows() ? "On Windows\n\tvcpkg install mpir:x64-windows mpfr:x64-windows\n" :
             "On your platform, please consider opening a pull request to add support to build from source.\n")
        )
    end
    println("HypergeometricFunctions built from source.")
else
    println("HypergeometricFunctions using precompiled binaries.")
end
