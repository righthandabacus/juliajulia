using ArgParse

using UnicodePlots # for diagnostic

# Allow argparse to understand complex number
function ArgParse.parse_item(::Type{Complex}, x::AbstractString)
    return Complex(eval(parse(x)))
end

# Entry point
function main()
    parsed_args = parse_commandline()
    max_z = parsed_args["max"]
    min_z = parsed_args["min"]
    pix_w = parsed_args["dim"][1]
    pix_h = parsed_args["dim"][2]
    mandelbrot = parsed_args["mandelbrot"]
    c0 = parsed_args["C"]
    filename = parsed_args["filename"]
    println("Plot $pix_w x $pix_h from $min_z to $max_z with c0 = $c0 to $filename (Mandelbrot=$mandelbrot)")
    if mandelbrot
        plot_mandelbrot(pix_w, pix_h, min_z, max_z, c0, filename)
    else
        plot_julia(pix_w, pix_h, min_z, max_z, c0, filename)
    end
end

# Parse command line arguments
function parse_commandline()
    s = ArgParseSettings(
            description = "Plotting Julia set"
            ,version = "0.1"
        )
    @add_arg_table s begin
         "--dim", "-D"
            help = "Width and height of the output in pixels"
            nargs = 2
            arg_type = Int
            default = [1024, 1024]
            range_tester = (x -> x>0)
         "--mandelbrot"
            help = "Plot Mandelbrot set instead of Julia set"
            action = :store_true
         "-C"
            help = "If Julia, the constant c in formula z^2+c; if Mandelbrot, the value of z_0 used. In form of real and imaginary parts"
            arg_type = Complex
            default = -1.2+0.156im # nice values: -0.8+0.16im, -1.2+0.156im, -0.7+0.27015im
         "--min", "-m"
            help = "Lower left complex coordinate of output"
            arg_type = Complex
            default = -1.5-1.5im
         "--max", "-M"
            help = "Upper right complex coordinate of output"
            arg_type = Complex
            default = 1.5+1.5im
         "filename"
            help = "Output filename"
            required = true
    end
    return parse_args(ARGS, s)
end

# Given z,c return the iteration count just before going beyond the escape
# radius er, on formula z := z^2 + c, with iteration count bounded by nmax.
# Points inside the basin will report Inf
function count_escape(z, c, er, nmax::Int64)
    for n = 1:nmax
        if abs(z) > er
            return n-1
        end
        z = z^2 + c
    end
    return Inf
end

# Given z,c return the distance just before going beyond the escape radius er
# using the Hubbard-Donady potential function, on Julia set with formula
# z := z^2 + c, with iteration count bounded by nmax
function juliadistance(z, c, er, nmax::Int64)
    dz = 1 + 0im
    for n = 1:nmax
        dz = 2z*dz
        z = z^2 + c
        if abs(z) > er
            return abs(z)*log(abs(z))/abs(dz)
        end
    end
    return 0
end

# Given z,c return the distance just before going beyond the escape radius er
# using the Hubbard-Donady potential function, on Mandelbrot set with formula
# z := z^2 + c, with iteration count bounded by nmax
function mandelbrotdistinace(z, c, er, nmax::Int64)
    dz = 0
    for n = 1:nmax
        dz = 2z*dz + 1
        z = z^2 + c
        if abs(z) > er
            return abs(z)*log(abs(z))/abs(dz)
        end
    end
    return 0
end

# Convert HSV color into RGB
# @param h is in degree, [0,360], out of range will be mod 360
# @param s,v are in intensity [0,1]
# @return 3-tuple of RGB, each in intensity [0,1]
function hsv2rgb(h, s, v)
    h = mod(h, 360) / 60 # mod to convert h into [0,360], then find the 60-degree sector
    h_int = floor(h)     # color sector, 60-degree each
    h_frac = h - h_int   # offset in sector

    p = v * (1 - s)
    q = v * (1 - s * h_frac)
    t = v * (1 - s * (1-h_frac))

    if h_int == 0 (return [v, t, p])
    elseif h_int == 1 (return [q, v, p])
    elseif h_int == 2 (return [p, v, t])
    elseif h_int == 3 (return [p, q, v])
    elseif h_int == 4 (return [t, p, v])
    elseif h_int == 5 (return [v, p, q])
    else (return [0, 0, 0])
    end
end

# Convert complex number to HSV colour
# source: https://en.wikibooks.org/wiki/Color_Theory/Color_gradient
function z2hsv(z)
    # express z as z == zmod * exp(zarg)
    zarg = mod(angle(z), 2pi) # in radians, from 0 to +2pi
    zmod = abs(z)

    # hue, depend on arg(z), express in degree
    hue = zarg / 2pi * 360

    # let zmod in [e^n, e^{n+1}), then find k as intercept fraction
    range_end = exp(ceil(log(zmod)))
    range_begin = rangee/e
    k = (zmod-range_begin)/(range_end-range_begin)

    # saturation, depends on zmod, within [0,1]
    sat = k<0.5 ? k*2 : 1-(k-0.5)*2
    sat = 1 - (1-sat)^3
    sat = 0.4 + sat * 0.6

    # value, depends on zmod, within [0,1]
    val = k>0.5 ? k*2 : 1-(k-0.5)*2
    val = 1 - val
    val = 1 - (1-val)^3
    val = 0.6 + val*0.4

    return [hue, sat, val]
end

# Print histogram to screen to help investigate
function print_histogram(img, fn=x -> x)
    print(histogram(map(fn, vec(img)), bins=25))
end

# Save image into PNG file
# @param img is 2D array of pixels, each pixel is a [R,G,B] array with each color in [0,1] intensity
# @param outf is the output filename, suffix .png will be added if necessary
function save_pngfile(img, outf::AbstractString)
    # Write as PGM file in binary, color mode
    pgmname = tempname() * ".pgm"
    s = open(pgmname, "w")
    write(s, "P6\n")
    n,m = size(img)
    write(s, "$m $n 255\n")
    Pct2Byte(x) = UInt16(floor(x*255)) % UInt8
    for i=1:n, j=1:m
        rgb = map(Pct2Byte, img[i,j])
        write(s, rgb[1])
        write(s, rgb[2])
        write(s, rgb[3])
    end
    close(s)
    # Use ImageMagick to convert the PGM file into PNG
    outf = strip(outf)
    if !endswith(uppercase(outf), ".PNG")
        outf = outf * ".png"
    end
    run(`convert $pgmname $outf`)
    println("Wrote to $outf")
    rm(pgmname)
end

# Scale array M to integers between 0 and 359. Inf values replaced by 0.
function scale360(M)
    maxM = maximum(filter(x -> x != Inf, M)) + eps()
    for i in eachindex(M)
        M[i] = M[i] == Inf ? 0 : floor(360*M[i]/maxM)
    end
    return M
end

# Transform boolean into black (false) and white (true) pixels
function black_and_white(M)
    for i in eachindex(M)
        M[i] = M[i]?[1,1,1]:[0,0,0]
    end
    return M
end

# Generic function to plot
# @param pix_w, pix_h are the pixel size of the image
# @param min_z, max_z are the coordinates of the two corners, in complex number
# @param c0 is the parameter of the formula
# @param filename is the filename of PNG output
# @param compute(z,c0) is a function to compute for coorindate z, with parameter c0
# @param transform(M) is a function to transform each pixel of bitmap M into RGB format, with values in [0,1]
function generic_plot(pix_w, pix_h, min_z, max_z, c0, filename, compute, transform)
    min_x = real(min_z)
    min_y = imag(min_z)
    max_x = real(max_z)
    max_y = imag(max_z)
    M = Array(Any, pix_h, pix_w)
    t_0 = time()
    for y=1:pix_h, x=1:pix_w
        z = complex(x*(max_x-min_x)/pix_w+min_x, max_y - y*(max_y-min_y)/pix_h)
        M[y,x] = compute(z, c0)
    end
    t_1 = time()
    save_pngfile(transform(M), filename)
    println("Finished in $(t_1-t_0) seconds.")
end

# Plot black-and-white image of the fill-in Mandelbrot set
function plot_mandelbrot_bw(pix_w, pix_h, min_z, max_z, c0, filename)
    compute(z,c) = count_escape(c,z,4,500) != Inf
    transform = black_and_white
    generic_plot(pix_w, pix_h, min_z, max_z, c0, filename, compute, transform)
end

# Plot black-and-white image of the fill-in Julia set
function plot_julia_bw(pix_w, pix_h, min_z, max_z, c0, filename)
    compute(z,c) = count_escape(z,c,4,500) != Inf
    transform = black_and_white
    generic_plot(pix_w, pix_h, min_z, max_z, c0, filename, compute, transform)
end

function _old()
    # Julia set parameters: nice values -0.8+0.16im, -1.2+0.156im, -0.7+0.27015im
    c0 = -1.2+0.156im

    # Image parameters
    min_x = -0.3
    max_x = 0.3
    min_y = 0.0
    max_y = 0.6
    h = 1024; w = 1024
    pgm_name = "juliaset.pgm"

    M = Array(Float64, h, w)
    N = Array(Float64, h, w)

    t_0 = time()
    for y=1:h, x=1:w
        # scale (x,y) to [-1.5,+1.5]
        z = complex(x*(max_x-min_x)/w+min_x, max_y - y*(max_y-min_y)/h)
        M[y,x] = mandelbrotjuliacount(z, c0, 4.0, 360)
        N[y,x] = juliadistance(z, c0, 2.0, 360)
    end
    t = time()
    #Mx = filter(x -> x != Inf, M)
    #print("min dist = $(minimum(Mx)), max dist = $(maximum(Mx)), mean = $(sum(Mx)/length(Mx)))\n")
    #scale360(M)
    #scale360(N)
    save_pgmfile(scale360(M+N), pgm_name)
    print("Written $pgm_name\nFinished in $(t-t_0) seconds.\n")
end

function plot_mandelbrot(pix_w, pix_h, min_z, max_z, c0, filename)
    plot_mandelbrot_bw(pix_w, pix_h, min_z, max_z, c0, filename)
end

function plot_julia(pix_w, pix_h, min_z, max_z, c0, filename)
    plot_julia_bw(pix_w, pix_h, min_z, max_z, c0, filename)
end

# TODO normalised iteration count
# https://en.wikipedia.org/wiki/Mandelbrot_set#Continuous_.28smooth.29_coloring
# also:
# http://lodev.org/cgtutor/juliamandelbrot.html
main()

# vim:set et ts=4 sw=4:
