include("mydefs.jl")

using WAV,Winston,Distributions

sample_rate = 44100 ÷ 1

function play(fun::Function,dur::Float64,env::Function=(x->1.0),sr::Int64=sample_rate)
    list = [ env(t)*fun(t) for t=0:1/sr:dur ]
    wavplay(list,sr)
    return(list)
end

a3 = 220;
a4 = 440;
c4 = 220 * 2^(1/4);

sin1(t::Float64) = sin(2 * pi * t)

# function fm(sig::Function,mod::Function,width::Number)::Function
#     arg = function(t::Float64)
#         int,_ = quadgk(mod,0.0,t;)
#         t + width * int
#     end
#     function(t::Float64)
#         sig(arg(t))
#     end
# end

function fm_lin(sig::Function,
                mod::Array{Array{Float64,1},1},
                width::Number)::Function
    intfun = interpolation_lin_int(mod)
    arg = function(t::Float64)
        #println(t)
        t + width * intfun(t)
    end
    function(t::Float64)
        sig(arg(t))
    end
end

function sineaddm(sig::Function,mod::Function,width::Number)::Function
    function(t::Float64)
        sig(t)+width*mod(t)
    end
end
function sineam(sig::Function,mod::Function,width::Number)::Function
    function(t::Float64)
        sig(t)*width*mod(t)
    end
end

function decibelsToAmp(d::Float64,min::Number=-1000)::Float64
    d>min ? 10^(d/20) : 0.0
end
function ampToDecibel(a::Float64)::Float64
    20 * log(10,abs(a))
end

function make_env(lst,dur,start=0.0)
    function(t::Float64)
        if t >= start + dur
            return(0.0)
        else
            r = (t-start)/dur
            filt = filter(x -> x[1] <= r,lst)
            ind = length(filt)
            a,b = lst[ind],lst[ind+1]
            (r-a[1])/(b[1]-a[1])*(b[2]-a[2])+a[2]
        end
    end
end

function make_note_env(attack,decay,sustain,kill)
    function(dur::Float64,start::Float64=0.0)
        make_env([ [0 0],
                   [attack[1]/dur attack[2]],
                   [(attack[1]+decay)/dur sustain],
                   [1-kill/dur sustain],
                   [1 0] ],
                 dur, start )
    end
end

recorder_env = make_note_env([0.2 1],0,1,0.2)

function voice(amp,v::Function)
    fmfun = fm(sin1,v,1)
    function(t::Float64)
        amp(t) * fmfun(t)
    end
end
function voice(amp,v::Array{Array{Float64,1},1})
    fmfun = fm_lin(sin1,v,1)
    function(t::Float64)
        amp(t) * fmfun(t)
    end
end

function lcdiff(x1::Float64,x2::Float64)::Float64
    lst = [x1-x2,2*x1-x2,x1-2*x2]
    sort!(lst,by=abs)
    -lst[1]
end

function force_one(con::Float64,x1::Float64,x2::Float64)::Float64
    g = force(con,lcdiff(x1,x2),Inf)
    x2>x2 ? -g : g
end

function force(con::Float64,g::Float64,lim::Float64=40)::Float64
    g<lim ? con * (400 - (g-20)^2) : 0
end

function force(con::Float64,lst::Array{Float64,1})::Array{Float64,1}
    num = length(lst)
    res = [ 0.0 for _=1:num ]
    for i=1:num-1
        for j=i+1:num
            f = -force_one(con,lst[i],lst[j])
            res[i] += f
            res[j] -= f
        end
    end
    res
end

function random(num::Int64,time::Float64,fcon::Float64,range=100.0:0.1:1760.0)
    amps = reverse(sort(rand(0.5:0.1:5,num)))
    hts = sort(rand(range,num))
    oldhts = copy(hts)
    evolve(hts,amps,time,fcon,range)
end

function evolve(start::Array{Float64,1},
                startamps::Array{Float64,1},
                time::Float64,
                fcon::Float64,
                rcon::Float64,
                tcon::Float64)
    num = length(start)
    hts,amps = sort_start(copy(start),copy(startamps))
    hts,amps = [hts...],[amps...]
    #oldhts = copy(hts) .* exp.(rcon * rand(-1.:2:1,num)) ./ tcon
    vels = exp.(rcon * rand(-1.:2:1,num))
    println("start vels:",vels)
    tab,amptab = Array{Float64,1}[],Array{Float64,1}[]
    fvec = rand(Normal(fcon,0.1*fcon),num)
    t = 0
    inc = tcon
    while t <= time + inc
        htst = vcat([t],hts)
        #println(htst)
        push!(tab,htst)
        #println(length(amptab))
        ampt = vcat([t],amps)
        #println(ampt)
        push!(amptab,ampt)
        #println(length(amptab))
        # fs = force(fcon/(num-1),hts) ./ amps
        # #gs = (hts ./ start) .^ -gcon
        # gs = map(x -> 1 ./ log(2,x), hts ./ 20) .+
        #      map(x -> 1 ./ log(2,x), hts ./ 20000)
        # vels = (hts ./ oldhts) .* (2 .^ fs) .* (2 .^ (gcon * gs))
        fs = fvec .* spring_force(start,hts)
        #println("fs = ",fs)
        vels .*= exp.(fs)
        #println("vels = ",vels)
        oldhts = copy(hts)
        vfact = vels .^ inc
        #println("v.fact = ",vfact)
        hts .*= vfact
        hts = map(x->max(20,min(20000,x)),hts)
        #println(hts)
        amps .*= sqrt.(oldhts ./ hts)
        amps = map(x->max(0.1,min(10.0,x)),amps)
        t += inc
    end
    println("end: ",tab[end][2:end])
    (tab,amptab)
end

midiToFreq(m::Int64) = 220*2^((m-57)/12.0)

function interpolate_wave(lst::Array{Array{Float64,1},1},
                     amplst::Array{Array{Float64,1},1})
    len = length(lst[1]) - 1
    lists = [ Array{Float64,1}[] for _=1:len ]
    amplists = [ Array{Float64,1}[] for _=1:len ]
    rem = lst
    remamp = amplst
    while length(rem)>0
        cur = rem[1]
        curamp = remamp[1]
        rem = rem[2:end]
        remamp = remamp[2:end]
        t = cur[1]
        hts = cur[2:end]
        amps = curamp[2:end]
        curvals = [ [t,h] for h in hts ]
        curamps = [ [t,h] for h in amps ]
        map(i -> push!(lists[i],curvals[i]),1:len)
        map(i -> push!(amplists[i],curamps[i]),1:len)
    end    
    #funs = map(interpolation_lin,lists)
    ampfuns = map(interpolation_lin,amplists)
    vfuns = map(z -> voice(z...),zip(ampfuns,lists))
    #return amplists[1],ampfuns[1],vfuns[1]
    return function(t::Float64)
        [ f(t) for f in vfuns ] |> sum
    end
end

clarinet(f) = [f*[1,3,5,7,9,11,13],[1,0.75,0.5,0.14,0.5,0.12,0.17]]

cis5 = 440.0*2^(1/3)
cis4 = 220.0*2^(1/3)
clA3 = clarinet(220.0)
clE4 = clarinet(220.0*2^(7/12))
clCis5 = clarinet(cis5)
clCis4 = clarinet(cis4)
clstart = map(z->vcat(z...),zip(clA3,clCis5,clE4))

wplay(lst) = wavplay(lst,sample_rate)

function cycle(start::Array{Array{Float64,1},1},
               d::Float64,
               fcon::Float64,
               rcon::Float64,
               tcon::Float64)
    ev = evolve(start...,d,fcon,rcon,tcon);
    int = interpolate_wave(ev...)
    play(int,d),ev
end
    
function cycle(start::Array{Array{Float64,1},1},
               d::Float64,
               rcon::Float64,
               tcon::Float64)
    fun(c) = evscore(start[1],
                     evolve(start...,d,c,rcon,tcon)[1][end][2:end])
    fcon = findbest(fun,0.01)
    println("fcon = ",fcon)
    ev = evolve(start...,d,fcon,rcon,tcon)
    int = interpolate_wave(ev...)
    play(int,d),ev
end

function findbest(f,tol)
    t = 1.0
    d = 1.0
    bnd = false
    s = 0
    while d > tol || s == -1
        s = f(t)
        #println(s)
        bnd = bnd || s == -Inf
        d *= !bnd ? 2.0 : 0.5
        t += s == -Inf ? -d : d
        #println(bnd," ,",t)
        #sleep(3)
    end
    return t
end

## w is in octaves
function lfo(f::Float64,m::Float64,w::Float64,p::Float64=0.)::Function
    function(t::Float64)
        m + w * sin(2*pi*f*t + p)
    end
end

function waver(f::Float64,o::Float64,wr::Float64)
    lf = lfo(wr,0.,o)
    function(t::Float64)
        f * 2^lf(t)
    end
end

function spring_force(nats::Array{Float64,1},
                      nows::Array{Float64,1})::Array{Float64,1}
    nats1 = vcat([20],nats,[20000])
    nows1 = vcat([20],nows,[20000])
    nats2,nows2 = log.(nats1),log.(nows1)
    nats3 = nats2[2:end] .- nats2[1:end-1]
    nows3 = nows2[2:end] .- nows2[1:end-1]
    diffs = nats3 .- nows3
    diffs[1:end-1] .- diffs[2:end]
end

function sort_start(fs,as)
    z = collect(zip(fs,as))
    sort!(z,by=first)
    collect(zip(z...))
end

function makestart(st,rat::Array{Float64,1})
    len = length(rat)
    res = [Float64[],Float64[]]
    amprat = 0.95
    i,f,a = 1,st,2.0
    while f < 10000
        push!(res[1],f)
        push!(res[2],a)
        i = i==len ? 1 : i+1
        f *= rat[i]
        a *= amprat
    end
    res
end

# function evscore(start,last)
#     rats = last ./ start
#     rats = log.(rats) ./ log(2)
#     rats = 1200 * abs.(rats)
#     score = sum(rats)
#     rats = last ./ [ 20 for _ in rats ]
#     rats = log.(rats) ./ log(2)
#     clst = 1200 * min(abs.(rats)...)
#     score += clst==0.0 ? -Inf : clst
#     rats = last ./ [ 20000 for _ in rats ]
#     rats = log.(rats) ./ log(2)
#     clst = 1200 * min(abs.(rats)...)
#     score += clst==0.0 ? -Inf : clst
#     score
# end

function evscore(start,last)
    i = findfirst(x->x==20.0 || x==20000.0,last)
    if i != nothing
        return -Inf
    else
        rats = last ./ start
        rats = log.(rats) ./ log(2)
        i = findfirst(x->abs(x)<0.8,rats)
        return i == nothing ? 1 : -1
    end
end

function plotev(ev)
    data = array_transpose(ev[1])
    len = length(data)
    for i=2:len
        data[i] = map(x->log(2,x),data[i])
    end
    args = [data[1],data[2]]
    for i=3:len
        push!(args,data[1],data[i])
    end
    plot(args...)
end

harms = [ (i+1)/i for i=1:30 ]

# function rand_stereo(lst)
#     len = length(lst[1][1]) - 1
#     pans = rand(0:0.01:1,len)
#     amps = 
# end

#curpl,curev = cycle(clstart,60.0,0.5,1.7);
#startocts = [[ 100.0 * 2^i for i=0:7 ],[ 1.0 - 0.1*i for i=0:7 ]]
#curpl,curev = cycle(makestart(100.,15/8),60.,0.512,20.);
#curpl,curev = cycle(makestart(100.,5/2),30.,0.5256,4.,10.);plotev(curev)
#curpl,curev = cycle(makestart(100.,15/8),30.,0.5070,1.,20.);plotev(curev)
#curpl,curev = cycle(makestart(196,[3/2]),20.,9.7,2.,0.1);plotev(curev)
#curpl,curev = cycle(makestart(264,[15/8,16/9]),10.,9.,2.,0.1);plotev(curev)

;
