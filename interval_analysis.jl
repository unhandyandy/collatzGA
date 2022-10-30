

include("mydefs.jl")
using MIDI,Printf,LinearAlgebra
#dot =  LinearAlgebra.dot

function get_max_time(notes)
    convert(Int64,max(map(n -> n.position + n.duration,notes)...))
end

struct intinst
    low::Int64
    int::Int64
    pos::Int64
    dur::Number
end

function make_intinst(n1::Note,n2::Note)
    s1,d1,s2,d2 = map(Int64,(n1.position,n1.duration,n2.position,n2.duration))
    st = max(s1,s2)
    dr = min(s1 + d1,s2 + d2) - st
    p1,p2 = map(Int64,(n1.pitch,n2.pitch))
    int = abs(p1 - p2)
    return intinst(min(p1,p2),int,st,dr)
end

function normalize(totdur::Int64,int::intinst)
    newlow = int.low % 12
    newdur = int.dur / totdur
    return intinst(newlow,int.int,int.pos,newdur)
end

function extract_notes(mf::MIDIFile)
    res = MIDITrack()
    for t in mf.tracks
        addnotes!(res,getnotes(t))
    end
    return res
end

struct incr_data
    time::Int64
    start::Int64
    last::Int64
    ints::Array{intinst,1}
    current::Array{Note,1}
end

function incr_data(st::Int64)
    return incr_data(st,st,st,intinst[],Note[])
end

function add_note(d::incr_data,n::Note)
    t = n.position
    nts = filter(o -> o.position + o.duration >= t,d.current)
    newdata = map(o -> make_intinst(n,o),nts)
    return incr_data(t,
                     d.start,
                     max(d.last,n.position + n.duration),
                     vcat(d.ints,newdata...),
                     vcat(nts,n))
end

function extract_data(mf::MIDIFile)
    nts = extract_notes(mf)
    t = Int64(getnotes(nts)[1].position)
    id = incr_data(t)
    for n in getnotes(nts)
        id = add_note(id,n)
    end
    return id
end

get_length(d::incr_data) = d.last - d.start

function print_tab(tab,normdur=true)
    cols = ["m2", "M2", "m3", "M3", "P4", "TT", "P5", "m6", "M6", "m7", "M7"]
    rows = ["C", "C♯", "D", "E♭", "E", "F", "F♯", "G", "A♭", "A", "B♭", "B"]
    dmpr = Array{Any,2}(undef,13,12)
    for c=1:11
        dmpr[1,c+1] = cols[c]
        for r=1:12
            dmpr[r+1,1] = rows[r]
            if normdur
                dmpr[r+1,c+1] = @sprintf("%.4f",tab[r,c])
            else
                dmpr[r+1,c+1] = @sprintf("%d",tab[r,c])
            end
        end
    end
    display(dmpr)
end


function summarize_data(d::incr_data,normdur=true)
    if normdur
        totdur = get_length(d)
    else
        totdur = 1
    end
    norms = map(i -> normalize(totdur,i),d.ints)
    dm = zeros(12,11)
    for i in norms
        low,int = i.low,i.int
        if 0 < int < 12 
            dd = i.dur
            r = 1 + low
            c = int
            dm[r,c] += dd
        end
    end
    print_tab(dm)
    return dm
end

function cents(rat)
    return 1200 * log(2,rat)
end

function cents_to_ratio(c)
    return 2 ^ (c/1200)
end

syncom_r = 81/80
pytcom_r = 3^12/2^19
syncom_c = cents(syncom_r) 
pytcom_c = cents(pytcom_r)
p5_c = cents(3/2)
diesis = cents(2 / 1.25^3)
diesisG = cents( 1.2^4 / 2)

function circle_to_scale(ints)
    scale = Dict{Int64,Float64}()
    scale[0] = 0
    scale[7] = prev = ints[1] + p5_c
    scale[2] = prev += ints[2] + p5_c
    scale[9] = prev += ints[3] + p5_c
    scale[4] = prev += ints[4] + p5_c
    scale[11] = prev += ints[5] + p5_c
    scale[6] = prev += ints[6] + p5_c
    scale[1] = prev += ints[7] + p5_c
    scale[8] = prev += ints[8] + p5_c
    scale[3] = prev += ints[9] + p5_c
    scale[10] = prev += ints[10] + p5_c
    scale[5] = prev += ints[11] + p5_c
    scalelist = []
    for i=0:11
        push!(scalelist,scale[i] % 1200)
    end
    return scalelist
end

function export_scale(scl,fn,target="tbl",cfreq=261626)
    # cfreq = 261626
    # cfreq = 246761
    name = fn * "." * target
    path = "/home/dabrowsa/Music/scala-22-pc64-linux/" * target * "/" * name
    if target=="tbl"
        freqs = [ ]
        for i=1:12
            newfreq = cfreq * cents_to_ratio(scl[i])
            push!(freqs, Int64(round(newfreq))) 
        end
        allfreqs = ""
        for i=0:127
            ref = (i % 12) + 1
            octs = 5 - div(i,12)
            newfreq = freqs[ref] / 2.0^octs
            newfreq = Int64(round(newfreq))
            allfreqs *= string(newfreq) * "\n"
        end
        #println(allfreqs)
        println(fn)
        write(path,allfreqs)
    elseif target=="scl"
        name = last(split(fn,"/"))
        txt="!\n!\n" * name * "\n12\n!\n"
        for i=2:12
            txt *= string(scl[i]) * "\n"
        end
        txt *= "2/1\n"
        println(txt)
        write(path,txt)
    end
end

function norm_vec(v)
    l = length(v)
    m = sum(v) / l
    res = map(x -> x - m,v)
    mm = dot(res,res)
    return sqrt(l / mm) .* res
end

function norm_tab(tab)
    res = Array{Float64,1}[]
    for i=1:11
        push!(res, norm_vec(tab[:,i]))
    end
    return [ res[i][p] for p=1:12,i=1:11 ]
end

jt = [1,16/15,9/8,6/5,5/4,4/3,25/18,3/2,8/5,5/3,16/9,15/8]
jt_scl = map(cents,jt)


function make_scale_cross_table(scl)
    sclex = vcat(scl,1200 .+ scl)
    init = [ sclex[n+i] - sclex[n] - jt_scl[i + 1] for n=1:12,i=1:11 ]
    return init
end

function compare_vecs(v1,v2)
    return dot(v1,v2)/sqrt(dot(v1,v1)*dot(v2,v2))
end

function compare_mats(m1,m2)
    sq = x -> x*x
    return [ compare_vecs(m1[:,i],
                          (i!=4) ? m2[:,i] : abs.(m2[:,i]))
             for i=1:11 ]
end

function correlation(v1,v2,absQ=false)
    l = length(v1)
    if l != length(v2)
        println("vectors of unequal len")
        return
    end
    m1,m2 = sum(v1)/l,sum(v2)/l
    s1,s2 = v1 .- m1, v2 .- m2
    if absQ
        corr = compare_vecs(s1,abs.(s2))
    else
        corr = compare_vecs(s1,s2)
    end
    d = sqrt(1 - corr^2)
    p = d^(l - 3) / (2 * sqrt(pi)) * (d + (l - 2) * (1 - abs(corr)))
    return corr, p
end


function dabcirc2(x)
    y = (pytcom_c - 5*x) / 3
    return [ 0,0,-x,-x,-x,-x,-x,-y,0,-y,-y ]
end

function eckersly(x)
    y = (pytcom_c - 5*x) / 4
    return [ -x,-x,-x,-x,-x,0,0,0,-y,-y,-y ]
end

function beats_to_cents5(f,b)
    d = (3*f + b)/2
    return cents(d/f * 2/3), d, cents(d/f)
end
function beats_to_cents4(f,b)
    d = (4*f + b)/3
    return cents(d/f * 3/4), d, cents(d/f)
end

function beats_to_scale(btlst,stfrq=261.626)
    scl = zeros(12)
    s = 0  ## 0-based
    d = stfrq
    newdc = 0
    for _=1:11
        b = btlst[s+1]
        # if b == nothing; println("break"); break end
        if b == nothing; break end
        _,d,dc = beats_to_cents5(d,b)
        if d/stfrq >= 2 d /= 2 end
        newdc = (newdc + dc) % 1200
        s = (s + 7) % 12
        scl[s+1] = newdc
        #println("s: ",s,"b: ",b,", d: ",d,", new: ",newdc)
    end
    s = 5 ## 0-based
    d = stfrq
    newdc = 0
    for  _=1:11
        b = btlst[s+1]
        if b == nothing break end
        if s < 5; b *= 2 end
        #println(b)
        _,d,dc = beats_to_cents4(d,-b)
        if d/stfrq >= 2 d /= 2 end
        newdc = (newdc + dc) % 1200
        scl[s+1] = newdc
        s = (s + 5) % 12
        #println("s: ",s,"b: ",b,", d: ",d,", new: ",newdc)
    end
    return scl
end

function prom_maj(i)
    t = div(i-1,2)
    return (t + (1 - (i % 2)) * 3) % 12 + 1
end
        
function get_selfs(data,tab,modes="Mm",int="M3")
    selftab = Any["T"];
    cents = Any["C"]
    start = ('M' in modes) ? 1 : 2
    inc = 3 - length(modes)
    last = ('m' in modes) ? 24 : 23
    if int=="M3"
        for i=start:inc:last
            push!(selftab,data_wtc1[i][prom_maj(i),4])
            push!(cents,tab[prom_maj(i),4])
        end
    elseif int=="P5"
        for i=start:inc:last
            push!(selftab,data_wtc1[i][div(i+1,2),7])
            push!(cents,tab[div(i+1,2),7])
        end
    end        
    # println(selftab)
    # println(cents)
    rows = ["","C", "c","C♯","c♯", "D","d", "E♭","e♭", "E","e", "F","f", "F♯","f♯", "G","g", "A♭","a♭", "A","a", "B♭","b♭", "B","b"]
    if length(modes)==1
        if ('M' in modes)
            rows = ["","C", "C♯", "D", "E♭", "E", "F", "F♯", "G", "A♭", "A", "B♭", "B"]
        else
            rows = ["","c", "c♯", "d", "e♭", "e", "f", "f♯", "g", "a♭", "a", "b♭", "b"]
        end
    end
    tabpr = zip(rows,selftab,cents)
    println("corr, p: ",correlation(selftab[2:end],cents[2:end]))
    return collect(tabpr)
end
function get_selfs(data)
    selftab = Any[];
    for i=1:24
        push!(selftab,data_wtc1[i][prom_maj(i),4])
    end
    tabpr = zip(["C", "c","C♯","c♯", "D","d", "E♭","e♭", "E","e", "F","f", "F♯","f♯", "G","g", "A♭","a♭", "A","a", "B♭","b♭", "B","b"],
                selftab)
    return collect(tabpr)
end

function freqs_to_beats(f1::Float64,f2::Float64)
    h1,h2 = f1,f2
    while true
        if abs(h1 - h2) <= 20; break end
        if h1 < h2
            h1 += f1
        else
            h2 += f2
        end
    end
    return h2 - h1
end

function freqs_to_beats(fs::Array{Float64,1})
    beats = []
    len = length(fs)
    for i=1:len-1
        for j=i+1:len
            push!(beats,freqs_to_beats(fs[i],fs[j]))
        end
    end
    return beats
end

function freqs_to_beats(fs::Array{Float64,1},chord::Array{Int64,1})
    frqs = Float64[]
    for i in chord
        ind = (i % 12) + 1
        oct = div(i - 60,12)
        oct -= i<60 ? 1 : 0
        push!(frqs,fs[ind] * 2.0^oct)
    end
    return freqs_to_beats(frqs)
end

function beats_p(st,btl,p)
    cents = beats_to_scale(btl,st)
    c1,c2 = cents[p],cents[(p+7) % 12]
    f1,f2 = st * 2.0^(c1/1200), st * 2.0^(c2/1200)
    #println(f1,", ",f2)
    return freqs_to_beats(f1,f2)
end

function smooth_circle(crc::Array{Float64,1},thirds::Array{Float64,1})
    v1 = [1,-1,1,-1,1,-1,1,-1,1,-1,1,-1]
    v2 = [1,1,-1,-1,1,1,-1,-1,1,1,-1,-1]
    v3 = [1,-1,-1,1,1,-1,-1,1,1,-1,-1,1]
    lim = 8.0
    #fun(t,v,dir) = norm(-et_crc + v + t .* dir)
    #fun(t,v,dir) = max(abs.(et_crc + v + t .* dir)...)
    fun(t,v,dir) = beat_ratio_dist(thirds,v + t .* dir)
    start = crc
    #while lim > 0.1
    t = find_min_local(t -> fun(t,start,v1),-lim,lim,0.01)
    start = start + t .* v1
    t = find_min_local(t -> fun(t,start,v2),-lim,lim,0.01)
    start = start + t .* v2
    t = find_min_local(t -> fun(t,start,v3),-lim,lim,0.01)
    start = start + t .* v3
    #lim /= 2
    #end
    return start
end

function thirds_to_crcet(mat)
    res = zeros(Float64,12)
    thirds = zeros(Float64,12)
    for j=1:4
        k,m = j + 4,j + 8
        m1,m2 = mat[2,j],mat[1,j]
        m3 = -m1 - m2
        del = m1/2
        res[j] -= del
        res[k] += del
        del += m2
        res[m] += del
        thirds[j] = m1
        thirds[k] = m2
        thirds[m] = m3
    end
    crcet = [ res[i+1] - res[i] for i=1:11 ]
    push!(crcet,res[1] - res[12])
    # crcet = smooth_circle(crcet,thirds)
    crcet = fixed_point(c -> smooth_circle(c,thirds),
                        crcet,
                        (x,y) -> max(abs.(x-y)...),
                        0.1)
    return crcet
end

function thirds_to_scl(mat)
    crcet = thirds_to_crcet(mat)
    crc = crcet + et_crc
    scl = circle_to_scale(crc)
    return scl
end

third_to_fifth(c::Number,r::Int64) = 1200 * log(2,5/(2*r) + 3/2 - 2/r * 2^(1/3+c/1200)) - 700

function beat_ratio(mc::Number,dc::Number)::Float64
    mr = 2^(1/3 + mc/1200)
    dr = 2^(7/12 + dc/1200)
    return (5 - 4 * mr) / (2 * dr - 3)
end

function beat_ratio_dist(mc::Number,dc::Number)::Float64
    #lim = 6
    r = abs(beat_ratio(mc,dc))
    #if r < lim
    err = r - floor(r)
    return 1/r * 2.0 ^ min(err,1 - err)
    #else
     #   return 0.0
    #end
end

function beat_ratio_dist(mc::Array{Float64,1},dc::Array{Float64,1})::Float64
    ds = [ beat_ratio_dist(mc[i],dc[i]) for i=1:12 ]
    return sum(ds)
end
    
M3_et = - ones(2,4) * (400 - cents(5/4))
    M3_mat_01 = [2  0 -4 4;
                 0 -4  0 2]
    M3_mat_02 = [2  2 -4 2;
                 0 -3  0 2]
    M3_mat_03 = [2  2 -4 2;
                 0 -4  0 2]
    M3_mat_04 = [2  2 -5 1;
                 0 -4  1 2]
    M3_mat_05 = [2  2.5 -5 1;
                 0 -5    1 2]
    M3_mat_06 = [2  2.5 -5.50 1;
                 0 -5    1.25 2]
    M3_mat_07 = [2  2.5 -6.0 1;
                 0 -5    1.5 2]
    M3_mat_08 = [2  2.5 -6.5 1;
                 0 -5    1.75 2]
# broadwoodU = M3_et + [17.69 17.69 17.69 17.69;
#                       5.69 6.69 10.69 14.69]
broadwoodU = [ 4  4  4 4;
              -8 -7 -3 1]

P5_to_M3_mat = [ 1 1 1 1 0 0 0 0 0 0 0 0;
                 0 1 1 1 1 0 0 0 0 0 0 0;
                 0 0 1 1 1 1 0 0 0 0 0 0;
                 0 0 0 1 1 1 1 0 0 0 0 0;
                 0 0 0 0 1 1 1 1 0 0 0 0;
                 0 0 0 0 0 1 1 1 1 0 0 0;
                 0 0 0 0 0 0 1 1 1 1 0 0;
                 0 0 0 0 0 0 0 1 1 1 1 0;
                 0 0 0 0 0 0 0 0 1 1 1 1;
                 1 0 0 0 0 0 0 0 0 1 1 1;
                 1 1 0 0 0 0 0 0 0 0 1 1;
                 1 1 1 0 0 0 0 0 0 0 0 1 ]

beats09 = [-1,nothing,-1,-1,-1,0,0,0,-2,-2,-2,-2]
beats10 = [0,-2,-2,-2,-2,-2,nothing,-1,-1,-1,0,0]
beats11 = [0,-1,-1,-1,-1,-1,nothing,-2,-2,-2,0,0]
beats12 = [-2,nothing,-2,-2,-2,0,0,0,-1,-1,-1,-1]

et_scl = [ 100 * i for i=0:11 ]
et_crc = [ -pytcom_c/12 for _=1:12 ]

dab2_scl = circle_to_scale(dabcirc2(1/6 * pytcom_c))
dab3_scl = circle_to_scale(dabcirc2(1/5 * pytcom_c))
dab4_scl = circle_to_scale(dabcirc2(1/5 * syncom_c))
dab5_scl = circle_to_scale(dabcirc2(1/6 * syncom_c))

dab6_scl = [ 0,
             102.84423990334108,
             201.56980965731668,
             297.74876216792086,
             401.31415196404737,
             498.04499913461245,
             600.8892390379538,
             699.6148087919292,
             798.7448774472318,
             901.441980810682,
             998.7189386886016,
             1101.4126612120663 ]

dab7_scl = [ 0.0,           
             102.14605963264967,
             201.56980965731668,
             298.73302660704576,
             400.6159716933562,
             498.04499913461245,
             600.1910587672623,
             699.6148087919292,
             799.7274659830084,
             901.441980810682,
             998.7189386886016, 
             1100.7144809413749 ] 

odonn1_crc = [ -1/6 * pytcom_c,
               -1/6 * pytcom_c,
               -1/6 * pytcom_c,
               0,
               -1/6 * pytcom_c,
               0,
               -1/12 * pytcom_c,
               -1/12 * pytcom_c,
               -1/12 * pytcom_c,
               -1/12 * pytcom_c,
               0 ]
odonn1_scl = circle_to_scale(odonn1_crc)
eck5_scl = circle_to_scale(eckersly(5))

lehman2_crc = [ -1/6 * pytcom_c,
                -1/6 * pytcom_c,
                -1/6 * pytcom_c,
                -1/6 * pytcom_c,
                0,
                0,
                0,
                -1/12 * pytcom_c,
                -1/12 * pytcom_c,
                -1/12 * pytcom_c,
                1/12 * pytcom_c,
                -1/6 * pytcom_c ]

zapf_scl_ratios = [      1/1     
       132449/125000
           28/25    
       148693/125000
         1257/1000  
       166967/125000
       353531/250000
          187/125   
        24803/15625 
          419/250   
       222789/125000
        29461/15625 ]
zapf_scl = map(cents,zapf_scl_ratios)

dab9_freqs = [  246.761
                261.97 
                277.356
                293.717
                311.026
                329.681
                349.154
                370.142
                391.956
                415.034
                439.575
                465.538 ]

dab7_freqs = [  246.761           
                261.758           
                277.231           
                293.235           
                311.01            
                329.015           
                349.011           
                369.642           
                391.647           
                415.347
                439.353           
                466.015            ]

lehman2_freqs = [  246.761           
                   260.55            
                   276.355           
                   292.788           
                   309.498           
                   328.643
                   347.400
                   369.306           
                   390.825           
                   413.597
                   438.686
                   463.2   ]

et_freqs = [  261.6255653005986 
              277.1826309768721 
              293.6647679174076 
              311.1269837220809 
              329.6275569128699 
              349.2282314330039 
              369.9944227116344 
              391.99543598174927
              415.3046975799451 
              440.0             
              466.1637615180899 
              493.8833012561241 ]

prelude03_01 = [ 54,72,81 ]
prelude03_02 = [ 54,74,81 ]
prelude03_03 = [ 54,63,71 ]
prelude03_04 = [ 54,61,70 ]

prelude08_01 = [ 46,49,54 ]
prelude08_02 = [ 47,51,54 ]

    M3_01_scl = [ 0,101,200,300,400,505,600,704,802,899,996,1100 ]
    M3_01_crc_et = [4,-4,-1,1,0,0,1,1,-2,-4,9,-5]

    
;


#=

lsstr = read(`ls /home/dabrowsa/Music/midi/wtc1/`,String)
fns = filter(s -> match(r"r-pf..\.mid",s)!=nothing,split(lsstr,"\n"))
wtc1 = []
for f in fns
       push!(wtc1,readMIDIfile("/home/dabrowsa/Music/midi/wtc1/" * f))
end
data_wtc1 = map(m -> summarize_data(extract_data(m)),wtc1);
data_all_wtc1 = reduce(+,data_wtc1)
data_all_wtc1_normed = norm_tab(data_all_wtc1)


julia> print_tab(data_all_wtc1)
13×12 Array{Any,2}:
 #undef   "m2"      "M2"      "m3"      "M3"      "P4"      "TT"      "P5"      "m6"      "M6"      "m7"      "M7"    
    "C"   "0.0173"  "0.0831"  "0.4164"  "0.2814"  "0.3700"  "0.1787"  "0.3845"  "0.3641"  "0.3548"  "0.2093"  "0.0838"
    "C♯"  "0.0282"  "0.0957"  "0.3666"  "0.3189"  "0.2808"  "0.1735"  "0.4400"  "0.2843"  "0.3947"  "0.1420"  "0.0519"
    "D"   "0.0104"  "0.0722"  "0.4094"  "0.2806"  "0.2988"  "0.1920"  "0.3319"  "0.3444"  "0.4120"  "0.1607"  "0.0517"
    "E♭"  "0.0245"  "0.1328"  "0.4330"  "0.2568"  "0.4014"  "0.1870"  "0.4692"  "0.3227"  "0.4991"  "0.2214"  "0.0709"
    "E"   "0.0293"  "0.0847"  "0.3894"  "0.2648"  "0.2970"  "0.1571"  "0.2339"  "0.2637"  "0.4601"  "0.1424"  "0.0438"
    "F"   "0.0218"  "0.0767"  "0.5601"  "0.2827"  "0.3952"  "0.1718"  "0.3859"  "0.3646"  "0.4902"  "0.2022"  "0.0536"
    "F♯"  "0.0095"  "0.1227"  "0.5189"  "0.2894"  "0.2609"  "0.1916"  "0.3561"  "0.3330"  "0.4187"  "0.1644"  "0.0364"
    "G"   "0.0287"  "0.0803"  "0.4363"  "0.2117"  "0.3422"  "0.2020"  "0.3936"  "0.4068"  "0.4240"  "0.1591"  "0.0389"
    "A♭"  "0.0244"  "0.1066"  "0.3984"  "0.2757"  "0.3736"  "0.1918"  "0.3908"  "0.3475"  "0.4857"  "0.1912"  "0.0588"
    "A"   "0.0154"  "0.0880"  "0.4171"  "0.3375"  "0.4228"  "0.2064"  "0.3015"  "0.2915"  "0.4332"  "0.1915"  "0.0588"
    "B♭"  "0.0180"  "0.1182"  "0.4608"  "0.2665"  "0.4059"  "0.1759"  "0.4375"  "0.3887"  "0.3950"  "0.2341"  "0.0808"
    "B"   "0.0289"  "0.1056"  "0.4084"  "0.3310"  "0.3052"  "0.1708"  "0.3162"  "0.3174"  "0.4091"  "0.1598"  "0.0576"




julia> print_tab(data_all_wtc1_durs,false)
13×12 Array{Any,2}:
 #undef   "m2"    "M2"     "m3"     "M3"     "P4"     "TT"     "P5"     "m6"     "M6"     "m7"     "M7"  
    "C"   "2082"  "8187"   "45058"  "28434"  "40197"  "19049"  "39941"  "38994"  "35473"  "22385"  "8491"
    "C♯"  "2947"  "10994"  "37754"  "36585"  "30168"  "18595"  "46403"  "27928"  "42475"  "13824"  "5988"
    "D"   "1191"  "6786"   "42102"  "26681"  "27916"  "22535"  "29751"  "34214"  "42867"  "15951"  "5408"
    "E♭"  "2551"  "15393"  "51637"  "29373"  "44629"  "18918"  "59139"  "33981"  "54821"  "25061"  "9332"
    "E"   "2828"  "8317"   "38216"  "26298"  "26628"  "14952"  "22573"  "26063"  "44360"  "13627"  "4530"
    "F"   "2556"  "8035"   "64335"  "26630"  "46411"  "18699"  "41524"  "40797"  "54250"  "22518"  "5088"
    "F♯"  "1225"  "12945"  "47825"  "34083"  "27858"  "20449"  "35741"  "30531"  "45834"  "14812"  "3965"
    "G"   "3410"  "7302"   "49136"  "20830"  "35003"  "20780"  "38359"  "47815"  "42559"  "16554"  "3896"
    "A♭"  "2600"  "11536"  "41207"  "31821"  "39577"  "20309"  "42422"  "32818"  "54345"  "19255"  "6657"
    "A"   "1390"  "8383"   "43220"  "31338"  "36392"  "20615"  "27690"  "29965"  "41575"  "16852"  "5921"
    "B♭"  "1921"  "13001"  "54766"  "26978"  "51391"  "17430"  "50193"  "41751"  "43861"  "25563"  "8818"
    "B"   "2868"  "11738"  "39817"  "34401"  "28737"  "17315"  "29630"  "31227"  "43769"  "14938"  "4790"

=#

