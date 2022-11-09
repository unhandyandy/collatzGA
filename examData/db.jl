
using CSV,DataFrames,Statistics,Interpolations,LinearAlgebra

## registrar format
function makeNameDict(csvf)
    len = length(csvf)
    d = Dict{String, Dict{String,String}}()
    for i=2:len
        #println(i)
        entry = csvf[i][3]
        if !ismissing(entry)
            entry = split(entry,"@")[1]
            un = lowercase(entry)
            idstr = lpad(csvf[i][1],10,"0")
            d[un] = Dict("id"=>idstr,"name"=>csvf[i][2])
        end
    end
    d
end

function makeIDDict(csvf)
    len = length(csvf)
    d = Dict{String, Dict{Union{String,Missing},Union{String,Missing}}}()
    for i=1:len
        id = lpad(csvf[i][1],10,"0")
        if !ismissing(id)
            un = csvf[i][3]
            if !ismissing(un)
                un = split(un,"@")[1]
                un = lowercase(un)
            end
            name = csvf[i][2]
            d[id] = Dict("username"=>un,"name"=>name)
        end
    end
    d
end



checkCode(c,codes) = occursin(c,codes)

function version(row,cD,codes)
    c = row[14+cD]
    if ismissing(c) || c==""
        missing
    else
        checkCode(ec(row),codes) ? (c % 4) + 1 : missing
    end
end

row(i) = df[i,:]

readEntry(x) = ismissing(x) ? "" : x

username(row)=lowercase(join([readEntry(row[i]) for i=1:8]))

ec(row)=join([readEntry(row[i]) for i=14:19], " ")

answer(v,n,anskeys) = ismissing(v) ? missing : anskeys[v][n+1]

function response(row,n)
    r = row[n+19]
    ismissing(r) ? missing :
        length(r)==1 ? r : missing
end

function scoreone(row,v,n,anskeys)
    @bp
    a = ismissing(v) ? missing : uppercase(answer(v,n,anskeys))
    r = response(row,n)
    ismissing(v) ? missing :
        ismissing(r) ? 0 :
        occursin(r,a) ? 1 : 0
end

function scorerow(row,nP,v,anskeys)
    sum([scoreone(row,v,n,anskeys) for n=1:nP])
end

function classmean(ss)
    sl = filter(x->!ismissing(x),ss)
    mean(sl)
end

function extractData(frm,term,exam,codeDigit,numberPrompts,codes,anskeys,nameDict,vmap)
    len = size(frm)[1]
    nf = DataFrame()
    nf.term = String7[]
    nf.exam = String15[]
    nf.numberPrompts = Int16[]
    # nf.realname=Union{String,Missing}[]
    nf.IUIDnumber=Union{String15,Missing}[]
    nf.username=Union{String15,Missing}[]
    nf.version=Union{Int8,Missing}[]
    atype = Union{String1,Missing}
    stype = Union{Int8,Missing}
    for n=1:numberPrompts
        nf.n1=atype[]
        qname = string("Q",n)
        rename!(nf,:n1=>qname)
        nf.n2=stype[]
        sname = string("S",n)
        rename!(nf,:n2=>sname)
    end
    nf.score=Union{Int8,Missing}[]
    for i=1:len
        row = frm[i,:]
        uname = username(row)
        iuiddict = get(nameDict,uname,missing)
        iuid = ismissing(iuiddict) ? missing : iuiddict["id"]
        v = version(row,codeDigit,codes)
        newrow = [term exam numberPrompts iuid uname v]
        quests = Matrix{Any}(repeat([missing],1,numberPrompts*2))
        for n=1:numberPrompts
            res = response(row,n)
            scr = scoreone(row,v,n,anskeys)
            j = ismissing(v) ? 2*n : vmap(v,n)*2
            i = j - 1
            quests[i:j] = [res scr]
        end
        newrow = [newrow quests]
        tot = scorerow(row,numberPrompts,v,anskeys)
        newrow = [newrow tot]
        push!(nf,newrow)
    end
    nf
end

function extractData2(frm,term,exam,codeDigit,numberPrompts,codes,anskeys,nameDict,vmap)
    numberVersions = 4
    numCols = 6 + numberPrompts*(numberVersions*2 + 1) + 1
    len = size(frm)[1]
    nf = DataFrame()
    nf.term = String7[]
    nf.exam = String15[]
    nf.numberPrompts = Int16[]
    # nf.realname=Union{String,Missing}[]
    nf.IUIDnumber=Union{String15,Missing}[]
    nf.username=Union{String15,Missing}[]
    nf.version=Union{Int8,Missing}[]
    atype = Union{String1,Missing}
    stype = Union{Int8,Missing}
    for n=1:numberPrompts
        nf.n1=Union{String,Missing}[]
        cname = string("problem",n)
        rename!(nf,:n1=>cname)
        for v=1:numberVersions
            nf.n1=atype[]
            qname = string("V",v,"Q",n)
            rename!(nf,:n1=>qname)
            nf.n2=stype[]
            sname = string("V",v,"S",n)
            rename!(nf,:n2=>sname)
        end
        # nf.n1 = atype[]
        # qname = string("V?Q",n)
        # rename!(nf,:n1=>qname)
    end
    nf.score=Union{Int8,Missing}[]
    probs = getProbNames()
    for i=1:len
        row = frm[i,:]
        j,_ = size(nf)
        j += 1
        uname = username(row)
        iuiddict = get(nameDict,uname,missing)
        iuid = ismissing(iuiddict) ? missing : iuiddict["id"]
        v = version(row,codeDigit,codes)
        newrow = [[term exam numberPrompts iuid uname v] repeat([missing],numCols - 6)']
        push!(nf,newrow)
        if ismissing(v)
            continue
        end
        for n=1:numberPrompts
            pccol = string("problem",n)
            nf[j,pccol] = probs[n]
            vn = vmap(v,n)
            colq = Symbol(string("V",v,"Q",vn))
            cols = Symbol(string("V",v,"S",vn))
            res = response(row,n)
            scr = scoreone(row,v,n,anskeys)
            nf[j,colq] = res
            nf[j,cols] = scr
        end
        tot = scorerow(row,numberPrompts,v,anskeys)
        nf[j,:score] = tot
    end
    nf
end

function compareVers(d)
    nP = d[1,:numberPrompts]
    numSamps,_ = size(d)
    sig4 = p->8*sqrt(p*(1-p)/numSamps)
    ## cols with scores for each prompt
    colnums = collect(8:2:6+2*nP)
    ## create vector of success rates for all prompts
    ps = mapcols(c->mean(skipmissing(c)),select(d,colnums))
    ps = ["total"; values(ps[1,:])...]
    sigs = sig4.(ps[2:end])
    sigs = ["4sig"; sigs]
    ## compare versions
    vergrps = groupby(d,:version)
    ## len = size(d)[2]
    cpVers = combine(vergrps, colnums .=> mean)
    ## transpose dataframe:
    cpVers = DataFrame([[names(cpVers)]; collect.(eachrow(cpVers))], [:column; Symbol.(axes(cpVers, 1))])
    cpVers.total = ps
    cpVers[!,"4sig"] = sigs
    ## exam means
    show(cpVers[2:end,[collect(2:5)...,7,8]])
    println()
    combine(cpVers[2:end,:], [collect(2:5)...,7] .=> mean)
end

getAve(d) = skipmissing(d[!,end]) |> mean

getIUIDs() =  CSV.File("iuids.csv",select=[2,14,15])

## codeDigit is 0-based
function getData(;dir,term,exam,codeDigit,numberPrompts)
    cd(dir)
    if isfile("all_reads.csv")
        run(`rm all_reads.csv`)
    end
    run(`bash -c 'cat  *_reads.csv > all_reads.csv'`)
    allstr=read("all_reads.csv",String)
    allstr = replace(allstr,"BLANK []"=>"","BLANK"=>"")
    csv = CSV.File(IOBuffer(allstr),header=0)
    global df = DataFrame(csv)

    codes=read("codes.csv",String)

    anskeys=CSV.File("keys.csv",header=0)
    iuids = 0
    if !isfile("iuids.csv")
        if isfile("iuids.XLS")
            run(`ssconvert iuids.XLS iuids.cvs`)
            sleep(2)
        else
            iuids = -1
        end
    end
    iuids = iuids==-1 ? -1 : getIUIDs()

    nameDict = iuids==-1 ? Dict() : makeNameDict(iuids)
    vmap = (v,n) -> n
    if isfile("versionMap.csv")
        vermap = CSV.File("versionMap.csv",header=0)
        vmap = (v,n) -> vermap[v][n]
    end

    data = extractData(df,term,exam,codeDigit,numberPrompts,codes,anskeys,nameDict,vmap);
    data2 = extractData2(df,term,exam,codeDigit,numberPrompts,codes,anskeys,nameDict,vmap);
    CSV.write(term * "_" * exam * ".csv",data2)
    data
end

function studentScore(uname,data)
    un = data[uname.== data.username,end]
    length(un)==0 ? missing : un[1]
end

function studentProgress(d1,c1,d2,c2)
    len,_ = size(d2)
    cf1 = makeCurveFun(c1)
    cf2 = makeCurveFun(c2)
    nD = DataFrame()
    nD.student = String[]
    nD.S1 = Union{Int8,Missing}[]
    nD.C1 = Union{Float32,Missing}[]
    nD.G1 = Union{String1,Missing}[]
    nD.S2 = Union{Int8,Missing}[]
    nD.C2 = Union{Float32,Missing}[]
    nD.G2 = Union{String1,Missing}[]
    for s=1:len
        uname = d2[s,:username]
        if !ismissing(uname)
            s1 = studentScore(uname,d1)
            s2 = d2[s,:score]
            if !ismissing(s1) & !ismissing(s2) 
                c1 = cf1(s1)
                g1 = gradeFromCurve(c1)
                c2= cf2(s2)
                g2 = gradeFromCurve(c2)
                push!(nD,[uname,s1,c1,g1,s2,c2,g2])
            end
        end
    end
    nD
end

function removeMissing(df)
    df[.!ismissing.(df.S1) .& .!ismissing.(df.S2),:]
end

function regressMat(mat::Matrix)
    len,_ = size(mat)
    m1 = mat[:,1]
    m = hcat(ones(len),m1)
    m2 = mat[:,2]
    proj = (m'*m)^-1*m'
    predict = m*proj*m2
    v1 = var(m2)
    v2 = var(m2 - predict)
    err = sqrt(v2/v1)
    b,a = proj*m2
    r = cor(m1,m2)
    println("aX + b")
    println("a = ",a)
    println("b = ",b)
    println("err = ",err)
    println("r = ",r)
    (a,b,r,err)
end

function regress(mat,vec)
    minv = (mat'*mat)^-1
    coeffs = minv*mat'*vec
    preds = mat*coeffs
    rel = cor(preds,vec)
    coeffs,rel
end

function regressProg(prog::DataFrame)
    gradeG = groupby(prog,:G1)
    for k in Base.keys(gradeG.keymap)
        g = k[1]
        i = gradeG.keymap[k]
        mat = [gradeG[i].C1 gradeG[i].C2]
        println("group: ",g)
        regressMat(mat)
    end
    gradeG
end

function makeCurveFun(rawscores)
    curve = linear_interpolation(rawscores, [60,70,80,90,100]; extrapolation_bc=Line())
    f(x) = ismissing(x) ? missing : curve(x)
    f
end

function gradeFromCurve(c)
    ismissing(c) ? missing :
        c>=90 ? "A" :
        c>=80 ? "B" :
        c>=70 ? "C" :
        c>=60 ? "D" : "F"
end

# function addGrades(df,curvevec)
#     cfun = makeCurveFun(curvevec)
#     ss = df.score
#     cs = cfun.(ss)
#     gs = gradeFromCurve.(cs)
#     hcat(df,DataFrame("curve"=>cs,"grade"=>gs))
# end
    
function regressCols(df,col1,col2)
    len,_ = size(df)
    mat = rand(0,2)
    for r=1:len
        x1 = df[r,col1]
        x2 = df[r,col2]
        if !ismissing(x1) & !ismissing(x2)
            mat = [mat;x1 x2]
        end
    end
    regressMat(mat)
end

function entropy(df,col)
    cv = skipmissing(df[!,col])
    a = minimum(cv)
    b = maximum(cv)
    len = b-a+1
    counts = zeros(1,len)
    for i=1:len
        v = a+i-1
        counts[i] = count(cv .== v)
    end
    n = sum(counts)
    counts = counts / n
    ent = 0
    for i = 1:len
        ent -= counts[i]*log(2,counts[i])
    end
    (ent,counts...)
end

function calcEigMax(qmat)
    _,nP = size(qmat)
    vec = abs.(eigvecs(cor(qmat))[:,end])
    vec/sum(vec)*nP
end


function qStats(df)
    nP = df[1,:numberPrompts]
    # vv = sort([[entropyC(df,n)...] for n=1:nP])
    # reduce(hcat,vv)'
    qdf = DataFrame()
    # qdf.Q = Int8[]
    qdf.reliability = Float32[]
    qdf.cut = Int8[]
    qdf.info = Float32[]
    qdf.qinfo = Float32[]
    qdf.success = Float32[]     
    qdf.discr = Float32[]
    qdf.delTotal = Float32[]
    qdf.weights = Float32[]
    wts = calcEigMax(makeSmat(df))
    scrs = df.score
    sigma = std(skipmissing(scrs))
    for n=1:nP
        col = Symbol(string("S",n))
        c,e,enew = entropyC(df,n)
        r,_ = calcRelQ(df,n)
        _,_,p = entropy(df,col)
        ss = df[!,[col,:score]]
        ss = filter(r->!ismissing(r[col]) & !ismissing(r[:score]),ss)
        crl = cor(ss[!,col],ss[!,:score])
	del = crl*sigma*sqrt(1/p-1) - 1
        push!(qdf,[r,c,e,enew,p,crl,del,wts[n]])
    end
    qdf 
end

# function makeQMat(df)
#     nP = df[1,:numberPrompts]
#     nS,_ = size(df)
#     res = Matrix{Union{Float32,Missing}}(zeros(0,nS))
#     for i=1:nP
#         sym = Symbol(string("S",i))
#         res = vcat(res,df[!,sym]')
#     end
#     res
# end

function covQMat(mat)
    n,_ = size(mat)
    res = zeros(Float32,n,n)
    for i=1:n
        for j=1:i
            v1 = mat[i,:]
            v2 = mat[j,:]
            inds = ismissing.(v1) .| ismissing.(v2)
            inds = .!inds
            res[i,j]=res[j,i]=v1[inds]'*v2[inds]
        end
    end
    res
end

removeMissing(df) = df[.!ismissing.(df.score),:]

function checkSameSignCol(mat)
    r,c = size(mat)
    for i=1:c
        if sum(abs.(mat[:,i]))==abs(sum(mat[:,i]))
            println(i)
        end
    end
end

function covApprox(mat,vals,n)
    len,_ = size(mat)
    res = zeros(len,len)
    for i=len - n + 1:len
        res += vals[i]*mat[:,i]*mat[:,i]'
    end
    res
end

function entP(p::Number)
    if in([0,1])(p)
        0
    else
        -p*log(2,p)-(1-p)*log(2,1-p)
    end
end

function entP(ps::Vector)
    psn = ps[ps .!= 0]
    length(psn)>0 ? - sum(psn.*log2.(psn)) : 0
end

function entP(ps::Matrix)
    nr,nc = size(ps)
    ent = 0
    for r=1:nr
        row = ps[r,:]
        p = sum(row)
        if p>0
            ent += p*entP(row/p)
        end
    end
    ent
end

logit(p) = log(p/(1-p))
function invLogit(x)
    ex = exp(x)
    ex/(1+ex)
end

logit2(p) = log2(p/(1-p))

## amount of info obtained about score on question q
## by knowledge of whether total score >= c
function entropyC(df,q,c)
    len,_ = size(df)
    cntrs = zeros(2,2)
    col = Symbol(string("S",q))
    ss = df[!,col]
    scores = df[!,:score]
    for i=1:len
        s,t = ss[i],scores[i]
        if !ismissing(t)
            if t<c
                cntrs[1,s+1] += 1
            else
                cntrs[2,s+1] += 1
            end
        end
    end
    m = sum(cntrs[1,:])
    n = sum(cntrs[2,:])
    l = sum(cntrs[:,2])
    p,q,r,s = cntrs[1,1]/m,cntrs[2,2]/n,l/len,m/len
    e1,e2 = entP(p),entP(q)
    e0 = entP(r)
    enew =  s*e1 + (1-s)*e2
    e0 - enew, enew
end

function entropyC(df,q::Number)
    bestc = -1
    beste = -1
    bestnew = -1
    num = df[1,:numberPrompts]
    for c = 1:num
        e,enew = entropyC(df,q,c)
        if e>beste
            beste = e
            bestc = c
	    bestnew = enew
        end
    end
    (bestc,beste,bestnew)
end

function cliqueMatches(corMat,list,cor)
    n,_ = size(corMat)
    if length(list) == 0
        return [i for i=1:n]
    end
    len = length(list)
    mat = corMat[:,list].>=cor
    lbits = reduce(.&,[mat[:,i] for i=1:len])
    findall(in(1),lbits)
end

function checkClique(corMat,list,cor)
    rows = cliqueMatches(corMat,list,cor)
    res = true
    for r in list
        res &= in(rows)(r)
    end
    res
end

function findCliques(corMat,list,cor)
    res = []
    rows = cliqueMatches(corMat,list,cor)
    if length(list) == length(rows)
        push!(res,sort(rows))
    else
        for r in rows
            if !in(list)(r)
                res = vcat(res,findCliques(corMat,[list...,r],cor))
            end
        end
    end
    unique(res)
end

findCliques(corMat,cor) = findCliques(corMat,[],cor)

function entropyExcessQ2(df,q1,q2)
    num = df[1,:numberPrompts]
    len,_ = size(df)
    sym1 = Symbol(string("S",q1))
    sym2 = Symbol(string("S",q2))
    l1 = df[!,sym1]
    l2 = df[!,sym2]
    cntr = zeros(2,2)
    for i=1:len
        s1,s2 = l1[i],l2[i]
        if !ismissing(s1) & !ismissing(s2)
            cntr[s1+1,s2+1] += 1
        end
    end
    pmat = cntr/sum(cntr)
    p = sum(pmat[:,1])
    q = sum(pmat[1,:])
    entP(p) + entP(q) - entP(pmat) 
end

function makeEntExcMat(df)
    num = df[1,:numberPrompts]
    [entropyExcessQ2(df,i,j) for i=1:num,j=1:num]
end

function makeSCols(qnums)
    cols = Symbol[]
    nq = length(qnums)
    for q in qnums
        push!(cols, Symbol(string("S",q)))
    end
    cols
end

function calcQEntropy(df,qs)
    cols = Symbol[]
    nq = length(qs)
    for q in qs
        push!(cols, Symbol(string("S",q)))
    end
    nr = 2^nq
    nP = df[1,:numberPrompts]
    cnts = zeros(Float32,nr,nP+1)
    dft = removeMissing(df)
    nS,_ =size(dft) 
    pows = [2^i for i=0:nq-1]
    for row in eachrow(dft)
        bits = Vector(row[cols])
        r = sum(pows .* bits) + 1
        cnts[r,row[:score]+1] += 1
    end
    cnts /= nS
    entP(cnts)
end

function findOptEntQs(df,tol=1)
    nP = df[1,:numberPrompts]
    qs = Int8[]
    ent = Inf
    cands = collect(1:nP)
    while ent>tol
        ent,ind = findmin(x->calcQEntropy(df,[qs...,x]),cands)
        q = cands[ind]
        push!(qs,q)
        cands = filter(!in(qs),cands)
    end
    qs,ent
end

function calcQSD(df,qs)
    cols = Symbol[]
    nq = length(qs)
    for q in qs
        push!(cols, Symbol(string("S",q)))
    end
    nP = df[1,:numberPrompts]
    dft = removeMissing(df)
    nS,_ =size(dft)
    unit = ones(Float32,nS)
    mat = map(x->dft[!,x],cols)
    mat = Matrix{Float32}(reduce(hcat,[unit,mat...]))
    scores = dft.score
    spectrum = mat'*scores
    mcov = mat'*mat
    minv = mcov^-1
    coeffs = minv*spectrum
    preds = mat*coeffs
    errs = scores - preds
    std(errs),coeffs
end

function makeSmat(df)
    dft = removeMissing(df)
    nP = dft[1,:numberPrompts]
    cols = makeSCols(collect(1:nP))
    vec = [dft[!,c] for c in cols]
    reduce(hcat,vec)
end


function makeSvLPTable(df,q)
    col = Symbol(string("S",q))
    nP = df[1,:numberPrompts]
    success = zeros(nP+1)
    studs = zeros(nP+1)
    for r in eachrow(df)
        s,t = r[:score],r[col]
        if !ismissing(s)
            studs[s+1] += 1
            success[s+1] += t
        end
    end
    nums = collect(0:nP)
    tab = success ./ studs
    l = .!isnan.(tab)
    [nums[l] tab[l] studs[l]]
end


function calcRelQ(df,q)
    svp = makeSvLPTable(df,q)
    len,_ = size(svp)
    mat = [svp[:,1] ones(len)]
    vec = svp[:,2]
    coeffs,rel = regress(mat,vec)
    rel,coeffs
end

S4222_midterm_fractions = [2,4,6,7,10,14,15,16,18,24]
S4222_midterm_PCPr = [1,2,8,9,11,13,18,19,21,23]
S4222_midterm_SetTheory = [3,22,25]

function makeSqRoot(df::DataFrame)
    smat = makeSmat(df)
    scor = cor(smat)
    evals,evecs = eigen(scor)
    for i=1:length(evals)
       evecs[:,i] *= sqrt(evals[i])
    end
    -evecs'
end

function makeSqRoot(mat::Matrix)
    evals,evecs = eigen(mat)
    for i=1:length(evals)
        evecs[:,i] *= sqrt(evals[i])
    end
    -evecs'
end
 

sortVec(v) = [r for r in eachrow([v 1:length(v)])] |> sort

charToInt(c) = Int(c[1]) - 64
intToChar(n) = Char(n + 64)

function countAnswers(df,q,nA,nV=4)
    acol = Symbol(string("Q",q))
    counter = zeros(Int,nA,nV)
    for row in eachrow(df)
        v = row[:version]
        a = row[acol]
        if ismissing(v) || ismissing(a); continue end
        r = charToInt(a)
        counter[r,v] += 1
    end
    [intToChar.(1:nA) counter]
end

function getProbNames()
    if isfile("probs.txt")
        probStr = read("probs.txt",String)
        split(probStr,"\n")
    else
        repeat([""],50)
    end
end

midtermBasis = Dict(1=>"math aptitude",2=>"intuition",3=>"algebra")

## convert -delete 0-1 -density 300 a.pdf[2448X2200+0+0] -crop 100%x50% -alpha opaque out.png
## ristretto prob-{11,8,21,12}.png&
## feh --zoom 35 -g 800x400 out-{11,8,21,9,2,15,5}.png&

S4222_final_PCPr = [6,7,8,11,12]
#1: general aptitude
#2: permutations
#3: blocks
#4: P v. C v. power
#5: set theory


function calcB(a,r)
    c2 = 2*r*a
    c3 = a^2 - 1
    (-c2 + sqrt(c2^2-4*c3))/2
end

function formCorrs2(vec)
    v = invLogit.(vec)
    r = v[1]
    as = v[2:end]
    cormat = [1 r; r 1]
    bs = map(a -> calcB(a,r),as)
    coeffmat = [as'; bs']
    cormat, coeffmat
end

function formCorrs(vec,nvars;renormalize=true)
    len1 = Int(round(nvars*(nvars-1)/2))
    v = invLogit.(vec)
    rs = v[1:len1]
    as = v[len1+1:end]
    cormat = Matrix{Float64}(I(nvars))
    for r=nvars:-1:2
        for c = r-1:-1:1
            cormat[r,c] = cormat[c,r] = pop!(rs)
        end
    end
    len2 = length(as)
    nc = Int(round(len2/nvars))
    coeffmat = reshape(as,nc,nvars)' |> Matrix
    if renormalize
        for c=1:nc
            v = coeffmat[:,c]
            m = v'*cormat*v
            coeffmat[:,c] /= sqrt(m)
        end
    end
    cormat, coeffmat, renormalize ? formVec(cormat, coeffmat) : vec
end

mutable struct MatMask
    n::Int8
    mat::Matrix{Bool}
    len::Int8
end

function MatMask(n)
    empty = zeros(n,n).==1
    MatMask(n,empty,0)
end

function toggle(mm::MatMask,i,j)
    old = mm.mat[i,j]
    mm.mat[i,j] = mm.mat[j,i] = !old
    mm.len += mm.mat[i,j]==1 ? 1 : -1
    mm
end

function formVec(crm,cfm)
    nvars,_ = size(crm)
    nr,nc = size(cfm)
    res = Float64[]
    for r=2:nvars
        for c=1:r-1
            push!(res,logit(crm[r,c]))
        end
    end
    for r=1:nr
        for c=1:nc
            push!(res,logit(cfm[r,c]))
        end
    end
    res
end


function calcError(target,cormat::Matrix,coeffmat::Matrix)
    nq,_ = size(target)
    diff = target - coeffmat'*cormat*coeffmat
    sum(diff.^2)/(nq^2)
end

function calcError(target,vec::Array,nvars::Int)
    crm,cfm,newvec = formCorrs(vec,nvars,renormalize=false)
    calcError(target,crm,cfm),newvec
end

function calcGradient(target,vec,nvars)
    del = 0.01
    len = length(vec)
    res = zeros(len)
    err0,_ = calcError(target,vec,nvars)
    for i=1:len
        v = copy(vec)
        v[i] += del
        err,_ = calcError(target,v,nvars)
        res[i] = (err - err0)/del
    end
    res
end

function improveOne(target,vec,nvars)
    err,v = calcError(target,vec,nvars)
    g = calcGradient(target,v,nvars)
    m = norm(g)
    v - err/(2*m)*g,err
end

function estimateParams(target,v0,tol,nvars)
    v = copy(v0)
    olderr = Inf
    newerr = 0
    cnt = 0
    while abs(olderr - newerr)>tol || cnt<40
        cnt += 1
        olderr = newerr
        v,newerr = improveOne(target,v,nvars)
        println(newerr)
    end
    formCorrs(v,nvars,renormalize=false)[1:2]...,newerr
end
function estimateParams(target,v0,nvars)
    v = copy(v0)
    olderr = Inf
    newerr = 1000000
    while newerr<olderr
        olderr = newerr
        v,newerr = improveOne(target,v,nvars)
        # println(newerr)
    end
    formCorrs(v,nvars,renormalize=false)[1:2]...,newerr
end

function makeParamVec(nq,nv)
    l1 = Int(nv*(nv-1)/2)
    l2 = nq*nv
    logit.(rand(l1+l2))
end

function randFit(target,nq,nv)
    v = makeParamVec(nq,nv)
    res = estimateParams(target,v,nv)
    display(res[1])
    display(res[2])
    display(res[3])
    res
end

makeProjRows(mat) = mat'*(mat*mat')^-1*mat

## add one variable 
function formCorrs3(vec,precrm,precfm)
    nv,_ = size(precrm)
    _,nq = size(precfm)
    as = invLogit.(vec[1:nv])
    bs = invLogit.(vec[nv+1:end])
    crm = [1 as';as precrm]
    cfm = [bs';precfm]
    crm,cfm
end

function calcError3(target,vec::Array,precrm,precfm)
    crm,cfm = formCorrs3(vec,precrm,precfm)
    calcError(target,crm,cfm)
end

function calcGradient3(target,vec,precrm,precfm)
    del = 0.01
    len = length(vec)
    res = zeros(len)
    err0 = calcError3(target,vec,precrm,precfm)
    for i=1:len
        v = copy(vec)
        v[i] += del
        err = calcError3(target,v,precrm,precfm)
        res[i] = (err - err0)/del
    end
    res
end

function improveOne3(target,vec,precrm,precfm)
    crm,cfm = formCorrs3(vec,precrm,precfm)
    err = calcError(target,crm,cfm)
    g = calcGradient3(target,vec,precrm,precfm)
    m = norm(g)
    vec - err/(2*m)*g,err
end

function estimateParams3(target,v0,precrm,precfm)
    v = copy(v0)
    olderr = Inf
    newerr = 1000000
    while newerr<olderr
        olderr = newerr
        v,newerr = improveOne3(target,v,precrm,precfm)
        # println(newerr)
    end
    formCorrs3(v,precrm,precfm)...,newerr
end

function perturb(mat)
    newmat = copy(mat)
    nr,nc = size(mat)
    for r = 1:nr
        c = rand(1:nc)
        prev = newmat[r,c]
        newmat[r,c] = prev==1 ? 0 : 1
    end
    newmat
end

function invert_dict(d)
    Dict(values(d) .=> keys(d))
end

function all_answers(df,q,permmap)
    col = Symbol(string("Q",q))
    dfver = groupby(df,[:version,col],skipmissing=true,sort=true)
    keyvec = keys(dfver.keymap)
    # revmap = invert_dict(keymap)
    ## len = dfver.ngroups
    pm = permmap[q]
    nv,na = size(pm)
    res = zeros(Int32,nv,na)
    ## res = [ [revmap[i] size(dfver[i])[1]] for i=1:len]
    for r=1:nv
        function f(r,k)
            ans = string(intToChar(k))
            key = (r,ans)
            if key in keyvec
                size(dfver[key])[1]
            else
                0
            end
        end
        res[r,:] = map(k->f(r,k),pm[r,:])
    end
    ##reduce(vcat,res)
    header = reshape(string.(intToChar.(1:na)),(1,na))
    display([header; res])
    res
end

function getCurve()
    curve = DataFrame(CSV.File("curveRef.csv"))
    function rawToCurve(r)
        if ismissing(r)
            return missing
        else
            c = curve.Curve[curve.Raw.==r]
            return c[1]
        end
    end
end

function scoreToLetter(s)
    if ismissing(s)
        return missing
    end
    return s<60 ? "F" :
        s<70 ? "D" :
        s<80 ? "C" :
        s<90 ? "B" : "A"
end

function addGrades(df)
    rtc = getCurve()
    df.curve = rtc.(df.score)
    df.grade = scoreToLetter.(df.curve)
end

function readPerms(enm)
    fnm = string(enm,"_perms_table.lua")
    permvec = readlines(fnm)
    start = false
    res = Any[]
    for l in permvec
        # println(res)
        start = start || l=="-- Table: {2}"
        if start
            if l=="{"
                push!(res,Int8[])
            elseif l=="},"
                continue
            else
                m = match(r"   ([^,]),",l)
                if !isnothing(m)
                    push!(res[end],parse(Int8,m[1]))
                end
            end
        end
    end
    res
end

function perm_divide(p1,p0)
    len = length(p0)
    res = zeros(Int8,len)
    for i=1:len
        res[i] = findfirst(in(p1[i]),p0)
    end
    res                    
end

function combinePerms(p0,ps...)
    nv = 1 + length(ps)
    len = length(p0)
    res = Matrix{Int8}[]
    for i=1:len
        na = length(p0[i])
        push!(res,zeros(nv,na))
        res[i][1,:] = collect(1:na)
        for j=2:nv
            perm = ps[j-1][i]
            res[i][j,:] = perm_divide(p0[i],perm)
        end
    end
    res
end

function getSemGrades(fn,course,term)
    grades = DataFrame(CSV.File(fn))
    grades.UNIV_ID = lpad.(grades.PRSN_UNIV_ID,10,"0")
    select!(grades,Not(:PRSN_UNIV_ID))
    grades = grades[grades.COURSE.==course .&& grades.TERM.==term,:]
    unique(grades)
end
    
function findInDF(df,col,val)
    ## println(val)
    if ismissing(val)
        missing
    else
        fa = findall(==(val), df[:,col])
        ## println(fa)
        if length(fa)==1
            df[only(fa), :]
        else
            missing
        end
    end
end

function makeGradeFinder(grades,nd)
    function (un)
        id = nd[un]["id"]
        row = findInDF(grades,:UNIV_ID,id)
        ismissing(row) ? missing : row[:GRADE]
    end
end

function makeMTFinder(mtdf,idd)
    function (id,col)
        un = idd[id]["username"]
        row = findInDF(mtdf,:username,un)
        ismissing(row) ? missing : row[col]
    end
end


function extractDFWframe(sem_grade_df,mt_df)
    addGrades(mt_df)
    iuids = getIUIDs()
    idd = makeIDDict(iuids)
    mtf = makeMTFinder(mt_df,idd)
    sem_grade_df.USERNAME = (x->idd[x]["username"]).(sem_grade_df.UNIV_ID)
    sem_grade_df.MT_RAW = (x->mtf(x,:score)).(sem_grade_df.UNIV_ID)
    sem_grade_df.MT_CURVE = (x->mtf(x,:curve)).(sem_grade_df.UNIV_ID)
    sem_grade_df.MT_GRADE = (x->mtf(x,:grade)).(sem_grade_df.UNIV_ID)
    dfw_df = select(sem_grade_df,
                    :TERM,
                    :COURSE,
                    :UNIV_ID,
                    :USERNAME,
                    :MT_RAW,
                    :MT_CURVE,
                    :MT_GRADE,
                    :GRADE=>:SEM_GRADE)
    term = dfw_df[1,:TERM]
    course = dfw_df[1,:COURSE]
    fn = "$(course)_$(term)_DFW.csv"
    CSV.write(fn,dfw_df;quotestrings=true)
    dfw_df
end
        
function invertPermutation(vec)
    len = length(vec)
    res = zeros(Int8,len)
    for i=1:len
        res[vec[i]] = i
    end
    res
end

function invertVerMap(vm)
    len = length(vm)
    res = []
    for i=1:len
        push!(res,invertPermutation(values(vm[i])))
    end
    res
end
