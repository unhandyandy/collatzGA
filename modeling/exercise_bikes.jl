
include("modeling.jl")

using Random,StatsFuns


function arrival_rate(t)
    return (t <=30) ? 1 : (t <= 60) ? 1.5 : 2/3 
end

function rand_arrivals(t)
    return round(arrival_rate(t) * randexp())
end

usertime = 8
waitline = 3

## num = # bikes
function make_state(num::Int64)
    dict(:tm => 0,
         :waiting => Int64[],
         :waitingtimes => Int64[],
         :arr => 0,
         :use => 0,
         :disapp => 0,
         :free => num,
         :bikes => zeros(Int64,num))
end

function update_bikes(bks::Array{Int64,1},numwts::Int64,nws::Int64)::Array{Int64,1}
    bksnew = map(b -> max(b - 1,0),bks)
    free = count(b -> b==0,bksnew)
    if numwts + nws >= free
        return sort(map(b -> b==0 ? usertime : b,bksnew))
    else
        diff = free - numwts - nws
        return vcat(zeros(Int(diff)),
                    filter(b -> b>0,bksnew),
                    usertime * ones(Int64(numwts + nws)))
    end
end

function update_waiting(exuse::Int64,wts::Array{Int64,1},arrs::Int64)::Array{Int64,1}
    numwts = length(wts)
    numwtsnew = min(numwts + arrs - exuse, waitline)
    newwts = vcat(wts,zeros(Int64(arrs)))
    newwts = newwts[exuse + 1:end]
    newwts .+= 1
    return newwts[1:Int64(numwtsnew)]
end

function new_wtimes(old,new)
    newwts = count(w -> w>1,new)
    return old[1:(length(old) - newwts)]
end

function make_updater(num)
    dict(:tm => :(tmOld + 1),
         :arr => :(Int64(rand_arrivals(tm))),
         :bikes => :(update_bikes(map(Int64,bikesOld),length(waitingOld),arr)),
         :free => :(count(b -> b==0,bikes)),
         :use => :(useOld + count(b -> b==8,bikes)),
         :waiting => :(update_waiting(use - useOld,waitingOld,arr)),
         :disapp => :(disappOld + max(length(waitingOld) + arr - (use - useOld + 3),0)),
         :waitingtimes => :(vcat(waitingtimesOld,new_wtimes(waitingOld,waiting))))
end

init_state = make_state(5)

rules = make_updater(5)

parts = partition_rules(rules)

function parse_data(d::dict)
    newd = dict()
    wts = d[:waitingtimes]
    nw = length(wts)
    newd[:numwaits] = nw
    newd[:avewait] = sum(wts)/nw
    newd[:use] = d[:use]
    newd[:disapp] = d[:disapp]
    return convert(DataFrame,newd)
end

## 0.5 * prob client is disappointed at least 4 times in month
function prob_cancel(d::Float64)
    return 0.5 * (1 - binomcdf(10,d,3))
end

function calc_aves(n,frame)
    use = df_get_entry(frame,:use,:mean)
    disnum = df_get_entry(frame,:disapp,:mean)
    disprb = disnum/(disnum + use)
    prcncl = prob_cancel(disprb)
    cancst = prcncl * 270 * 30
    bkscst = 100 * (n - 5)
    cstttl = bkscst + cancst
    newdict =  Dict("Bikes" => n,
                    "Prob Dspt" => disprb,
                    "Prob Cncl" => prcncl,
                    "Cncl Cost" => cancst,
                    "Bike Cost" => bkscst,
                    "Total Cost" => cstttl)
    return convert(DataFrame,newdict)
end

#data = [ run_till(init_state,parts,:(tm==90)) for _=1:100 ]
partsN = [ partition_rules(make_updater(n)) for n=5:12 ]

function make_data_n(num)
    return [ run_till(make_state(n),partsN[n - 4],:(tm==90)) for n=5:12,_=1:num ]
end

# pd = [ parse_data(dataN[n,j]) for n=1:8, j=1:100 ]
# vpd = map(r -> vcat(r...),[ pd[i,:] for i=1:8 ])
# aves = [ calc_aves(n+4,describe(vpd[n])) for n=1:8 ]
# vaves = vcat(aves...)


;

# │ Row │ Bike Cost │ Bikes │ Cncl Cost │ Prob Cncl   │ Prob Dspt │ Total Cost │
# ├─────┼───────────┼───────┼───────────┼─────────────┼───────────┼────────────┤
# │ 1   │ 0         │ 5     │ 2436.02   │ 0.300744    │ 0.393579  │ 2436.02    │
# │ 2   │ 100       │ 6     │ 1477.28   │ 0.18238     │ 0.305366  │ 1577.28    │
# │ 3   │ 200       │ 7     │ 805.28    │ 0.0994173   │ 0.238918  │ 1005.28    │
# │ 4   │ 300       │ 8     │ 367.165   │ 0.045329    │ 0.181545  │ 667.165    │
# │ 5   │ 400       │ 9     │ 151.539   │ 0.0187085   │ 0.137212  │ 551.539    │
# │ 6   │ 500       │ 10    │ 57.9472   │ 0.00715398  │ 0.103263  │ 557.947    │
# │ 7   │ 600       │ 11    │ 19.7619   │ 0.00243973  │ 0.0762517 │ 619.762    │
# │ 8   │ 700       │ 12    │ 5.69878   │ 0.000703553 │ 0.0543694 │ 705.699    │
