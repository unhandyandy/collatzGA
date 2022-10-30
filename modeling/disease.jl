
include("modeling.jl")

using Random,StatsFuns

randbinom(n,p) = binominvcdf(n,p,rand())

infection_prob = 0.2
rand_latent() = binominvcdf(5,0.5,rand())
rand_infect() = binominvcdf(5,0.5,rand())
rand_symp() = binominvcdf(10,0.3,rand())
rand_contacts() = binominvcdf(5,0.2,rand())

function make_state(pop,sick,qrate)
    days = rand(1:5,sick)
    lats = [ rand_latent() for _=1:sick ]
    # println(days)
    # println(lats)
    lats -= days 
    inf = filter(x -> x<=0,lats)
    lats = filter(x -> x>0,lats)
    inf = [ i + rand_infect() for i in inf ]
    symp = filter(x -> x<=0,inf)
    inf = filter(x -> x>0,inf)
    symp = [ s + rand_symp() for s in symp ]
    imm = count(x -> x<=0,symp)
    symp = filter(x -> x>0,symp)
    return dict(:tm => 0,
                :healthy => pop - sick,
                :latent => lats,
                :infectious => inf,
                :symptomatic => symp,
                :quarantine => [],
                :quarrate => qrate,
                :immune => imm,
                :newlat => 0,
                :newinf => 0,
                :newsymp => 0,
                :newquar => 0)
end

function update_healthy(h,i)
    contacts = [ rand_contacts() for _=1:h ]
    inf = map(n -> randbinom(n,i/(i+h)),contacts)
    probs = map(m -> 1 - (1 - infection_prob)^m,inf)
    return count(x -> x===false,map(randBool,probs))
end

function update_latent(newsick,lat)
    newlat = [ rand_latent() for _=1:newsick ]
    newlat = vcat(lat .- 1,newlat)
    return filter(x -> x>0,newlat)
end

function update_infectious(moves,old)
    new = [ rand_infect() for _=1:moves ]
    new = vcat(old .- 1,new)
    return filter(x -> x>0, new)
end
    
function update_symp(moves,old,qrate)
    n = randbinom(moves,1 - qrate)
    new = [ rand_symp() for _=1:n ]
    new = vcat(old .- 1,new)
    return filter(x -> x>0,new)
end

function update_quar(moves,old)
    new = [ rand_symp() for _=1:moves ]
    new = vcat(old .- 1,new)
    return filter(x -> x>0,new)
end

rules = dict(:tm => :(tmOld + 1),
             :healthy => :(update_healthy(healthyOld,
                                          length(infectiousOld) + length(symptomaticOld))),
             :newlat => :(healthyOld - healthy),
             :latent => :(update_latent(newlat,latentOld)),
             :newinf => :(newlat + length(latentOld)-length(latent)),
             :infectious => :(update_infectious(newinf,infectiousOld)),
             :newsymp => :(newinf + length(infectiousOld) - length(infectious)),
             :symptomatic => :(update_symp(newsymp,symptomaticOld,quarrateOld)),
             :newquar => :(newsymp + length(symptomaticOld) - length(symptomatic)),
             :quarantine => :(update_quar(newquar,quarantineOld)),
             :immune => :(immuneOld + newquar + length(quarantineOld) - length(quarantine)),
             :quarrate => :quarrateOld)


parts = partition_rules(rules)

function extract_info(d,i,q)
    nd = dict(:infected => d[:immune],
              :tm => d[:tm],
              :uninfected => d[:healthy],
              :init_infect => i,
              :qrate => q)
    return convert(DataFrame,nd)
    #return nd
end

tab(i,q) = data[(data[:qrate].==q) .& (data[:init_infect].==i),:]

fix_q(q) = [ describe(tab(i,q))[1:5] for i=5:5:20 ]
fix_i(i) = [ describe(tab(i,q))[1:5] for q=0:0.1:1.0 ]

#data = vcat([ extract_info(run_till(make_state(1000,i,q),parts,:(healthy + immune == 1000)),i,q) for i=5:5:20 for q=0:0.1:1.0 for _=1:1000 ]...)

#[ describe(tab(i,q)) for i=5:5:20, q=0:0.1:1.0 ]

;
