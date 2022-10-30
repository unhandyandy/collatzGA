
using DataFrames,Statistics#,Plots,StatPlots

include("../mydefs.jl")

#using Memoize


## vvvv from ptb on stackoverflow vvvv
dictreplace!(ex, s, v) = ex

dictreplace!(ex::Symbol, s, v) = s == ex ? v : ex

function dictreplace!(ex::Expr, s, v)
    for i=1:length(ex.args)
        ex.args[i] = dictreplace!(ex.args[i], s, v)
    end
    ex
end

dictreplace(ex::Expr, s, v) = dictreplace!(copy(ex), s, v)
dictreplace(ex::Symbol, s, v) = dictreplace!(ex, s, v)

function dictreplace_all(expr, kys, dsym)
    for k in kys
        expr = dictreplace(expr, k, :($(dsym)[$(QuoteNode(k))]))
    end
    expr
end
## ^^^^ from ptb on stackoverflow ^^^^

function get_atoms(atoms::Base.KeySet{Symbol},ex::Expr)
    args = ex.args
    ats = []
    for a in args
        if typeof(a)==Symbol && a in atoms
            push!(ats,a)
        elseif typeof(a)==Expr
            push!(ats,get_atoms(atoms,a)...)
        end
    end
    return ats       
end

function get_atoms(atoms::Base.KeySet{Symbol},sym::Symbol)
    return (sym in atoms) ? [sym] : []
end

dict = Dict{Symbol,Any}
ranks = dict()

function get_rank(rules,sym::Symbol)
    #println(sym)
    return get!(ranks,sym) do
        atoms = keys(rules)
        as = get_atoms(atoms,rules[sym])
        length(as)==0 ? 0 : 1 + max(map(a -> get_rank(rules,a),as)...)
    end
end
        
    
function partition_rules(rules)::Array{dict,1}
    maxrnk = max([get_rank(rules,s) for s in keys(rules)]...) + 1
    atoms = collect(keys(rules))
    atomsall = vcat(atoms,map(add_old,atoms))
    #println(atoms)
    parts = [ Dict() for _=1:maxrnk ]
    for a in atoms
        r = ranks[a] + 1
        #println(a)
        parts[r][a] = @eval function (args) $(dictreplace_all(rules[a], atomsall, :args)) end
    end
    return parts
end

function add_old(s::Symbol)
    return Meta.parse(string(s) * "Old")
end

function make_old_rules(d::dict)
    ks = keys(d)
    newd = Dict{Symbol,Any}()
    for k in ks
        newd[add_old(k)] = d[k]
    end
    return newd
end

# function update_values(curvals::dict,parts::Array{dict,1})
#     oldvals = make_old_rules(curvals)
#     for (a,v) in oldvals
#         ex = Expr(:(=),a,v)
#         #println(ex)
#         eval(ex)
#     end
#     newvals = dict()
#     for rules in parts
#         for (a,v) in rules
#             ex = Expr(:(=),a,v)
#             #println(ex)
#             eval(ex)
#             newvals[a] = eval(a)
#         end
#     end
#     return newvals
# end

function update_values(curvals::dict,parts::Array{dict,1})
    oldvals = make_old_rules(curvals)
    newvals = dict()
    allvals = merge(newvals,oldvals)
    for rules in parts
        for (a,v) in rules
            #println(v)
            #println(allvals)
            val = v(allvals)
            newvals[a] = allvals[a] = val
        end
    end
    return newvals
end    

function check_cond(dict,cond)
    @read_dict(dict)
    return eval(cond)
end

function run_till(state::Dict{Symbol,Any},parts,cond)
    st = copy(state)
    while !check_cond(st,cond)
        st = update_values(st,parts)
    end
    return st
end

# function make_data(init,parts,cond,n)
#     data = [ convert(DataFrame,run_till(init,parts,cond)) for _=1:n ]
#     return vcat(data...)
# end

