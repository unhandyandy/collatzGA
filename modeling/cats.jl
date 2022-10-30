
rand_letter() = rand_element(["C","A","T","S"],[3,1,3,3])

init_state = Dict{Symbol,Any}(:c => false, :a => false, :t => false, :s => false,
                              :tm => 0,:letter => "")

rules = Dict{Symbol,Expr}(:tm => :(tmOld + 1),
                          :letter => :(rand_letter()),
                          :c => :(cOld || letter == "C"),
                          :a => :(aOld || letter == "A"),
                          :t => :(tOld || letter == "T"),
                          :s => :(sOld || letter == "S"))

parts = partition_rules(rules)

data = [ convert(DataFrame,
                 run_till(init_state,parts,:(c && a && t && s)))
         for _=1:1000 ]

frame = vcat(data...)

mean(frame[:tm])
