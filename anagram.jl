f(s,k) = if k==0
    [[]]
else
    ps = []
    for c in s
        cur = map(x->[c x...],f(setdiff(s,[c]),k-1))
        #println(cur)
        push!(ps,cur...)
    end
    ps
end
