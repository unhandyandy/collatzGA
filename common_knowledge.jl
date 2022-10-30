
univ0 = [(x,y,z) for x=1:9 for y=1:9 for z=1:9]
univ1 = filter(t -> t[1]!=t[2] && t[2]!=t[3] && t[1]!=t[3],univ0)
cond(t) = sum(t) == 2*max(t...)
univ2 = filter(t -> cond(t),univ1)

posses_one(u,i,j,x) = length(unique(map(t->t[j],filter(t -> t[i]==x,u))))
posses(u,i,j) = [ posses_one(u,i,j,x) for x=1:9 ]
filter_posses(u,i,j,p) =  filter( t -> posses_one(u,i,j,t[i])==p,u)
possA = filter_posses(univ2,1,2,8)
