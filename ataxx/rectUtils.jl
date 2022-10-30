

function neighbor1Locs(r,c,s)
    res = zeros(Int16,2,0)
    for i=r-1:r+1
        if 1<=i<=s
            for j=c-1:c+1
                if 1<=j<=s
                    if i!=r || j!=c
                        res = [res [i; j]]
                    end
                end
            end
        end
    end
    return res
end

# function neighbor1LocsB(r,c,s)
#     res = mapslices(x -> x + [r;c],relnbrs1,dims=1]
#     filter(res,max(res(:,i)

function neighbor2Locs(r,c,s)
    res = zeros(Int16,2,0)
    for k in [-2,2]
        for l in [-2,2]
            d1 = -sign(k)
            d2 = -sign(l)
            d0 = d1*d2>0
            for m = 0:3
                if d0
                    del = m*d1
                    rn = r+k+del
                    cn = c+l
                    if 1<=cn<=s && 1<=rn<=s
                        res = [res [rn;cn]]
                    end
                else
                    del = m*d2
                    rn = r+k
                    cn = c+l+del
                    if 1<=cn<=s && 1<=rn<=s
                        res = [res [rn;cn]]
                    end
                end
            end
        end
    end
    return res
end

"""
Returns locations at distance n from given one
r = row #
c = col #
s = size of board
n = distance = max(|r1-r2|,|c1-c2|)
"""
function neighborNLocs(r,c,s,n)
    res = zeros(Int16,2,0)
    # start at one of the four corners of the nxn square around (r,c)
    for k in [-n,n]
        for l in [-n,n]
            d1 = -sign(k)
            d2 = -sign(l)
            d0 = d1*d2>0
            # proceed CCW from corner along edge
            for m = 0:2*n-1
                if d0
                    del = m*d1
                    rn = r+k+del
                    cn = c+l
                    if 1<=cn<=s && 1<=rn<=s
                        res = [res [rn;cn]]
                    end
                else
                    del = m*d2
                    rn = r+k
                    cn = c+l+del
                    if 1<=cn<=s && 1<=rn<=s
                        res = [res [rn;cn]]
                    end
                end
            end
        end
    end
    return res
end


## Slower!
# function neighborNLocs(r,c,s,n)
#     s1 = r-n>0
#     s2 = c+n<=s
#     s3 = r+n<=s
#     s4 = c-n>0
#     x1 = max(1,c-n)
#     x2 = min(s,c+n)
#     y1 = max(1,r-n)
#     y2 = min(s,r+n)
#     res = zeros(Int16,2,0)
#     w = x2-x1+1
#     h = y2-y1+1
#     if s1
#         if s4
#             nv = Int16[y1*ones(w-1)';
#                        collect(x1+1:x2)';]
#             res = [res nv]
#         else
#             nv = Int16[y1*ones(w)';
#                        collect(x1:x2)';]
#             res = [res nv]
#         end
#     end
#     if s2
#         if s1
#             nv = Int16[collect(y1+1:y2)';
#                        x2*ones(h-1)']
#             res = [res nv]
#         else
#             nv = Int16[collect(y1:y2)';
#                        x2*ones(h)']
#             res = [res nv]
#         end
#     end
#     if s3
#         if s2
#             nv = Int16[y2*ones(w-1)';
#                        collect(x1:x2-1)';]
#             res = [res nv]
#         else
#             nv = Int16[y2*ones(w)';
#                        collect(x1:x2)';]
#             res = [res nv]
#         end
#     end
#     if s4
#         if s3
#             nv = Int16[collect(y1:y2-1)';
#                        x1*ones(h-1)']
#             res = [res nv]
#         else
#             nv = Int16[collect(y1:y2)';
#                        x1*ones(h)']
#             res = [res nv]
#         end
#     end
#     return res
# end
@test
"""
Check whether loc is a legal square
loc = [row;col]
"""
function onBoardQ(loc,numRows,numCols)
    return 0<loc[1]<=numRows && 0<loc[2]<=numCols
end
function  onBoardQ(loc,size)
    return  onBoardQ(loc,size,size)
end

"""
List of all legal moves in direction dir from given loc
"""
# function oneLine(pos,loc,dir,numRows,numCols,captureQ=false)
#     res = [
