

## params = current parameter vector
## vals = current values assigned to positions
## scores = scores of positions, one per row
## errors = vals - scores*params

## compute delta of params for one row of scores
function delForRow(scoreRow,err)
    m2 = norm(scoreRow)^2
    return err * scoreRow/m2
end
    
function delTotal(scores,errors)
    (r,c) = size(scores)
    del = zeros(c,1)
    for i = 1:r
        d = delForRow(scores[i,:],errors[i])
        del += d/r
    end
    return del
end

function updateCoeffs(scores,vals,oldcoeffs)
    (r,c) = size(scores)
    errors = vals - scores*oldcoeffs
    i = 0
    prms = oldcoeffs
    while (max(abs.(errors)...) > 0.001 && i <= 100)
        i += 1
        del = delTotal(scores,errors)
        prms += del
        errors = vals - scores*prms
    end
    return prms
end
