
function pointPercent( del )
	 dist1 = Normal( del, 1.414 )
	 dist2 = Normal( -del, 1.414 )
	 whitescore = ccdf( dist1,  0.91 ) + 0.5 * ( cdf( dist1, 0.91 ) - cdf( dist1, -2.03 ) )
	 blackscore =  cdf( dist2, -2.03 ) + 0.5 * ( cdf( dist2, 0.91 ) - cdf( dist2, -2.03 ) )
	 ( whitescore + blackscore ) / 2
end  