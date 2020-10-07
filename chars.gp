Tlifts(N,l) = {
    if(l == 0 || eulerphi(N)%l != 0,
        return([m | m <- [0..N-1], gcd(m,N) == 1])
    ); my(T = List());
    for(m = 0, N-1,
        if(gcd(m,N) == 1 && znorder(Mod(m,N))%l != 0,
            listput(T,m)
        )
    ); return(T)
}

Tlift(N,eps,l, Lifts = []) = {
    if(l == 0 || znorder(Mod(eps,N))%l != 0, return(eps));
    if(Lifts == [], Lifts = Tlifts(N,l));
    for(i = 1,#Lifts,
        if(isprimepower(znorder(Mod(eps/Lifts[i],N))*l),
            return(Lifts[i])
        )
    )
}

get_chars(G,eps,l = 0) = {
    my(chars = List(),eps1,eps10,eps2,eps20);
    my(N = G.mod, M = if(l == 0,N,N/l^valuation(N,l)));
	my(Lifts = Tlifts(N,l), Teps = Tlift(N,eps,l,Lifts));

	for(i1 = 1,#Lifts,
		eps1 = Lifts[i1]; for(i2 = 1,#Lifts,
			eps2 = Lifts[i2];
			if(Mod(eps1*eps2,N) == Mod(Teps,N),
				listput(chars,[eps1,eps2])
			)
		)
    );
	my(C = List(),c1,c2,G1,G2);
	for(i = 1,#chars,
		[eps1,eps2] = chars[i];
		c1 = znconreyconductor(G,eps1,&eps10);
		c2 = znconreyconductor(G,eps2,&eps20);
		if(M%(if(#c1 == 1,c1,c1[1])*if(#c2 == 1,c2,c2[1])) == 0,
			if(#c1 == 1,G1 = G, G1 = znstar(c1,1); c1 = c1[1]);
			if(#c2 == 1,G2 = G, G2 = znstar(c2,1); c2 = c2[1]);
			listput(C,[[znconreyexp(G1,eps10),c1],[znconreyexp(G2,eps20),c2]])
		)
	); return(C)
}