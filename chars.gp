read("res_fields.gp");

Tlifts(N,l) = return([m | m <- [0..N-1], gcd(m,N) == 1 && znorder(Mod(m,N))%l != 0]);
addhelp(Tlifts,"Compute the subgroup of Dirichlet characters modulo N that are Teichmuller lifts modulo l, i.e. the characters of prime-to-l order.");

Tlift(N,eps,l, Lifts = []) = {
	if(l == oo,return(eps));
    if(Lifts == [], Lifts = Tlifts(N,l));
	my(ord_quo); for(i = 1,#Lifts,
		ord_quo = znorder(Mod(eps,N)/Mod(Lifts[i],N));
		if(isprimepower(l*ord_quo), return(Lifts[i]))
    )
}
addhelp(Tlift,"Compute the Teichmuller lift modulo l of the Dirichlet character modulo N, eps given by its Conrey label. The list of Teichmuller lifts modulo l can be given as an optional argument.");

get_kl(G,epsl,lambda,nthroot) = {
    my(l = lambda.l);
    if(l == 2 || charorder(G,epsl) == 1, return(0)); \\ kl = 0 if l = 2 or epsl has prime-to-l order
	
/* Compute a generator of (Z/lZ)* in (Z/NZ)* */

    my(N = G.mod, gen = lift(znstar(l).gen[1]));
    while(gcd(N,gen) != 1, gen = gen+l);
	
/* Find kl such that epsl(gen) =: e = gen^kl mod(lambda) */
	
    my(e = modpr(chareval(G,epsl,gen,nthroot),lambda));
    my(n = gen);
    for(kl = 1,l-2,
        if(e == n,return(kl));
        n = lift(Mod(n*gen,l))
    )
}
addhelp(get_kl,"Let G=znstar(N,1) and epsl be a character modulo N of conductor a power of l. Compute the only integer kl between 0 and l-2, such that epsl is congruent to the kl-th power of the cyclotomic character modulo lambda.");

red_params(N,k,G,eps,l,nthroot,lambda) = {
    my(epsl);
    if(l == oo || (N == 2 && l == 2), epsl = 1, epsl = znchardecompose(G,eps,l)); \\Compute the l-part of eps
    my(eps0 = chardiv(G,eps,epsl));												  \\Compute the prime-to-l part of eps
	eps0 = znconreyexp(znstar(N,1),zncharinduce(G,eps0,N));						  \\and its Conrey label modulo N
    
/* Compute the pairs (m1,m2) */

	my(M); if(l == oo,M = [[0,k-1]],
		M = List();
		my(kl = get_kl(G,epsl,lambda,nthroot)); \\ Compute kl
		for(m1 = 0,l-2,
			for(m2 = m1,l-2,
				if((m1+m2-k-kl+1)%(l-1) == 0 && (2*N%l == 0 || m2+l*m1 <= k-1),listput(M,[m1,m2]))
			)
		)
	);

/* Compute the pairs (eps1,eps2) */

	my(E = List());
    my(Lifts = if(l == oo,[m | m <- [0..N-1], gcd(m,N) == 1],Tlifts(N,l)));	\\Compute the Teichmuller lifts mod l
    eps0 = Tlift(N,eps0,l,Lifts); 										    \\Compute the Teichmuller lift modulo l of the prime-to-l part of eps

    my(G1,G2,eps1,eps2,c1,c2,o1,o2,fa,bool,i);
    for(i1 = 1,#Lifts,
        for(i2 = i1,#Lifts,
            if(Mod(Lifts[i1]*Lifts[i2],N) == Mod(eps0,N),
			
                [G1,eps1] = znchartoprimitive(G,Lifts[i1]); c1 = G1.mod;		\\Compute the primitive characters
                [G2,eps2] = znchartoprimitive(G,Lifts[i2]); c2 = G2.mod;		\\associated to eps1 and eps2
                fa = factor(N/(c1*c2)); bool = (l == oo || c1*c2%l != 0); i = 1;\\Check the conditions on the conductors
				
                while(bool == 1 && i <= #fa[,1],
                    bool = (fa[i,1] == l || (fa[i,2] >= 0 && fa[i,2] <= 2)); i++
                ); if(bool,
					o1 = charorder(G1,eps1); o2 = charorder(G2,eps2);
					if(type(eps1) != "t_INT",eps1 = znconreyexp(G1,eps1));
					if(type(eps2) != "t_INT",eps2 = znconreyexp(G2,eps2));
					listput(E,[[eps1,c1,o1],[eps2,c2,o2]])
				)
            )
        )
    );
    
/* Return only the non-redondant quadruplets, i.e. (eps1,eps2,m1,m2) and (eps2,eps1,m1,m2) if m1 != m2, and (eps1,eps2,m1,m2) if m1 = m2 */

    my(R = List()); 
	for(i = 1,#M,
		[m1,m2] = M[i]; if(m1 != m2,
			for(j = 1,#E, 
				[eps1,eps2] = E[j]; listput(R,[eps1,eps2,m1,m2]);
				if(eps1 != eps2, listput(R,[eps2,eps1,m1,m2]))
			),
			for(j = 1,#E,
				[eps1,eps2] = E[j];
				if(eps1[2] > eps2[2],listput(R,[eps1,eps2,m1,m2]),listput(R,[eps2,eps1,m1,m2]))
			)
		)
	); return(Vec(R))
}
addhelp(params,"Let G = znstar(N,1) and eps be a Dirichlet character modulo N represented by its Conrey label. Compute the set R_(N,k,eps)(lambda). nthroot must be such that chareval(G,eps,n,nthroot) = nthroot[1]^chareval(G,eps,n).");

mfnorm(af,Pf,Pfcyclo,aE,PEcyclo) = {
	my(n = poliscyclo(PEcyclo)/poliscyclo(Pfcyclo));
	my(t = variable(Pfcyclo),y = variable(Pf),x = variable(PEcyclo));
	
	my(P = subst(liftall(Pf),t,x^n), af = subst(liftall(af),t,x^n));
	my(a); if(type(aE) == "t_VEC",
		a = af*(af-liftall(aE[1]))*(af-liftall(aE[2])),
		a = af-liftall(aE)
	); return(polresultant(PEcyclo,polresultant(P,a,y),x))
}

big_primes(N,vf,Pf,Pfcyclo,vE,PEcyclo,C = 0,r = 1) = {
	my(pr = Mat(vE)[,1],af,aE);
	my(L); if(r > 1,L = 0,L = mfnorm(C,Pf,Pfcyclo,0,PEcyclo));
	
	my(i = 1,p); while(L != 1 && i <= #pr,
		p = pr[i]; if(r%p != 0,
			if(N%p == 0,
				af = mapget(vf,p); aE = mapget(vE,p);
				L = gcd(L, mfnorm(af,Pf,Pfcyclo,aE,PEcyclo)),
				L = gcd(L, mfnorm(mapget(vf,p),Pf,Pfcyclo,vecsum(mapget(vE,p)),PEcyclo))
			)
		); i++
	); return(factor(L)[,1])
}