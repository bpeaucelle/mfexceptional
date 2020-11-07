read("chars.gp");
read("res_fields.gp");
read("bounds.gp");
read("primecoefs.gp");

mfreducible(f) = {
	my([N,k,eps,Pf,Pfcyclo] = mfparams(f));
	my(reducible = List());
	my(faN = factor(N));

/* Compute G = znstar(N,1) and the Conrey label modulo N of eps */

	my([G,eps] = znchar(eps));
	if(G.mod != N, eps = zncharinduce(G,eps,N); G = znstar(N,1));
	eps = znconreyexp(G,eps); my(eps_ord = znorder(Mod(eps,N)));
	my([G0,eps0] = znchartoprimitive(G,eps));

/* Compute the absolute extension Kf/Q */

	my([Pfabs,nthroot,n] = rnfequation(Pfcyclo,Pf,1)); \\Compute an absolute polynomial defining Kf. nthroot is the generator of Q(eps) and y-n*nthroot is a root of Pf.
	my(nthroot_abs = [nthroot,eps_ord]);
	my(nthroot_rel = [Mod(t,Pfcyclo),eps_ord]);
	
/* Compute the possible parameters for the small primes */

	my(params = List());
	my(small_pr = [p | p <- [2..max(k+1,N*eulerphi(N))], isprime(p) && (p <= k+1 || (N*eulerphi(N))%p == 0)]);
	my(l,pr);
	for(i = 1,#small_pr,
		l = small_pr[i]; pr = factor_ideal(Pfabs,l);
		if(pr == "dede",listput(reducible,[l,"dede"]),
			for(j = 1,#pr,
				listput(params,[pr[j],red_params(N,k,G,eps,l,nthroot_abs,pr[j])])
			)
		)
	);
	
/* Compute the bounds for the small primes */

	my(bounds_map = Map(),r,B,Bmax,eps1,eps2,m1,m2);
	my(a2 = mfcoef(f,2));
	for(i = 1,#params,
		[pr,L] = params[i]; l = pr.l;
		for(j = 1,#L,
			[eps1,eps2,m1,m2] = L[j];
			[B,r] = get_Bred([N,faN],k,l,a2,Mod(eps1[1],eps1[2]),Mod(eps2[1],eps2[2]),m1,m2);
			Bmax = max(Bmax,B);
			mapput(bounds_map,[pr,[eps1,eps2,m1,m2]],[B,r])
		)
	); Bmax = max(Bmax,get_Bred([N,faN],k,oo,a2,Mod(0,1),Mod(0,1),0,k-1)[1]);
	
/* Generate all the needed coefficients of f and compute C */

	my(vf = mfprimecoefs(f,max(Bmax,N)),n,vE,nthrootE);	
	my(C = -bernoulli(k,G0,eps0,nthroot_rel)/(2*k));
	for(i = 1,#fa[,1],
		af = mapget(vf,fa[i,1]);
		C = C*af*(af-fa[i,1]^(k-1)*chareval(G0,eps0,fa[i,1],nthroot_rel))
	);
	
/* Compute the big primes */
	
	my(big_params = red_params(N,k,G,eps,oo));	\\Compute the parameters for the big primes
	my(big_pr = List(),nb_pr = #params, nb_ideals = List(),L);
	
	for(i = 1,#big_params,
		[eps1,eps2,m1,m2] = big_params[i]; eps1 = Mod(eps1[1],eps1[2]); eps2 = Mod(eps2[1],eps2[2]);	\\Get the parameters
		[B,r] = get_Bred([N,faN],k,oo,a2,eps1,eps2,m1,m2);	\\Compute the bound
		
		n = lcm(znorder(eps1),znorder(eps2));
		if(n <= eps_ord, nthrootE = nthroot_rel, nthrootE = [Mod(t,polcyclo(n)),n]);	\\Compute the root of unity for the coefficients of E
		vE = mfprimecoefs([k,eps1,eps2,nthrootE],B);	\\Compute the coefficients of E
		
		L = big_primes(N,vf,Pf,Pfcyclo,vE,nthrootE[1].mod,C*(eps1 == Mod(0,1)),r);	\\Compute the big primes
		
		for(j = 1,#L,
			l = L[j]; if(l > k+1 && (N*eulerphi(N))%l != 0,	\\Select only the primes that are big
				pr = setsearch(big_pr,l);	\\Search if l has already been encounter
				if(pr == 0,
					pr = setsearch(big_pr,l,1); listinsert(big_pr,l,pr);				\\If no add it to big_pr,
					l = factor_ideal(Pfabs,l); r = nb_pr+sum(m = 1,pr-1,nb_ideals[m]);	\\compute the prime ideals in Kf above it,
					if(l == "dede",
						listput(reducible,[l,"dede"]); listinsert(nb_ideals,0,pr),
						
						for(m = 1,#l,
							listinsert(nb_ideals,#l,pr);
							listinsert(params,[l[m],List([[[lift(eps1),eps1.mod],[lift(eps2),eps2.mod],0,k-1]])],r+m)	\\and add the parameters to params.
						)
					),
					
					r = nb_pr+sum(m = 1,pr-1,nb_ideals[m]);											\\If yes,
					for(m = 1,nb_ideals[pr],
						listput(params[r+m][2],[[lift(eps1),eps1.mod],[lift(eps2),eps2.mod],0,k-1])	\\just add the parameters in the right spot.
					)
				)
			)
		)
	);
	
	my(lambda,p);
	for(i = 1,#params,
		[lambda,p] = params[i];
		for(j = 1,#p,
			[eps1,eps2,m1,m2] = p[j];
			[B,r] = mapget(bound_map,[lambda,p[j]]);
			
		)
	);
	return(reducible)
}