read("params.gp");
read("res_fields.gp");
read("bounds.gp");
read("primecoefs.gp");

mfreducible(f,flag = 0) = {
	my(reducible = List());
	my([N,k,eps,Pf,Pfcyclo] = mfparams(f));
	my(faN = factor(N));
	my([G,eps] = znchar(eps), eps_ord = charorder(G,eps));										

/* Compute the absolute extension Kf/Q */

	my([Pfabs,nthroot,a] = rnfequation(Pfcyclo,Pf,1));  \\Compute an absolute polynomial defining Kf.
	my(nthroot_abs = [nthroot,eps_ord]);			    \\nthroot_abs is the generator of Q(eps) and y-a*nthroot_abs is a root of Pf.
	my(nthroot_rel = [Mod(t,Pfcyclo),eps_ord]);			\\keep a relative root of unity
	
/* Compute the possible parameters for the small primes */

	my(small_pr = [p | p <- [2..max(k+1,N*eulerphi(N))], isprime(p) && (p <= k+1 || (N*eulerphi(N))%p == 0)]); \\Compute the set of small primes
	my(params = List(),l,lambda,L,index_pr = List());
	for(i = 1,#small_pr,
		l = small_pr[i]; lambda = factor_ideal(Pfabs,l); \\ Compute the prime ideals above l in Kf/Q
		
		if(lambda == "index",
			if(flag == 0, listput(reducible,[l,"index"]), listput(index_pr,l)),
			for(j = 1,#lambda, listput(params, [[0,lambda[j]],red_params(N,k,G,eps,l,nthroot_abs,[0,lambda[j]])]) )
		)
	);
	
/* Compute the bounds for the small primes */

	my(bound_map = Map(),r,B,Bmax,eps1,eps2,m1,m2);
	my(a2 = mfcoef(f,2));
	for(i = 1,#params,
		[pr,L] = params[i]; l = pr.l;
		for(j = 1,#L,
			[eps1,eps2,m1,m2] = L[j];
			[B,r,kdash] = get_Bred([N,faN],k,l,a2,eps1,eps2,m1,m2);
			Bmax = max(Bmax,B);
			mapput(bound_map,[pr,[eps1,eps2,m1,m2]],[B,r,kdash])
		)
	); 
	for(i = 1,#index_pr,
		Bmax = max(Bmax,get_Bmax([N,faN],k,index_pr[i],a2));
	);	Bmax = max(Bmax,get_Bred([N,faN],k,oo,a2,Mod(0,1),Mod(0,1),0,k-1)[1]);

/* Generate all the needed coefficients of f and compute the constant C */

	my(vf = mfprimecoefs(f,max(Bmax,N)));	
	my(C = -bernoulli(k,G,eps,nthroot_rel)/(2*k));
	for(i = 1,#faN[,1],
		af = mapget(vf,faN[i,1]);
		C = C*af*(af - faN[i,1]^(k-1)*chareval(G,eps,faN[i,1],nthroot_rel))
	);
	
/* Compute the big primes */

	my(big_params = red_params(N,k,G,eps,oo));	\\Compute the parameters for the big primes
	my(big_pr = List(),nb_pr = #params, nb_ideals = List());
	
	my(eps1,eps2,m1,m2,B,r,o,nthrootE,vE,L,l,lambda,pr,s);
	for(i = 1,#big_params,
		[eps1,eps2,m1,m2] = big_params[i];	\\Get the parameters
		[B,r] = get_Bred([N,faN],k,oo,a2,eps1,eps2,m1,m2);	\\Compute the bound
		
		o = lcm(eps1[3],eps2[3]); if(eps_ord%o == 0, 		\\Compute the root of unity for the coefficients of E
			nthrootE = nthroot_rel, 
			nthrootE = [Mod(t,polcyclo(o,t)),o]
		); vE = mfprimecoefs([k,eps1,eps2,nthrootE],B);	\\Compute the coefficients of E

		L = big_primes(N,vf,Pf,Pfcyclo,vE,nthrootE[1].mod,C*(eps1 == [0,1,1]),r);	\\Compute the list of the big primes

		for(j = 1,#L,
			l = L[j]; if(l > k+1 && (N*eulerphi(N))%l != 0,	\\Select only the primes that are big
				pr = setsearch(big_pr,l);	\\Search if l has already been encounter
				if(pr == 0,
				
					pr = setsearch(big_pr,l,1); listinsert(big_pr,l,pr);					\\If no, add it to big_pr,
					lambda = factor_ideal(Pfabs,l); s = nb_pr+sum(m = 1,pr-1,nb_ideals[m]);		\\compute the prime ideals in Kf above it,
					
					if(lambda == "index",
						if(flag == 0, listput(reducible,[l,"index"]),listput(index_pr,l));
						listinsert(nb_ideals,0,pr),
						
						for(m = 1,#lambda,
							listinsert(nb_ideals,#lambda,pr); mapput(bound_map,[[0,lambda[m]],[eps1,eps2,0,k-1]],[B,r,k]);
							listinsert(params,[[0,lambda[m]],List([[eps1,eps2,0,k-1]])],s+m)	\\and add the parameters to params.
						)
					),
					
					s = nb_pr+sum(m = 1,pr-1,nb_ideals[m]);									\\If yes,
					for(m = 1,nb_ideals[pr],
						listput(params[s+m][2],[eps1,eps2,0,k-1]);								\\just add the parameters in the right spot.
						mapput(bound_map,[params[s+m][1],[eps1,eps2,0,k-1]],[B,r,k])
					)
				)
			)
		)
	);
	
/* Deal with the primes dividing the index */
	
	if(flag == 1 && #index_pr > 0, Kf = nfinit([Pfabs,Vec(index_pr)]));
	for(i = 1,#index_pr,
		lambda = idealprimedec(Kf,index_pr[i]);
		for(j = 1,#lambda,
			L = nfmodprinit(Kf,lambda[j]);
			p = red_params(N,k,G,eps,index_pr[i],nthroot_abs,[Kf,L]); listput(params,[[1,L],p]);
			for(m = 1,#p,
				[eps1,eps2,m1,m2] = p[m];
				[B,r,kdash] = get_Bred([N,faN],k,index_pr[i],a2,eps1,eps2,m1,m2);
				mapput(bound_map,[[1,L],[eps1,eps2,m1,m2]],[B,r,kdash])
			)
		)
	);
	
/* Embbed the coefficient of f in Kf/Q */

	my(vf_abs = Map()); vf = Mat(vf);
	for(i = 1,#vf[,1],
		mapput(vf_abs,vf[i,1],substvec(liftall(vf[i,2]),[y,t],[y-a*nthroot_abs[1],nthroot_abs[1]]))
	); C = substvec(liftall(C),[y,t],[y-a*nthroot_abs[1],nthroot_abs[1]]);
	
/* Check congruences */

	my(eps1,eps2,m1,m2,lambda,p,B,r,kdash,big,L);
	for(i = 1,#params,
		[lambda,p] = params[i];
		for(j = 1,#p,
			[eps1,eps2,m1,m2] = p[j];
			[B,r,kdash] = mapget(bound_map,[lambda,p[j]]);

			o = lcm(eps1[3],eps2[3]); if(eps_ord%o == 0, 	\\Compute the root of unity for the coefficients of E
				nthrootE = nthroot_abs, 
				nthrootE = [Mod(t,polcyclo(o,t)),o]
			);
			vE = mfprimecoefs([kdash,eps1,eps2,nthrootE],B);	\\Compute the coefficients of E
			big = (lambda.l > k+1 && (N*eulerphi(N))%lambda.l != 0);
			if(lambda[1] == 1,L = [Kf,lambda[2]],L = lambda);			
			bool = check_cong(L,N,Pfabs,vf_abs,nthrootE[1].mod,vE,m1,r,big,C*(eps1 == [0,1,1]));
			if(bool,
				listput(reducible,[lambda[2],[eps1[1..2],eps2[1..2],m1,m2]])
			)
		)
	);

	return([[Pfabs,nthroot_abs[1],a],Vec(reducible)])
}