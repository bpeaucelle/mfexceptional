read("mfred.gp");

charevalp(G,eps,p,nthroot) = {
	my(epsp = znchardecompose(G,eps,p));
	my(epsp_dash = chardiv(G,eps,epsp));
	my([G0,eps0] = znchartoprimitive(G,epsp_dash));
	return(chareval(G0,eps0,p,nthroot))
}

chars_dih(G,eps,flag = 0, flag_G = 0) = {
	my(N = G.mod, ceps = znconreyconductor(G,eps));	
	my(fa = factor(N));
	my(fa0 = [fa[i,1] | i <- [1..#fa[,1]], fa[i,2] == 1 && ceps%fa[i,1] != 0]);
	
/* Generate all the quadratic characters */
	my(chars = [n | n <- [0..N-1], gcd(n,N) == 1 && znorder(Mod(n,N)) == 2]);
	
/* Check the condition at primes such that v_p(N) = 1 and v_p(ceps) = 0 */
	my(TNeps = List(),ctheta,theta0,Gtheta,bool,j);
	for(i = 1,#chars,
		ctheta = znconreyconductor(G,chars[i],&theta0);
		if(type(ctheta) == "t_VEC", ctheta = ctheta[1]);
		
		bool = 1; j = 1;
		while(bool && j <= #fa0,
			if(ctheta%fa0[j] == 0, bool = 0, j++)
		);
		
		if(bool,
			Gtheta = znstar(ctheta,1);
			if(flag_G = 0,
				listput(TNeps,[znconreyexp(Gtheta,theta0),ctheta]),
				Gtheta = znstar(ctheta,1);
				listput(TNeps,[znconreyexp(Gtheta,theta0),Gtheta])
			)
		)
	);
	
	if(flag == 0,
		return(Vec(TNeps)),
		return([Vec(TNeps),fa,fa0])
	)
}


dih_bound(l,Gtheta,theta,e,params,flag2 = 1) = {
	my([N,k,c] = params,fa,fa0);
	if(#params == 2,
		fa = factor(N);
		fa0 = [fa[i,1] | i <- [1..#fa[,1]], fa[i,2] == 1 && c%fa[i,1] != 0],
		[fa,fa0] = params[4..5]
	);

	my(ctheta = Gtheta.mod,v2N,v2theta,v2 = 0);
	if(flag2 == 1,
		v2N = valuation(N,2); v2theta = valuation(ctheta,2);
		v2 = vecmax([0,2*v2theta-v2N,v2theta+valuation(c,2)-v2N])
	);
	if(l == oo, return(floor((2^v2)* N * k * prod(i = 1,#fa[,1],1+1/fa[i,1])/12)));

	my(P = List());
	for(i = 1,#fa0,
		if(fa0[i] != l && (chareval(Gtheta,theta,fa0[i],[-1,2])*p^(e*(l-1)/2))%l == -1,
			listput(P,fa0[i])
		)
	);
	my(a = if(N%l == 0,3,l+1));
	
	my(B = floor((2^v2)* N * (k+a*(1+e*(l-1)/2)) * prod(i = 1,#P,P[i]) * prod(i = 1,#fa[,1],1+1/fa[i,1])/12));
	return([B,P])
}

dih_big(vf,Gtheta,theta,lvl,B) = {
	my(L = 0,j = 1,ctheta = Gtheta.mod);
	my([N,k,ceps,epsp] = lvl);
	
	while(L != 1 && vf[j,1] <= B,
		p = vf[j,1]; af = vf[j,2];
		if(N%p != 0,
			L = gcd(L,norm(af*(1-chareval(Gtheta,theta,p,[-1,2]))))
		); 
		if(N%p == 0 && (N/ceps)%p != 0,
			if(ctheta%p != 0,
				L = gcd(L,af*(1-chareval(Gtheta,theta,p,[-1,2]))),
				L = gcd(L,norm(af*af - p^(k-1)*mapget(epsp,p)*charevalp(Gtheta,theta,p,[-1,2])))
			)
		); j++
	)
}

dih_cong(lambda,vf,Gtheta,theta,e,lvl,B,P) = {
	my([N,k,fa0,fa1,epsp] = lvl);
	my(l = lambda.l,af);
	my(bool = 1, j = 1);
	while(bool && j <= #vf[,1] && vf[j,1] <= B,
		if(flag == 1 || vf[j,1] != l,
			p = vf[j,1];
			af = modpr(vf[j,2],lambda);
		
			if(N%p != 0,
				bool = (af*(1-p^(e*(l-1)/2)*chareval(Gtheta,theta,p,[-1,2])) == 0),
				
				if(setsearch(fa0,p) && flag == 0,
					bool = (setsearch(P,p) || p%l == -1)
				); if(setsearch(fa1,p),
					if(Gtheta.mod%p != 0,
						bool = (af*(1-chareval(Gtheta,theta,p,[-1,2])) == 0),
						bool = (af*af - p^(k-1+e*(l-1)/2)*mapget(epsp,p)*charevalp(Gtheta,theta,p,[-1,2]) == 0)
					)
				)
			)
		); j++
	); return(bool)
}

mfdihedral(f,flag = 0) = {
	if(mfisCM(f) != 0, return("The form f has CM"));
	
	my([N,k,eps0,Pfrel,Pfcyclo] = mfparams(f));
	
	my([G0,eps0] = znchar(eps0),G = znstar(N,1));
	my(eps = znconreyexp(G,zncharinduce(G0,eps0,N)));
	my(ceps = znconreyconductor(G,eps)); if(type(ceps) == "t_VEC", ceps = ceps[1]);
	my([Pfabs,nthroot,a] = rnfequation(Pfcyclo,Pfrel,1));
	nthroot = [nthroot,charorder(G0,eps0)];
	
	my(a2 = mfcoef(f,2),flag2);
	if(a2 == 0 && N%2 == 0, flag2 = 1,flag2 = 0);

/* Computation of the primes possible twist (that are unramified at l) and of the various prime factors of N */

	my([TNeps,fa,fa0] = chars_dih(G,eps,1,1));
	my(fa1 = [p | p <- fa[,1], (N/ceps)%p != 0],eps2,c2);
	if(setsearch(fa1,2), 
		c2 = znconreyconductor(G0,znchardecompose(G0,eps0,2),&eps2);
		if(type(c2) == "t_VEC", c2 = c2[1])
	);

/* Computation of the small primes */

	my(small_pr = List(),Np1 = prod(i = 1,#fa[,1],fa[i,1]+1));
	forprime(l = 3,max(N+1,2*k-1),
		if(l < k-1 || l == 2*k-1 || l == 2*k-3 || N%l == 0 || Np1%l == 0,
			listput(small_pr,l)
		)
	);

/* Computation of the various bounds */

	my(bound_map = Map(),params = List());
	my(Bmax = dih_bound(oo,znstar(8,1),3,0,[N,k,ceps,fa,fa0],flag2),B,P);
	my(theta,ctheta,Gtheta,thetaN,o,bool,n);
	for(i = 1,#small_pr,
		l = small_pr[i]; paramsl = List();
		
	/* Bound for (theta,e) = (1,1) */
	
		[B,P] = dih_bound(l,znstar(1,1),0,1,[N,k,ceps,fa,fa0]);
		mapput(bound_map,[l,0,znstar(1,1),1],[B,P]);
		listput(paramsl,[0,znstar(1,1),1]);
		Bmax = B;

	/* Select the characters of TNeps that are compatible with the prime l */
		
		for(j = 1,#TNeps,
			[theta,Gtheta] = TNeps[j]; ctheta = Gtheta.mod;
			if(#fa1 > 0, thetaN = zncharinduce(Gtheta,theta,G));
			bool = (ctheta%l != 0); n = 1;
			while(bool && n <= #fa1,
				if(ctheta%fa1[n] == 0,
					o = charorder(G,znchardecompose(G,charmul(G,eps,thetaN),fa1[n]));
					if(isprimepower(l*o) == 0, bool = 0)
				); n++
			);
	
	/* Compute the bound for e = 0 and e = 1, if the character is compatible */
	
			if(bool,
				listput(paramsl,[theta,Gtheta,0]);
				[B,P] = dih_bound(l,Gtheta,theta,0,[N,k,ceps,fa,fa0]);
				mapput(bound_map,[l,theta,Gtheta,0],[B,P]);
				Bmax = max(Bmax,B);
				
				listput(paramsl,[theta,Gtheta,1]);
				[B,P] = dih_bound(l,Gtheta,theta,1,[N,k,ceps,fa,fa0]);
				mapput(bound_map,[l,theta,Gtheta,1],[B,P]);
				Bmax = max(Bmax,B)
			)
		); listput(params,[l,paramsl]);
	);
	
/* Compute the big primes */

	/* Select the character that are compatible with the big primes */
	
	my(TNeps_big = List());
	for(i = 1,#TNeps,	
		[theta,Gtheta] = TNeps[i]; ctheta = Gtheta.mod;
		bool = 1; n = 1;
		while(bool && n <= #fa0,
			if(chareval(Gtheta,theta,fa0[n],[-1,2]) != 1, bool = 0);
			n++
		); n = 1; thetaN = zncharinduce(Gtheta,theta,G);
		while(bool && n <= #fa1,
			if(ctheta%fa1[n] == 0,
				if(charorder(G,znchardecompose(G,charmul(G,thetaN,eps),fa1[n])) != 1,
					bool = 0
				)
			); n++
		); if(bool, 
			listput(TNeps_big,[theta,Gtheta]);
			B = dih_bound(oo,Gtheta,theta,0,[N,k,ceps,fa,fa0],flag2);
			mapput(bound_map,[oo,theta,Gtheta],B);
			Bmax = max(Bmax,B)
		)
	);
	
	/* Compute the needed coefficients of f */

	my(vf = mfprimecoefs(f,Bmax));
	my(vf_abs = Map(), vf = Mat(vf));
	for(i = 1,#vf[,1],
		mapput(vf_abs,vf[i,1],substvec(liftall(vf[i,2]),[y,'t],[y-a*nthroot[1],nthroot[1]]));
	); vf_abs = Mat(vf_abs);
	my(epsp = Map(matrix(#fa1,2,i,j,if(j == 1,fa1[i],charevalp(G,eps,fa1[i],nthroot)))));
	
	/* Compute the big primes for each character */

	my(big = List(),s,lg = #params);
	for(i = 1,#TNeps_big,
		[theta,Gtheta] = TNeps_big[i];
		B = mapget(bound_map,[oo,theta,Gtheta]);
		L = factor(dih_big(vf_abs,Gtheta,theta,[N,k,ceps,epsp],B))[,1];
		for(j = 1,#L,
			l = L[j];
			if(l >= k-1 && l != 2*k-1 && l != 2*k-3 && N%l != 0 && Np1%l != 0,
				s = setsearch(big,l); if(s != 0,
					listput(params[s][2],[theta,Gtheta,0]);
					mapput(bound_map,[l,theta,Gtheta,0],[B,[]]),
					
					s = setsearch(big,l,1);
					listinsert(big,l,s);
					mapput(bound_map,[l,theta,Gtheta,0],[B,[]]);
					listinsert(params,[l,List([[theta,Gtheta,0]])],s+lg)
				)
			)
		)
	);
	my(dihedral = List(),index = List(),p,pr);
	for(i = 1,#params,
		[l,p] = params[i];
		pr = factor_ideal(Pfabs,l);
		if(pr == "index", listput(index,l),
			
			for(j = 1,#pr,
				lambda = pr[j];
				for(n = 1,#p,
					[theta,Gtheta,e] = p[n];
					[B,P] = mapget(bound_map,[l,theta,Gtheta,e]);
					if(dih_cong([0,lambda],vf_abs,Gtheta,theta,e,[N,k,fa0,fa1,epsp],B,P),
						listput(dihedral,[lambda,[theta,Gtheta.mod,e]])
					)
				)
			)
		)
	); if(flag == 0,return([dihedral,index]));
	
	my(red = mfreducible(f)[2],dih = List());
	for(i = 1,#dihedral,
		lambda = dihedral[i][1];
		bool = 1; j = 1;
		while(bool && j <= #red,
			if(red[j][1] == lambda, bool = 0);
			j++
		); if(bool, listput(dih,dihedral[i]))
	); return([dih,index])
}