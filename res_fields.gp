dede(Q,l,fa = 0) = {
  if(fa == 0, fa = factormod(Q,l));
  my(G = 1);
  for(i = 1,#fa[,1], if(fa[i,2] >= 2,G *= fa[i,1]));
  my(F = Mod((factorback(lift(fa))-Q)/l,l));
  return(poldegree(gcd(F,G)) == 0);
}
addhelp(dede,"Apply the Dedekind criteria to determine if the prime number l divides the index of Q. The factorization of Q modulo l can be given as a third argument.")

factor_ideal(P,l,fa = 0) = {
	if(fa == 0, fa = factormod(P,l));
	if(dede(P,l,fa) == 0, return("dede"));
	
	my(pr = List(),pol,gen);
	for(i = 1,#fa[,1],
		pol = lift(fa[i,1]); gen = ffgen(fa[i,1]);
		listput(pr,[l,pol,gen])
	); return(Vec(pr))
}
addhelp(factor_ideal,"Compute the ideal decomposition of l in the number field defined by P in the case where l doesn't divide the index of P. The factorization of P modulo l can be given as a third argument.")

pr.l = {
	if(pr[1] == 0,
		return(pr[2][1]),
		return(pr[2].p)
	)
}

pr.poly = {
	if(pr[1] == 0, 
		return(pr[2][2]),
		
	)
}
pr.gen_mod = {
	if(pr[1] == 0,
		return(pr[2][3]),
		return(ffgen(Mod(nfmodpr(pr[1],0,pr[2]).mod,pr[2].p)))
	)
}

modpr(x,lambda) = {
	if(lambda[1] == 0,
		return(subst( lift(Mod(liftall(x),lambda.poly)) ,variable(lambda.poly),lambda.gen_mod)),
		if(type(x) == "t_VEC",return([nfmodpr(lambda[1],y,lambda[2]) | y <- x]),nfmodpr(lambda[1],x,lambda[2]))
	)
}
addhelp(modpr,"Given an element of a number field as a PolMod and a prime ideal output by factor_ideal, reduce x modulo pr. Can be apply to vectors as well.")

cong(l,N,pr,mapf,lambdaf,vf,mapE,lambdaE,vE,m1 = 0,r = 1,flag = 0,C = 0) = {
	my(bool = (flag == 0 || modpr(C,lambdaf) == 0), j = 1);

	while(bool && j <= #pr,
		if(r%pr[j] != 0 && (flag == 1 || pr[j] != l),
			af = ffmap(mapf,modpr(mapget(vf,pr[j]),lambdaf));
			aE = pr[j]^m1 * ffmap(mapE,modpr(mapget(vE,pr[j]),lambdaE));

			if(N%pr[j] == 0,
				a = af*(af-aE[1])*(af-aE[2]),
				a = af-aE[1]-aE[2]
			); bool = (a==0)
		); j++
	); return(bool)
}

check_cong(lambda,N,Pf,vf,PEcyclo,vE,m1,r,big = 0,C = 0) = {
	if(matsize(Mat(vE)) == [0,0],return(1));
	my(pr = Mat(vE)[,1],l = lambda.l,genf = lambda.gen_mod);

	if(Pf == PEcyclo,
		return(cong(l,N,pr,[genf,genf],lambda,vf,[genf,genf],lambda,vE,m1,r,big,C))
	);

	my(lambdaE = factor_ideal(PEcyclo,l));

	my(fa,g,mapf,mapE,pr,j,bool,af,aE,a);
	for(i = 1,#lambdaE,
		fa = factormod(lambdaE[i].poly,gen);
		[g,mapf] = ffextend(gen,fa[1,1]); mapE = ffembed(lambdaE[i].gen_mod,g);
	
		if(cong(l,pr,mapf,lambda,vf,mapE,lambdaE[i],vE,r,m1,big,C),return(1))
	); return(0)
}