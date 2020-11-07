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
	if(type(pr) != "t_VEC" || #pr  != 3, return("Not a prime ideal."));
	return(pr[1])
}
pr.poly = {
	if(type(pr) != "t_VEC" || #pr  != 3, return("Not a prime ideal."));
	return(pr[2])
}
pr.gen_mod = {
	if(type(pr) != "t_VEC" || #pr  != 3, return("Not a prime ideal."));
	return(pr[3])
}

modpr(x,pr) = {
	return(subst( lift(Mod(liftall(x),pr.poly)) ,variable(pr.poly),pr.gen_mod))
}
addhelp(modpr,"Given an element of a number field as a PolMod and a prime ideal output by factor_ideal, reduce x modulo pr. Can be apply to vectors as well.")

