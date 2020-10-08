bernoulli(k,G,eps,nthroot = 0) = {
  my(cond = G.mod);
  if(nthroot == 0,my(o = znorder(Mod(eps,cond)));nthroot = [Mod(t',polcyclo(o,'t)),o]);
  my(ber = 0,P = bernpol(k,var));
  for(i = 0,cond-1,
    ber = ber + subst(P,var,i/cond) * chareval(G,eps,i,nthroot)
  );
  return(cond^(k-1)*ber)
}
addhelp(bernoulli,"Compute the k-th Bernoulli number attached to the primitive character eps. nthroot can be set to [z,o] with z an o-th root of unity with o a multiple of the order of eps.");

mfcsteisenstein(k,G1,eps1,G2,eps2,nthroot) = {
    if(k >= 2 && eps1 != 0, return(0));
    if(k == 1 && eps1 != 0 && eps2 != 0, return(0));
    if(k == 2 && eps1 == 0 && Mod(eps2,G2.mod) == Mod(1,G2.mod), return((G2.mod-1)/24));
    if(k >= 2, return(-bernoulli(k,G2,eps2,nthroot)/(2*k)));
    if(k == 1 && eps1 == 0, return(-bernoulli(k,G2,eps2,nthroot)/(2*k)));
    return(-bernoulli(k,G1,eps1,nthroot)/(2*k))
}
addhelp(mfcsteisenstein,"Compute the constant coefficient of the Eisenstein series E_{k}^{eps1,eps2}.");

mfprimecoefs(f,B, flag = 0) = {
    my(pr = [p | p <- [0..B], isprime(p) || (flag == 1 && p == 0)]);
    if(type(f) == "t_VEC" && #f == 4,
        my(vf = mfcoefs(f,B));
        return(Map(matrix(#pr,2,i,j,if(j == 1,pr[i],vf[pr[i]+1]))))
    );
    if(type(f) == "t_Vec" && #f == 6,
        my([k,G1,eps1,G2,eps2,nthroot] = f);
        my(val1 = Map(matrix(G1.mod,2,i,j,if(j == 1,i-1,chareval(G1,eps1,i-1,nthroot)))));
        my(val2 = Map(matrix(G2.mod,2,i,j,if(j == 1,i-1,chareval(G2,eps2,i-1,nthroot)))));
        return(Map(matrix(#pr,2,i,j,if(j == 1,pr[i],if(pr[i] == 0,mfcsteisenstein(k,G1,eps1,G2,eps2,nthroot),[getmap(val1,pr[i]),(pr[i]^k)*getmap(val2,pr[i])])))))
    )
}
addhelp(mfprimecoefs,"There is two behaviours to this function. If f is a PARI modular form, it output a Map that give the Fourier coefficients of prime index of f up to B. If f is an Eisenstein series E_{k}^{eps1,eps2} represented by [k,G1,eps1,G2,eps2,nthroot], it output a Map that give for all the prime number up to B, [eps1(p),p^{k-1}eps2(p)]. If flag = 1, also give the constant coefficient.");