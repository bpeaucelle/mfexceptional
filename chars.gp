Tlifts(N,l) = return([m | m <- [0..N-1], gcd(m,N) == 1 && znorder(Mod(m,N))%l != 0]);
addhelp(Tlifts,"Compute the subgroup of Dirichlet characters modulo N that are Teichmuller lifts modulo l, i.e. the characters of prime-to-l order.");

Tlift(N,eps,l, Lifts = []) = {
    if(Lifts == [], Lifts = Tlifts(N,l));
    for(i = 1,#Lifts,
        if(isprimepower(znorder(Mod(eps/Lifts[i],N))*l),
            return(Lifts[i])
        )
    )
}
addhelp(Tlift,"Compute the Teichmuller lift modulo l of the Dirichlet character modulo N, eps given by its Conrey label. The list of Teichmuller lifts modulo l can be given as an optional argument.");

get_kl(G,epsl,ideal,nthroot) = {
    my(l = ideal.p);
    if(l == 2 || charorder(G,epsl) == 1, return(0));
    my(N = G.mod, gen = lift(znstar(l).gen[1]));
    while(gcd(N,gen) != 0, gen = gen+l);
    my(e = coefmod(chareval(G,epsl,gen,nthroot),ideal));
    my(n = gen);
    for(kl = 1,l-2,
        if(e == n,return(kl));
        n = lift(Mod(n*gen,l))
    )
}
addhelp(get_kl,"Let G=znstar(N,1) and epsl be a character modulo N of conductor a power of l. Compute the only integer kl between 0 and l-2, such that epsl is congruent to the kl-th power of the cyclotomic character modulo lambda.");

params(N,k,G,eps,l,nthroot,ideal) = {
    my(epsl);
    if(N == 2 && l == 2, epsl = 1, epsl = znconreyexp(G,znchardecompose(G,eps,l)));
    my(eps0 = lift(Mod(eps/epsl,N)));
    
    /* Compute the pairs (m1,m2) */
    my(kl = get_kl(G,epsl,ideal,nthroot), M = List());
    for(m1 = 0,l-2,
        for(m2 = m1,l-2,
            if((m1+m2-k-kl+1)%(l-1) == 0,listput(M,[m1,m2]))
        )
    );
    
    /* Compute the pairs (eps1,eps2) */
    my(Lifts = Tlifts(N,l),E = List());
    eps0 = Tlift(N,eps0,l,Lifts);
    my(G1,G2,eps1,eps2,c1,c2,fa,bool,i);
    for(i1 = 1,#Lifts,
        for(i2 = i1,#Lifts,
            if(Mod(Lifts[i1]*Lifts[i2],N) == Mod(eps0,N),
                [G1,eps1] = znchartoprimitive(G,Lifts[i1]); c1 = G1.mod;
                [G2,eps2] = znchartoprimitive(G,Lifts[i2]); c2 = G2.mod;
                fa = factor(N/(c1*c2)); bool = (c1*c2%l != 0); i = 1;
                while(bool == 1 && i <= #fa[,1],
                    bool = (fa[i,1] == l || (fa[i,2] >= 0 && fa[i,2] <= 2));
                    i++
                ); if(bool, listput(E,[[znconreyexp(G1,eps1),c1],[znconreyexp(G2,eps2),c2]]))
            )
        )
    );
    
    my(R = List()); 
    for(i = 1,#M,
        for(j = 1,#E,
            if(M[i][1] != M[i][2],
                listput(R,[E[j][1],E[j][2],M[i][1],M[i][2]]);
                listput(R,[E[j][2],E[j][1],M[i][1],M[i][2]]),
                if(E[j][1][2] > E[j][2][2],
                    listput(R,[E[j][1],E[j][2],M[i][1],M[i][2]]),
                    listput(R,[E[j][2],E[j][1],M[i][1],M[i][2]])
                )
            )
        )
    ); return(Vec(R))
}