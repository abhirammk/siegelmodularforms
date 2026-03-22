q = 'q;
y = 'y;
p = 'p;
\\' 

L(s,D) = lfun(lfuncreate(D),s);


h(r,N) = my(D = (-1)^r*N); if(D%4 == 0 || D%4 == 1, (-1)^(ceil(r/2))*(r-1)!*N^(r-1/2)*2^(1-r)*Pi^(-r)*L(D,r), 0);

sqfac(n) = my(l = divisors(n), t= listcreate(), len = length(l)); for(i=1,len, if(n%l[i]^2==0, listput(t,l[i]))); return(Vec(t));

T(r,f,X)= sumdiv(f,d,moebius(d)*kronecker(X,d)*d^(r-1)*sigma(f/d,2*r-1));


 H(r,N) = my(D=(-1)^r*N, f = coredisc(D,1)[2], X = coredisc(D,1)[1],t); if(D==0, t=zeta(1-2*r), if(N > 0, if(D%4 == 0 || D%4 ==1, t= L(1-r,X)*T(r,f,X),t=0), t=0)); return(t);

div(x) = if(x==0,[1],divisors(x));

comdiv(V) = my(l = length(V), k1 = div(V[1])); for(i=2,l, k1 = setintersect(k1, div(V[i]))); return(Vec(k1));

nonzero(V) = my(k=listcreate(), len = length(V),t); if(vecprod(V)!=0, t= V, for(i=1,len, if(V[i]!= 0, listput(k,V[i]));); t=Vec(k); ); return(t);

posdisc(n) = my(k= listcreate()); for(a=0,n, for(b=0,n, for(c=0,n, if(4*a*c -b^2 >= 0, listput(k,[a,b,c]))))); return(Vec(k));

\\ C(a,b,c,k) = my(abc=comdiv(nonzero([a,b,c])),r=k-1,\\ D=4*a*c-b^2,t);t=vecsum(apply(x->x^(k-1)*H(r,D/x^2),abc));\\ return(bestappr(-2*k*t/(bernfrac(k)*zeta(3-2*k)), 50000));

C(a,b,c,k) =
{
  if(a==0 && b==0 && c==0, return(1));   \\ constant term of Eisenstein series

  my(abc=comdiv(nonzero([a,b,c])),r=k-1,D=4*a*c-b^2,t);
  t=vecsum(apply(x->x^(k-1)*H(r,D/x^2),abc));
  return(bestappr(-2*k*t/(bernfrac(k)*zeta(3-2*k)), 50000));
}

SiegelEisensteinCoeffs(k,ord)= my(K=posdisc(ord),len=length(K));for(i=2,len,x=K[i][1];y=K[i][2];z=K[i][3];print(K[i],": ",C(x,y,z,k)));

kill(y);

qpyser(k,ord) = my(K=posdisc(ord),len=length(K), t = 0); for(i=2,len,n=K[i][1]; l=K[i][2]; m = K[i][3]; t = t+ C(n,l,m,k)*q^(n)*y^(l)*p^(m)); return(t);


c = -43867/(2^10*3^5*5^2*7*53);


Delta(n,l,m)=4*n*m-l^2;


divW(n,l,m)=
{
  my( g=gcd(gcd(abs(n), abs(l)),abs(m)));
  if(g== 0, return([1]));
  divisors(g)
};

Ar(n,l,m,r)=
{
  if(n==0 && l==0 && m==0 ,  return(0));

  my(D=Delta(n,l,m));
  my(divs=divW(n,l,m));
  my(s=0);

  for(i = 1,length(divs),
    d=divs[i];
    if(D%(d^2)==0,
      s+=d^(r-1)*H(r,D/d^2);
    );
  );

  s
};

Theta(n,l,m)=
{
  if(n==0, return(0));   \\ cusp condition

  my(S=Ar(n,l,m,9));

  for(n1 =0,n,
   for(m1=0,m,
    for(l1 =-20,20,

      if(n1==0 && m1==0 && l1==0, next);
      if(n1==n && m1==m && l1==l, next);

      n2 =n-n1;
      m2= m-m1;
      l2=l-l1;

      S+= Ar(n1,l1,m1,3)*Ar(n2,l2,m2,5);

    );
   );
  );

  round(lambda*S)
};


Theta_all(ord)=
{
  my(K= posdisc(ord));

  for(i= 1,length(K),

    n=K[i][1];
    l=K[i][2];
    m=K[i][3];

    print([n,l,m]," : ",Theta(n,l,m));
  );
};

coeff(p, cfs, vars) = {
  my(l1 =length(cfs), l2 = length(vars)); 
  if(l1!=l2, error("Length mismatch"));
  if(l1==0, error("Empty list"));
  my(c = p);
  for(i= 1, n,
    c = polcoef(c, cfs[i], vars[i]);
  );
  c
};

coeff(poly,exps,vars)=
{
  my(c= poly);
  for(i= 1,#vars,
    if(type(c)!="t_POL",
      return(if(exps[i]==0,c,0))
    );
    if(variable(c)!=variable(vars[i]),
      return(if(exps[i]==0,c,0))
    );
    c=polcoef(c,exps[i]);
  );
  c
};




eisensteinqser(wt,lim) = {
    my(s = 0);
    for(n=0,lim, 
    for(l = -lim, lim, 
    for(m = 0, lim, 
    s+= C(n,l,m,wt)*'q^n*'y^l*p^m;
    );););
    return(s)
};

\\ for(n = 0, 4, for(l = 0,4, for(m = 0, 4, print([n,l,m],":",lb*coeff(f10, [n,l,m],[q,y,p])))))




x10(lim) = c*(eisensteinqser(4,lim)*eisensteinqser(6,lim)- eisensteinqser(10,2*lim) - polcoef(eisensteinqser(4,lim)*eisensteinqser(6,lim)- eisensteinqser(10,2*lim),0,q));


detsemipos(v)=
{
  my(n=v[1],l=v[2],m=v[3]);
  (n>=0 && m>=0 && 4*n*m-l^2>=0)
}

decompose2(v)=
{
  my(n=v[1],l=v[2],m=v[3],res=List());

  for(n1=0,n,
    for(m1=0,m,
      for(l1=-abs(l)-2*sqrtint(n*m), abs(l)+2*sqrtint(n*m),

        my(v1=[n1,l1,m1]);
        my(v2=[n-n1,l-l1,m-m1]);

        if(detsemipos(v1) && detsemipos(v2),
          listput(res,[v1,v2]);
        );
      )
    )
  );

  Vec(res)
}

CF10(v) =
{
    my(a = v[1], b = v[2], c = v[3]);
    my(c10 = C(a,b,c,10));
    my(dec = decompose2(v));
    my(ll = -43867/2307916800);
    my(s = 0);
    for(i=1,length(dec),
      my(v1 = dec[i][1]);
      my(v2 = dec[i][2]);
      s += C(v1[1],v1[2],v1[3],4)*C(v2[1],v2[2],v2[3],6);
    );
    return(ll*(c10 - s));
}



