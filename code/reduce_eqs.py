# symbolic check of the two reduced equations E12, E23 (Lemma 3 of proof.md)
import sympy as sp
W,m,v=sp.symbols('W mu v')   # conj(W)=1/W, conj(mu)=1/mu, conj(v)=1/v on the unit circle
phi=W**4+W**3+W**2+W+1
mus=[1,W,m,m,W]
k=[sum(mus[j]*W**(j*d) for j in range(5))/5 for d in range(5)]
cj=lambda e: e.subs({W:1/W,m:1/m,v:1/v},simultaneous=True)
y=[(4+v)/5]+[(v-1)/5]*4
E12=sum(y[d]*W**(-d*d)*cj(k[d]) for d in range(5))
E23=sum(W**(2*d*d)*k[d]*cj(k[(3*d)%5]) for d in range(5))
def rednum(e):
    n,dn=sp.fraction(sp.together(sp.expand(e)))
    # clear powers of W in denominator then reduce mod phi
    return sp.rem(sp.expand(n),phi,W), sp.factor(dn)
n12,d12=rednum(E12); n23,d23=rednum(E23)
a=-4*W**4+4*W**3+4*W+6; b=2-2*W**4-2*W**3+2*W**2+5*W; c=4*W**4-4*W**3+6*W**2+4*W
c1=(-2*W**4*m + 4*W**3 - 2*W**2*m + 2*W**2 - 3*W*m + 4*W + 2*m); c0=2*W**4*m - 4*W**3 + 2*W**2*m - 2*W**2 + 8*W*m + 6*W + 8*m
print('den12',d12,' den23',d23)
for nm,num,target in [('E12',n12,c1*v+c0),('E23',n23,a*m**2+b*m+c)]:
    # find unit u in Q(W)*W^k with num = u*target mod phi
    P=sp.Poly(num,v,m); T=sp.Poly(sp.expand(target),v,m)
    mono=T.monoms()[0]; ratio=None
    for kk in range(5):
        for sgn in [1,-1]:
            for scal in [sp.Rational(p,q) for p in range(1,9) for q in [1,2,4,5,10,25]]:
                diff=sp.rem(sp.expand(num-sgn*scal*W**kk*target),phi,W)
                if sp.expand(diff)==0: ratio=(sgn*scal,kk)
    print(nm,'numerator = const*W^k*target :',ratio)
