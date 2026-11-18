#!/usr/bin/env python3
"""Exact verification (pure Python, fractions.Fraction only) of an orthogonal basis of
15 symmetric unitary 5x5 matrices.  Arithmetic in R = Q(w)[S]/(S^2 - sigma),
w = primitive 5th root of unity, sigma = 8w + 8w^4 - 7 = 4*sqrt5 - 11 (<0).
Embedding into C: w -> exp(2 pi i/5), S -> i*sqrt(11 - 4 sqrt5).
Complex conjugation corresponds to the ring automorphism w -> w^4, S -> -S."""
from fractions import Fraction as Fr
import itertools, cmath, math

class R:
    # element = (A, B) meaning A + B*S, A,B in Q[w]/(Phi5) stored as 5 coefficients of w^0..w^4
    __slots__=('a','b')
    def __init__(s,a,b=None):
        s.a=list(a); s.b=list(b) if b is not None else [Fr(0)]*5
    @staticmethod
    def pm(p,q):
        r=[Fr(0)]*5
        for i in range(5):
            if p[i]:
                for j in range(5):
                    if q[j]: r[(i+j)%5]+=p[i]*q[j]
        return r
    @staticmethod
    def padd(p,q): return [x+y for x,y in zip(p,q)]
    def __add__(s,o): o=lift(o); return R(R.padd(s.a,o.a),R.padd(s.b,o.b))
    __radd__=__add__
    def __neg__(s): return R([-x for x in s.a],[-x for x in s.b])
    def __sub__(s,o): return s+(-lift(o))
    def __rsub__(s,o): return lift(o)-s
    def __mul__(s,o):
        o=lift(o)
        A=R.padd(R.pm(s.a,o.a),R.pm(R.pm(s.b,o.b),SIGMA))
        B=R.padd(R.pm(s.a,o.b),R.pm(s.b,o.a))
        return R(A,B)
    __rmul__=__mul__
    def conj(s):
        return R([s.a[(-i)%5] for i in range(5)],[-s.b[(-i)%5] for i in range(5)])
    def iszero(s):
        # Q[x]/(x^5-1) element is 0 in Q(w) iff all coefficients are equal
        return len(set(s.a))==1 and len(set(s.b))==1
    def num(s):
        w=cmath.exp(2j*math.pi/5); Sv=1j*math.sqrt(11-4*math.sqrt(5))
        return sum(float(s.a[i])*w**i for i in range(5))+Sv*sum(float(s.b[i])*w**i for i in range(5))
def lift(x):
    if isinstance(x,R): return x
    return R([Fr(x),0,0,0,0])
SIGMA=[Fr(-7),Fr(8),Fr(0),Fr(0),Fr(8)]
def Wp(k): c=[Fr(0)]*5; c[k%5]=Fr(1); return R(c)
def elt(a_coeffs,b_coeffs):  # coefficients of w^0..w^3
    return R([Fr(x) for x in a_coeffs]+[Fr(0)],[Fr(x) for x in b_coeffs]+[Fr(0)])

# --- the exact parameters ---
# mu2 = -3/11 - 5/44 w + 3/11 w^2 + 7/22 w^3 + S*(2/11 + 7/44 w - 2/11 w^2 - 1/22 w^3)
mu2=elt([Fr(-3,11),Fr(-5,44),Fr(3,11),Fr(7,22)],[Fr(2,11),Fr(7,44),Fr(-2,11),Fr(-1,22)])
# v = -6/11 - 21/44 w - 9/44 w^2 - 27/44 w^3 + S*(4/11 + 25/44 w + 17/44 w^2 + 7/44 w^3)
v=elt([Fr(-6,11),Fr(-21,44),Fr(-9,44),Fr(-27,44)],[Fr(4,11),Fr(25,44),Fr(17,44),Fr(7,44)])
w=Wp(1)
Z=lambda: lift(0)
def circ_entry(mu,d):  # (1/5) sum_j mu_j w^{j d}
    t=Z()
    for j in range(5): t=t+mu[j]*Wp(j*d)
    return t*Fr(1,5)
mu =[lift(1),w,mu2,mu2,w]
mup=[lift(1),mu2,w,w,mu2]
V1=[[ (lift(1) if a==b else Z()) + (v-1)*Fr(1,5) for b in range(5)] for a in range(5)]
V2=[[ Wp(2*a*a+2*b*b)*circ_entry(mu ,(b-a)%5) for b in range(5)] for a in range(5)]
V3=[[ Wp(3*a*a+3*b*b)*circ_entry(mup,(b-a)%5) for b in range(5)] for a in range(5)]
Ms={}
for i,V in enumerate([V1,V2,V3],1):
    for k in range(5):
        Ms[(i,k)]=[[Wp(k*(a+b))*V[a][b] for b in range(5)] for a in range(5)]
labels=sorted(Ms)
ok=True
# unimodularity of parameters
for nm,z in [('mu2',mu2),('v',v)]:
    t=z*z.conj()-1; print(f'|{nm}|^2 - 1 == 0 :',t.iszero()); ok&=t.iszero()
# symmetry and unitarity
for L in labels:
    M=Ms[L]
    sym=all((M[a][b]-M[b][a]).iszero() for a in range(5) for b in range(5))
    uni=True
    for a in range(5):
        for b in range(5):
            t=Z()
            for c in range(5): t=t+M[a][c]*M[b][c].conj()
            if not (t-(1 if a==b else 0)).iszero(): uni=False
    ok&=sym and uni
    print(L,'symmetric:',sym,' unitary (M M^* = I):',uni)
# Hilbert-Schmidt orthogonality
bad=0
for P,Q in itertools.combinations_with_replacement(labels,2):
    A,B=Ms[P],Ms[Q]; t=Z()
    for a in range(5):
        for b in range(5): t=t+A[a][b]*B[a][b].conj()
    target=5 if P==Q else 0
    if not (t-target).iszero(): bad+=1; print('FAIL',P,Q)
print('Gram matrix tr(M_P M_Q^*) == 5*delta_PQ for all 120 pairs (incl. diagonal); failures =',bad)
ok&=bad==0

# ---- checks of the hand-proof lemmas (proof.md, Sec. 4) ----
def P(*cs): return R([Fr(x) for x in cs]+[Fr(0)]*(5-len(cs)))
A=P(6,4,0,4,-4); B=P(2,5,2,-2,-2); C=P(0,4,6,-4,4)
sqrt5=P(1,2,0,0,2)                      # 1+2(w+w^4) = sqrt5
Sel=R([Fr(0)]*5,[Fr(1),0,0,0,0])        # the element S
chk={}
chk['S^2 = 4*sqrt5-11']=(Sel*Sel-(sqrt5*4-11)).iszero()
chk['C = w^2*conj(A)']=(C-w*w*A.conj()).iszero()
chk['B = w^2*conj(B)']=(B-w*w*B.conj()).iszero()
chk['B^2-4AC = (5 w S)^2']=(B*B-A*C*4-(w*Sel*5)*(w*Sel*5)).iszero()
chk['|A|^2 = 80-20 sqrt5']=(A*A.conj()-(80-sqrt5*20)).iszero()
chk['|B|^2 = 45+20 sqrt5']=(B*B.conj()-(sqrt5*20+45)).iszero()
chk['2A*mu = -B + 5 w S']=(A*mu2*2-(-B+w*Sel*5)).iszero()
chk['A mu^2 + B mu + C = 0']=(A*mu2*mu2+B*mu2+C).iszero()
C1=P(0,4,2,4,0)+mu2*P(2,-3,-2,0,-2); C0=P(0,6,-2,-4,0)+mu2*P(8,8,2,0,2)
chk['C1 v + C0 = 0']=(C1*v+C0).iszero()
chk['|C0|^2 = |C1|^2']=(C0*C0.conj()-C1*C1.conj()).iszero()
for kk,val in chk.items(): print('lemma check:',kk,val); ok&=val
# floating sanity: embedding values, rank
import numpy as np
Mn=[np.array([[Ms[L][a][b].num() for b in range(5)] for a in range(5)]) for L in labels]
print('numeric mu2 =',mu2.num(),' v =',v.num())
print('numeric rank of the 15 matrices as vectors in C^25:',np.linalg.matrix_rank(np.array([m.ravel() for m in Mn])))
print('ALL EXACT CHECKS PASSED' if ok else 'SOME CHECK FAILED')
