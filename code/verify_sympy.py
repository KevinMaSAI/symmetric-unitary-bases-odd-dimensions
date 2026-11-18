#!/usr/bin/env python3
"""Second, independent exact check using sympy polynomial reduction modulo the ideal
(Phi5(W), S^2-(8W+8W^4-7)) (a Groebner basis, lex S>W). Conjugation: W->W^4, S->-S."""
import sympy as sp, itertools
W,S=sp.symbols('W S')
G=sp.groebner([W**4+W**3+W**2+W+1, S**2-(8*W+8*W**4-7)],S,W,order='lex',domain='QQ')
def red(e): return G.reduce(sp.expand(e))[1]
def cj(e): return sp.expand(e).subs({W:W**4,S:-S},simultaneous=True)
R=sp.Rational
mu2=R(-3,11)+R(-5,44)*W+R(3,11)*W**2+R(7,22)*W**3+S*(R(2,11)+R(7,44)*W-R(2,11)*W**2-R(1,22)*W**3)
v=R(-6,11)-R(21,44)*W-R(9,44)*W**2-R(27,44)*W**3+S*(R(4,11)+R(25,44)*W+R(17,44)*W**2+R(7,44)*W**3)
assert red(mu2*cj(mu2)-1)==0 and red(v*cj(v)-1)==0
mu=[1,W,mu2,mu2,W]; mup=[1,mu2,W,W,mu2]
k =[sum(mu[j]*W**(j*d) for j in range(5))/5 for d in range(5)]
kp=[sum(mup[j]*W**(j*d) for j in range(5))/5 for d in range(5)]
V1=sp.eye(5)+(v-1)/5*sp.ones(5,5)
V2=sp.Matrix(5,5,lambda a,b: W**((2*a*a+2*b*b)%5)*k[(b-a)%5])
V3=sp.Matrix(5,5,lambda a,b: W**((3*a*a+3*b*b)%5)*kp[(b-a)%5])
Ms=[sp.Matrix(5,5,lambda a,b: W**((t*(a+b))%5)*V[a,b]) for V in (V1,V2,V3) for t in range(5)]
Msc=[M.applyfunc(cj) for M in Ms]
ok=True
for i,M in enumerate(Ms):
    ok&= all(red(M[a,b]-M[b,a])==0 for a in range(5) for b in range(5))
    P=M*Msc[i].T
    ok&= all(red(P[a,b]-(1 if a==b else 0))==0 for a in range(5) for b in range(5))
print('symmetric & unitary:',ok)
fails=0
for i,j in itertools.combinations_with_replacement(range(15),2):
    t=red(sum(Ms[i][a,b]*Msc[j][a,b] for a in range(5) for b in range(5)))
    if t!=(5 if i==j else 0): fails+=1
print('Gram == 5I failures:',fails)
print('ALL PASSED' if ok and fails==0 else 'FAILED')
