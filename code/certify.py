#!/usr/bin/env python3
"""Rigorous (interval-arithmetic, Krawczyk) existence certificate for the "chirped circulant" construction.
For odd n, h=(n-1)/2, m=h+1, seeds i=0..h with chirp parameter c_i=i, unknown phases theta[i][r], r=1..h
(theta[i][0]=0). mu^i_l = exp(i*theta[i][fold(l)]), fold(l)=min(l,n-l).
k^i_d = (1/n) sum_l mu^i_l w^{l d},  w=exp(2 pi i/n).
Equations (i<j):  E_ij = sum_d w^{e_ij d^2} k^i_d conj(k^j_d) = 0,   e_ij=(c_i-c_j)*2^{-1} mod n.
Square real system F: R^{m h} -> R^{m h}.  Krawczyk: if K(X) is contained in int(X) then F has a unique zero in X.
Usage: python3 certify.py n [seed]        (search + certify)
       python3 certify.py n --from-center (re-certify archived centre)"""
import sys, numpy as np, time
from flint import arb, acb, arb_mat, ctx
from scipy.optimize import least_squares
n=int(sys.argv[1]); seed=int(sys.argv[2]) if len(sys.argv)>2 and sys.argv[2].isdigit() else 0
h=(n-1)//2; m=h+1; cs=list(range(m)); inv2=(n+1)//2
fold=[min(l,n-l) for l in range(n)]
pairs=[(i,j) for i in range(m) for j in range(i+1,m)]; N=m*h
# ---------- float model ----------
w=np.exp(2j*np.pi/n); Wld=np.array([[w**((l*d)%n) for d in range(n)] for l in range(n)])
d2=np.array([(d*d)%n for d in range(n)])
E=[[None]*m for _ in range(m)]
wt={p: w**(((cs[p[0]]-cs[p[1]])*inv2*d2)%n) for p in pairs}
def unpack(x): T=np.zeros((m,h+1)); T[:,1:]=x.reshape(m,h); return T
def F_float(x):
    T=unpack(x); MU=np.exp(1j*T[:,fold]); K=MU@Wld/n; r=[]
    for (i,j) in pairs: e=np.sum(wt[(i,j)]*K[i]*K[j].conjugate()); r+=[e.real,e.imag]
    return np.array(r)
def J_float(x):
    T=unpack(x); MU=np.exp(1j*T[:,fold]); K=MU@Wld/n
    dK=np.zeros((m,h,n),complex)   # dK[i,r-1,d] = d k^i_d / d theta_{i,r}
    for l in range(1,n): dK[:,fold[l]-1,:]+=1j*MU[:,l,None]*Wld[l][None,:]/n
    J=np.zeros((len(pairs)*2,N))
    for q,(i,j) in enumerate(pairs):
        a=(wt[(i,j)][None,:]*dK[i]*K[j].conjugate()[None,:]).sum(1); b=(wt[(i,j)][None,:]*K[i][None,:]*dK[j].conjugate()).sum(1)
        J[2*q,i*h:(i+1)*h]=a.real; J[2*q+1,i*h:(i+1)*h]=a.imag; J[2*q,j*h:(j+1)*h]=b.real; J[2*q+1,j*h:(j+1)*h]=b.imag
    return J
t0=time.time(); rng=np.random.default_rng(seed)
LOAD='--from-center' in sys.argv   # re-certify the archived centre in centers/cert_sol_n{n}.npy (no search)
for trial in ([] if LOAD else range(200)):
    s=least_squares(F_float,rng.uniform(0,2*np.pi,N),jac=J_float,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=5000)
    if np.abs(F_float(s.x)).max()<1e-12: break
if LOAD:
    import os; x=np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)),'centers',f'cert_sol_n{n}.npy')); trial='archived'
else:
    x=s.x
    for _ in range(5): x=x-np.linalg.solve(J_float(x),F_float(x))
print(f'n={n}: float solution found (trial {trial}), max|F|={np.abs(F_float(x)).max():.2e}, cond(J)={np.linalg.cond(J_float(x)):.2e}, {time.time()-t0:.1f}s',flush=True)
None if LOAD else np.save(f'cert_sol_n{n}.npy',x)
# ---------- interval model ----------
ctx.prec=128
pi2=arb.pi()*2
Wa=[acb(arb(2*k)/n).exp_pi_i() for k in range(n)]           # w^k, rigorous
def F_arb(X):   # X: list of arb balls
    T=[[arb(0)]+X[i*h:(i+1)*h] for i in range(m)]
    MU=[[acb(0,T[i][fold[l]]).exp() for l in range(n)] for i in range(m)]
    K=[[sum((MU[i][l]*Wa[(l*d)%n] for l in range(n)),acb(0))/n for d in range(n)] for i in range(m)]
    return MU,K
def eval_F(X):
    MU,K=F_arb(X); out=[]
    for (i,j) in pairs:
        e=sum((Wa[((cs[i]-cs[j])*inv2*d*d)%n]*K[i][d]*K[j][d].conjugate() for d in range(n)),acb(0)); out+=[e.real,e.imag]
    return out
def eval_J(X):
    MU,K=F_arb(X)
    dK=[[[sum((acb(0,1)*MU[i][l]*Wa[(l*d)%n] for l in range(1,n) if fold[l]==r),acb(0))/n for d in range(n)] for r in range(1,h+1)] for i in range(m)]
    J=[[arb(0)]*N for _ in range(2*len(pairs))]
    for q,(i,j) in enumerate(pairs):
        wd=[Wa[((cs[i]-cs[j])*inv2*d*d)%n] for d in range(n)]
        Kjc=[K[j][d].conjugate() for d in range(n)]
        for r in range(h):
            a=sum((wd[d]*dK[i][r][d]*Kjc[d] for d in range(n)),acb(0))
            b=sum((wd[d]*K[i][d]*dK[j][r][d].conjugate() for d in range(n)),acb(0))
            J[2*q][i*h+r]=a.real; J[2*q+1][i*h+r]=a.imag; J[2*q][j*h+r]=b.real; J[2*q+1][j*h+r]=b.imag
    return J
rad=1e-9
xt=[arb(float(v)) for v in x]                      # exact midpoint (binary64 values)
X=[arb(float(v),rad) for v in x]
Fx=eval_F(xt); JX=eval_J(X)
Y=np.linalg.inv(J_float(x)); Ya=arb_mat([[arb(float(v)) for v in row] for row in Y])
JXa=arb_mat(JX); I=arb_mat([[arb(1) if a==b else arb(0) for b in range(N)] for a in range(N)])
D=arb_mat([[arb(0,rad)] for _ in range(N)])        # X - xt
Kv=arb_mat([[v] for v in xt]) - Ya*arb_mat([[v] for v in Fx]) + (I-Ya*JXa)*D
ok=all(X[a].contains_interior(Kv[a,0]) for a in range(N))
maxdev=max(float((Kv[a,0]-xt[a]).abs_upper()) for a in range(N))
print(f'n={n}: Krawczyk box radius {rad:g}, max|K-x~|={maxdev:.3e}  ->  {"CERTIFIED: unique exact zero in box" if ok and maxdev<rad else "NOT certified"}  ({time.time()-t0:.1f}s)',flush=True)
