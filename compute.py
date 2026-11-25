import numpy as np, time
from core import *
out={}
for m in ['annealed','quenched','cosine']:
    t0=time.time(); kk=kicks(m)
    # band material, forward-tracked (exact), counted per 4-px slot of x
    nb=R//4
    for k,(a,b) in BANDS.items():
        yy=np.arange(np.ceil(a),np.floor(b)+1); xx=np.arange(0,R,0.5)
        X,Y=np.meshgrid(xx,yy); st=E.new_state(X.ravel(),Y.ravel())
        n0=X.size/nb
        rho=np.zeros((61,nb)); rl=np.zeros((61,nb))
        rho[0]=np.bincount((np.mod(st['x'],R)//4).astype(int),minlength=nb)/n0; rl[0]=rho[0]
        for N in range(1,61):
            E.step(st,kk[N-1],Tau,row0)
            S=E.sigma_max(st['a'],st['b'],st['c'],st['d']); lo=S<=2; xb=(np.mod(st['x'],R)//4).astype(int)
            rho[N]=np.bincount(xb,minlength=nb)/n0; rl[N]=np.bincount(xb[lo],minlength=nb)/n0
        out[f'{m}_{k}_rho']=rho; out[f'{m}_{k}_low']=rl
    # share of the picture in each S class, every cycle (4-px lattice: unbiased)
    edges=np.array([2,5,10,100]); sh=np.zeros((61,5)); sh[0,0]=1
    for N in range(1,61):
        _,S,_=picture_and_S(m,N,4)
        c=np.digitize(S.ravel(),edges,right=True); sh[N]=np.bincount(c,minlength=5)/S.size
        if N in (20,60): out[f'{m}_lnS_{N}']=np.log(S.ravel())
        out.setdefault(f'{m}_mean',np.zeros(61))[N]=np.log(S).mean(); out.setdefault(f'{m}_var',np.zeros(61))[N]=np.log(S).var()
    out[f'{m}_share']=sh
    print(m,round(time.time()-t0))
np.savez('stats.npz',**out)
