import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, gc
from core import *
plt.rcParams.update({'font.size':12})
Ns=[3,13,40]; nb=R//4; xc=np.arange(nb)*4+2
def band_states(mode,Nlist):
    kk=kicks(mode); pts=[]
    for k,(a,b) in [('green',BANDS['green'])]:
        yy=np.arange(np.ceil(a),np.floor(b)+1,1.0); X,Y=np.meshgrid(np.arange(0,R,0.5),yy); pts.append((X.ravel(),Y.ravel()))
    x=np.concatenate([p[0] for p in pts]); y=np.concatenate([p[1] for p in pts]); st=E.new_state(x,y); n0=x.size/nb; out={}
    for N in range(1,max(Nlist)+1):
        E.step(st,kk[N-1],Tau,row0)
        if N in Nlist: out[N]=(np.mod(st['x'],R).copy(),st['y'].copy())
    return out,n0
fig=plt.figure(figsize=(20,17.5))
g=fig.add_gridspec(4,3,height_ratios=[1,0.32,1,0.32],left=0.08,right=0.98,top=0.90,bottom=0.05,hspace=0.18,wspace=0.06)
NAME={'annealed':'ANNEALED','quenched':'QUENCHED'}
for i,m in enumerate(['annealed','quenched']):
    bs,n0=band_states(m,Ns)
    for j,N in enumerate(Ns):
        x,y=bs[N]; h=np.bincount((x//4).astype(int),minlength=nb)/n0; br=h>=2.0
        onb=br[(x//4).astype(int)]
        pic,S,_=picture_and_S(m,N,1)
        ax=fig.add_subplot(g[2*i,j]); ax.imshow(overlay(pic,S),interpolation='nearest',extent=[0,R,R,0])
        yy=np.mod(y,R); ax.plot(x[onb][::3],yy[onb][::3],',',color='#00e5ff',alpha=0.9)
        for b0 in np.nonzero(np.diff(np.r_[0,br.astype(int),0]))[0][::2]: pass
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0,R); ax.set_ylim(1500,420)
        ax.set_title(f'{NAME[m]}, cycle {N}: green band: the part IN a branch is cyan',fontsize=12)
        if j==0: ax.set_ylabel('row (direction)')
        ax2=fig.add_subplot(g[2*i+1,j]); ax2.semilogy(xc,h,'k-',lw=1); ax2.fill_between(xc,0.05,h,where=br,color='#00e5ff',step='mid')
        ax2.axhline(2,color='#00a0b0',ls='--',lw=1); ax2.axhline(1,color='k',ls=':',lw=1); ax2.set_ylim(0.05,30); ax2.set_xlim(0,R)
        ax2.set_xlabel('x (px)'); 
        if j==0: ax2.set_ylabel('green band\nper slot')
        ax2.set_title(f'branches at cycle {N} = slots with at least 2x (cyan): {100*br.mean():.0f}% of x',fontsize=10.5)
        del pic,S; gc.collect()
fig.suptitle('Seeing branches in the kick-drift picture itself',fontsize=17,fontweight='bold',y=0.975)
fig.text(0.08,0.935,'Under each picture: how much of the GREEN band sits above each slot of x; cyan slots are branches (at least 2x the even amount).\n'
 'In the picture, every bit of green band lying over a cyan slot is coloured cyan: that is the material making the branches. Look at its SHAPE.',fontsize=12.5,va='top')
fig.savefig('branches_in_picture.png',dpi=70)
