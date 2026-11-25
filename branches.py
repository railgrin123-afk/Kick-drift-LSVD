import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from scipy import ndimage as ndi
from core import *
plt.rcParams.update({'font.size':12})
SUB=20; NMAX=40; nb=R//4
def kick(st,dp):
    x=st['x']; xf=np.floor(x); i0=np.mod(xf.astype(np.int64),R); d0=dp[i0]; fp=dp[(i0+1)%R]-d0
    st['y']=st['y']+d0+fp*(x-xf); st['c']=st['c']+fp*st['a']; st['d']=st['d']+fp*st['b']
def fine_map(mode,band='red'):
    kk=kicks(mode); a0,b0=BANDS[band]
    yy=np.arange(np.ceil(a0),np.floor(b0)+1,2.0); X,Y=np.meshgrid(np.arange(0,R,0.5),yy); st=E.new_state(X.ravel(),Y.ravel())
    n0=st['x'].size/nb; H=[np.bincount((np.mod(st['x'],R)//4).astype(int),minlength=nb)/n0]; T=[0.0]
    for N in range(1,NMAX+1):
        kick(st,kk[N-1]); y=st['y']-row0
        for j in range(1,SUB+1):
            x=st['x']+(j/SUB)*Tau*y; H.append(np.bincount((np.mod(x,R)//4).astype(int),minlength=nb)/n0); T.append(N-1+j/SUB)
        st['x']=st['x']+Tau*y
    return np.array(H),np.array(T)
out={}
for m in ['annealed','quenched']:
    H,T=fine_map(m); bright=H>=2.0
    lab,n=ndi.label(bright,structure=np.ones((3,3)))
    info=[]
    for k in range(1,n+1):
        ii,jj=np.nonzero(lab==k); t0,t1=T[ii.min()],T[ii.max()]
        if t1-t0>=1.0: info.append((k,t0,t1,jj[ii.argmin()]*4+2,float(H[ii,jj].max())))
    out[m]=(H,T,lab,info)
    L=np.array([i[2]-i[1] for i in info])
    print(m,'branches lasting >=1 cycle:',len(info),' median life %.1f  longest %.1f'%(np.median(L),L.max()),
          ' alive at cycle 10/20/30/39:',[sum(1 for i in info if i[1]<=c<=i[2]) for c in (10,20,30,39)])
np.save('branch_info.npy',{m:(out[m][3]) for m in out},allow_pickle=True)
fig=plt.figure(figsize=(19,16))
g=fig.add_gridspec(2,2,height_ratios=[2.2,1],left=0.06,right=0.98,top=0.84,bottom=0.06,hspace=0.30,wspace=0.16)
cm=LinearSegmentedColormap.from_list('r',['#000000','#d62728','#ffffff'])
NAME={'annealed':'ANNEALED (new kick every cycle)','quenched':'QUENCHED (one kick, repeated)'}
for j,m in enumerate(['annealed','quenched']):
    H,T,lab,info=out[m]; ax=fig.add_subplot(g[0,j]); Hp=H.copy(); Hp[Hp<=0]=np.nan
    ax.imshow(Hp,origin='lower',aspect='auto',extent=[0,R,T[0],T[-1]],cmap=cm,norm=LogNorm(0.05,50),interpolation='nearest')
    keep=np.isin(lab,[i[0] for i in info])
    ax.contour(np.arange(nb)*4+2,T,keep.astype(float),[0.5],colors='cyan',linewidths=1.1)
    for k,t0,t1,x0b,mx in info: ax.plot(x0b,t0,'*',ms=11,mfc='yellow',mec='k',mew=0.8)
    ax.set_title(f'{NAME[m]}: red band, every moment stacked\ncyan outline = a BRANCH (at least 2x the even amount, lasting at least 1 cycle);  yellow star = where it starts',fontsize=12)
    ax.set_xlabel('x (px)'); ax.set_ylim(0,NMAX)
    if j==0: ax.set_ylabel('cycle (continuous)')
    ax2=fig.add_subplot(g[1,j]); L=np.array([i[2]-i[1] for i in info]); st_=np.array([i[1] for i in info])
    sc=ax2.scatter(st_,L,c=[min(i[4],15) for i in info],cmap='autumn_r',s=40,edgecolors='k',linewidths=.4)
    ax2.set_xlabel('cycle the branch starts'); ax2.set_ylabel('how long it lasts (cycles)'); ax2.set_xlim(0,NMAX); ax2.set_ylim(0,NMAX)
    ax2.set_title(f'{len(info)} branches; median lifetime {np.median(L):.1f} cycles, longest {L.max():.1f}',fontsize=12)
    cb=fig.colorbar(sc,ax=ax2,fraction=0.04); cb.set_label('peak brightness (x even)')
fig.suptitle('How to tell a branch: a bright ridge that keeps going. Each outline is one branch, found automatically.',fontsize=16,fontweight='bold',y=0.975)
fig.text(0.06,0.935,'Rule: a branch is a connected streak where the band is at least twice as concentrated as even spreading, lasting at least one cycle. Its start (star) is its cusp.\n'
         'Bottom: every branch found, by when it starts and how long it lasts.',fontsize=12.5,va='top')
fig.savefig('branches_found.png',dpi=72)
