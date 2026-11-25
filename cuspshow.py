import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from core import *
plt.rcParams.update({'font.size':12})
kk=kicks('annealed'); a0,b0=BANDS['red']; yc=0.5*(a0+b0)
def kick(st,dp):
    x=st['x']; xf=np.floor(x); i0=np.mod(xf.astype(np.int64),R); d0=dp[i0]; fp=dp[(i0+1)%R]-d0
    st['y']=st['y']+d0+fp*(x-xf); st['c']=st['c']+fp*st['a']; st['d']=st['d']+fp*st['b']
def drift(st,t):
    o=dict(st); o['x']=st['x']+t*Tau*(st['y']-row0); o['a']=st['a']+t*Tau*st['c']; o['b']=st['b']+t*Tau*st['d']; return o
# centre line (thin) and whole band (thick), seed points
D=64; x0=np.arange(R*D)/D; L=E.new_state(x0,np.full(x0.size,yc))
yy=np.arange(np.ceil(a0),np.floor(b0)+1,2.0); XX,YY=np.meshgrid(np.arange(0,R,1.0),yy); B=E.new_state(XX.ravel(),YY.ravel())
piece=(x0>850)&(x0<1050)
frames=[]; dens=[]; tt=[]; folds=[]
nb=R//4
for N in range(1,8):
    kick(L,kk[N-1]); kick(B,kk[N-1])
    for t in np.linspace(0,1,21)[1:] if N>1 else [1.0]:
        Lt=drift(L,t); h=np.bincount((np.mod(Lt['x'],R)//4).astype(int),minlength=nb)/(D*4)
        dens.append(h); tt.append(N-1+t)
        a_=Lt['a'][piece]; s_=np.nonzero(a_[:-1]*a_[1:]<=0)[0]; folds.append((N-1+t,np.mod(Lt['x'][piece][s_],R)))
        if 2<=N<=3: frames.append((N-1+t,Lt,drift(B,t)))
    L=drift(L,1.0); B=drift(B,1.0)
# exact cusp moment for the piece: first t where a hits 0 on the piece
# refine: during cycle 3's drift, a(t) = a + t*Tau*c is linear in t; find first t where min over the piece of a hits 0
L2=[f for f in frames if abs(f[0]-2.05)<1e-6][0][1]   # any frame in cycle 3 carries the post-kick state; rebuild from it
base={k:v.copy() for k,v in L2.items()}; base['x']=base['x']-0.05*Tau*(base['y']-row0); base['a']=base['a']-0.05*Tau*base['c']
ap=base['a'][piece]; cp=base['c'][piece]
tg=np.linspace(0,1,4001); mins=[np.min(ap+t*Tau*cp) for t in tg]
t_star=tg[np.argmax(np.array(mins)<=0)]; tc=2+t_star
star_state=drift(base,t_star)
B2=[f for f in frames if abs(f[0]-2.05)<1e-6][0][2]
bb={k:v.copy() for k,v in B2.items()}; bb['x']=bb['x']-0.05*Tau*(bb['y']-row0); bb['a']=bb['a']-0.05*Tau*bb['c']
frames.append((tc,star_state,drift(bb,t_star)))
print('cusp at cycle', tc)
# choose snapshots: before, at, after
want=[2.0, round(tc-0.15,2), tc, round(tc+0.2,2), 3.0]
frames.sort(key=lambda f:f[0])
fig=plt.figure(figsize=(19,15))
g=fig.add_gridspec(2,5,height_ratios=[1,1.15],left=0.05,right=0.98,top=0.86,bottom=0.06,hspace=0.32,wspace=0.12)
for j,w in enumerate(want):
    tv,Lt,Bt=min(frames,key=lambda f:abs(f[0]-w))
    ax=fig.add_subplot(g[0,j])
    bx=np.mod(Bt['x'],R); m=(bx>560)&(bx<900)
    ax.plot(bx[m],Bt['y'][m],',',color='#f2a0a0')
    lx=np.mod(Lt['x'][piece],R); ly=Lt['y'][piece]; a=Lt['a'][piece]
    ax.plot(lx,ly,'-',color='#a00000',lw=2.2)
    s=np.nonzero(a[:-1]*a[1:]<=0)[0]
    if abs(tv-tc)<1e-9:
        k=np.argmin(np.abs(a)); ax.plot(lx[k],ly[k],'*',ms=24,mfc='cyan',mec='k',mew=2)
    elif len(s): ax.plot(lx[s],ly[s],'o',ms=11,mfc='none',mec='k',mew=2)
    else:
        k=np.argmin(np.abs(a)); ax.plot(lx[k],ly[k],'o',ms=11,mfc='none',mec='0.5',mew=2,ls='')
    ax.set_xlim(600,880); ax.set_ylim(900,730); ax.set_xlabel('x (px)')
    if j==0: ax.set_ylabel('row (direction)')
    lab={0:'leaning (grey ring = steepest point)',1:'almost upright',2:'EXACTLY UPRIGHT AT ONE POINT\n= THE CUSP',3:'tipped over: two upright points\n= two caustics',4:'the overhang grows'}[j]
    ax.set_title(f'cycle {tv:.2f}\n{lab}',fontsize=12,fontweight='bold' if j==2 else 'normal')
fig.text(0.05,0.93,'Top: the red band (pale) and its centre line (dark red) while the light drifts from cycle 2 to cycle 3, shown at in-between moments.\n'
 'The band slowly stands up. The cusp (star) is the single instant and place where it first becomes exactly upright. Right after, it overhangs and has two upright points: the two caustics.',fontsize=12.5,va='top')
ax=fig.add_subplot(g[1,:])
H=np.array(dens); H[H<=0]=np.nan
cm=LinearSegmentedColormap.from_list('r',['#000000','#d62728','#ffffff'])
ax.imshow(H,origin='lower',aspect='auto',extent=[0,R,tt[0],tt[-1]],cmap=cm,norm=LogNorm(0.05,50),interpolation='nearest')
tv,Lt,Bt=min(frames,key=lambda f:abs(f[0]-tc)); a=Lt['a'][piece]; k=np.argmin(np.abs(a)); xs=np.mod(Lt['x'][piece][k],R)
ax.plot([xs],[tc],'*',ms=26,mfc='cyan',mec='k',mew=2)
ax.annotate('THE CUSP: the tip where this branch starts',xy=(xs,tc),xytext=(820,1.35),color='cyan',fontsize=13,fontweight='bold',arrowprops=dict(arrowstyle='->',color='cyan',lw=2))
for w in want: ax.axhline(w,color='w',ls=':',lw=.8)
fl=[(t,np.min(x),np.max(x)) for t,x in folds if t>=tc and t<=4.2 and len(x)>=2]
fl=np.array(fl); ax.plot(fl[:,1],fl[:,0],'c--',lw=2.2); ax.plot(fl[:,2],fl[:,0],'c--',lw=2.2)
ax.annotate('the two caustics (dashed) = the two edges of the branch',xy=(fl[-1,1],fl[-1,0]),xytext=(450,4.9),color='cyan',fontsize=12.5,fontweight='bold',arrowprops=dict(arrowstyle='->',color='cyan',lw=2))
ax.set_xlim(430,1100); ax.set_ylim(1.0,7); ax.set_xlabel('x (px)'); ax.set_ylabel('cycle (including in-between moments)')
ax.set_title('Bottom: the same, as the stacked map, but with in-between moments so the rows are thin. The bright V above the star is the branch; its two sides are the two caustics.',fontsize=12)
fig.suptitle('The cusp, shown directly (red band, annealed, seed 0)',fontsize=17,fontweight='bold',y=0.985)
fig.savefig('cusp_shown.png',dpi=78)
