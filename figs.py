import numpy as np, matplotlib, sys; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from core import *
plt.rcParams.update({'font.size':12})
D=np.load('stats.npz')
NAME={'annealed':'ANNEALED (new random kick every cycle)','quenched':'QUENCHED (one random kick, repeated)','cosine':'PERIODIC (sine kick, standard map)'}
BC={'red':'#d62728','green':'#2ca02c','blue':'#1f4fd6'}
cache={}
def get(m,N):
    if (m,N) not in cache: cache[(m,N)]=picture_and_S(m,N,1)
    return cache[(m,N)]
def show(ax,m,N,fade=False,title=None):
    pic,S,_=get(m,N); im=overlay(pic,S)
    if fade: im=np.where((S<=2)[...,None],im,0.94)
    ax.imshow(im,interpolation='nearest'); ax.contour(S,[2],colors='w' if not fade else 'k',linewidths=0.6)
    ax.set_xticks([]); ax.set_yticks([])
    if title: ax.set_title(title,fontsize=12)
from matplotlib.colors import ListedColormap
BLEND=ListedColormap(0.72*np.ones((256,3))+0.28*CMAP(np.linspace(0,1,256))[:,:3])
def colour_key(fig,rect):
    ax=fig.add_axes(rect); sm=ScalarMappable(norm=LogNorm(1,15),cmap=BLEND)
    cb=fig.colorbar(sm,cax=ax,orientation='horizontal'); cb.set_ticks([1,2,5,10,15]); cb.set_ticklabels(['1','2','5','10','15+'])
    cb.set_label('overlay colour = S, how stretched the material at that spot is   (grey-teal = barely, pale yellow = a lot; shown as it looks on the white picture)')
def fig1():
    Ns=[0,3,13,40]; modes=['annealed','quenched','cosine']
    fig,axs=plt.subplots(3,4,figsize=(20,16))
    fig.subplots_adjust(left=0.16,right=0.99,top=0.90,bottom=0.09,wspace=0.04,hspace=0.12)
    for i,m in enumerate(modes):
        for j,N in enumerate(Ns):
            _,S,_=get(m,N)
            show(axs[i,j],m,N,title=f'cycle {N}:  {100*(S<=2).mean():.1f}% of the picture has S ≤ 2')
        axs[i,0].text(-0.06,0.5,NAME[m].replace(' (','\n('),transform=axs[i,0].transAxes,ha='right',va='center',fontsize=13,fontweight='bold')
    fig.suptitle('Figure 1.  Your kick-drift picture with the stretch overlay, for the three kinds of kick',fontsize=17,fontweight='bold',y=0.975)
    fig.text(0.16,0.925,'Same starting picture (dots + red/green/blue bands) and same map in every row; only the kick differs. White line = the S = 2 boundary: inside it is low-stretch.',fontsize=12.5)
    colour_key(fig,[0.3,0.045,0.5,0.015])
    fig.savefig('F1_picture.png',dpi=72)
def fig2(m):
    Ns=[3,13,40]; xc=np.arange(R//4)*4+2
    fig=plt.figure(figsize=(19,22))
    g=fig.add_gridspec(3,3,height_ratios=[1.35,0.55,1.6],left=0.06,right=0.98,top=0.92,bottom=0.05,hspace=0.28,wspace=0.10)
    for j,N in enumerate(Ns):
        ax=fig.add_subplot(g[0,j]); show(ax,m,N,title=f'the picture at cycle {N}')
        ax2=fig.add_subplot(g[1,j])
        for k in BANDS: ax2.semilogy(xc,D[f'{m}_{k}_rho'][N],color=BC[k],lw=1)
        ax2.axhline(1,color='k',ls=':'); ax2.set_ylim(0.05,60); ax2.set_xlim(0,R)
        ax2.set_title(f'cycle {N}: how much of each band sits above each 4-px slot of x',fontsize=10.5)
        ax2.set_xlabel('x (px)')
        if j==0: ax2.set_ylabel('amount of band\n(1 = evenly spread)')
    for j,k in enumerate(BANDS):
        ax=fig.add_subplot(g[2,j]); cm=LinearSegmentedColormap.from_list(k,['#000000',BC[k],'#ffffff'])
        r=D[f'{m}_{k}_rho'].copy(); r[r<=0]=np.nan
        im=ax.imshow(r,origin='lower',aspect='auto',extent=[0,R,-0.5,60.5],cmap=cm,norm=LogNorm(0.05,50),interpolation='nearest')
        ax.set_facecolor('k')
        for N in Ns: ax.axhline(N,color='w',ls='--',lw=.8); ax.text(R*0.98,N+0.5,f'cycle {N}',color='w',ha='right',fontsize=9)
        ax.set_title(f'{k.upper()} band, every cycle stacked',color=BC[k],fontweight='bold'); ax.set_xlabel('x (px)')
        if j==0: ax.set_ylabel('cycle N (how far the light has gone)')
        else: ax.set_yticklabels([])
        cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=0.035,pad=0.08); cb.set_label('amount of band (1 = evenly spread)')
    fig.suptitle(f'Figure 2{"a" if m=="annealed" else "b"}.  From the picture to branched flow:  {NAME[m]}',fontsize=16,fontweight='bold',y=0.985)
    fig.text(0.06,0.955,'Top: your picture. Middle: ignore up/down and measure how much of each band sits above each slot of x (a spike = the band stands vertical there).\nBottom: do that for every cycle and stack the rows, cycle 0 at the bottom. Bright streaks = branches. Dashed lines = the cycles shown on top.',fontsize=12.5,va='top')
    colour_key(fig,[0.30,0.645,0.4,0.008])
    fig.savefig(f'F2{"a" if m=="annealed" else "b"}_bands_to_branches_{m}.png',dpi=72)
def fig3():
    modes=['annealed','quenched']; N=13
    fig=plt.figure(figsize=(17,25))
    g=fig.add_gridspec(4,2,height_ratios=[1,1,1.05,1.05],left=0.07,right=0.98,top=0.92,bottom=0.05,hspace=0.25,wspace=0.08)
    for j,m in enumerate(modes):
        show(fig.add_subplot(g[0,j]),m,N,title=f'{NAME[m]}\nthe picture at cycle {N}')
        show(fig.add_subplot(g[1,j]),m,N,fade=True,title='same picture, everything with S > 2 greyed out:\nonly the low-stretch material is left')
        tot=sum(D[f'{m}_{k}_rho'] for k in BANDS)/3; low=sum(D[f'{m}_{k}_low'] for k in BANDS)/3
        for i,(r,lab) in enumerate([(tot,'all the band material'),(low,'only the low-stretch band material (S ≤ 2)')]):
            ax=fig.add_subplot(g[2+i,j]); r=r.copy(); r[r<=0]=np.nan
            im=ax.imshow(r,origin='lower',aspect='auto',extent=[0,R,-0.5,60.5],cmap='magma',norm=LogNorm(0.05,50),interpolation='nearest')
            ax.set_facecolor('k'); ax.axhline(N,color='w',ls='--',lw=.8)
            ax.set_title(f'branched flow of the three bands: {lab}',fontsize=11.5); ax.set_xlabel('x (px)')
            if j==0: ax.set_ylabel('cycle N')
            cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=0.035,pad=0.16); cb.set_label('amount of band material (1 = evenly spread, black = none)')
    fig.suptitle('Figure 3.  Are the branches made of low-stretch material?',fontsize=17,fontweight='bold',y=0.985)
    fig.text(0.07,0.955,'Rows 1-2: the picture, then only its low-stretch part. Rows 3-4: branched flow of the bands, first with all their material, then with only the low-stretch part.\nIf the bottom picture still shows the branches, they are made of low-stretch material. If it goes black, they are not.',fontsize=12.5,va='top')
    fig.savefig('F3_lowstretch_branches.png',dpi=72)
def fig4():
    modes=['annealed','quenched','cosine']
    labs=['S ≤ 2 (low-stretch)','2 - 5','5 - 10','10 - 100','> 100']
    cols=[CMAP(0.0)[:3],CMAP(np.log(3.2)/np.log(15))[:3],CMAP(np.log(7)/np.log(15))[:3],CMAP(0.93)[:3],CMAP(1.0)[:3]]
    cols=[tuple(0.55*np.array(c)+0.45*np.array([1,1,1])) if i>0 else (0.0,0.45,0.55) for i,c in enumerate(cols)]
    fig=plt.figure(figsize=(19,13.5))
    g=fig.add_gridspec(2,3,height_ratios=[1,0.95],left=0.06,right=0.98,top=0.88,bottom=0.12,hspace=0.18,wspace=0.10)
    for j,m in enumerate(modes):
        show(fig.add_subplot(g[0,j]),m,13,title=f'{NAME[m]}\nthe picture at cycle 13')
        ax=fig.add_subplot(g[1,j]); sh=D[f'{m}_share']
        ax.stackplot(np.arange(61),100*sh.T,colors=cols,labels=labs,edgecolor='w',linewidth=0.3)
        ax.axvline(13,color='k',ls='--',lw=1); ax.set_xlim(0,60); ax.set_ylim(0,100); ax.set_xlabel('cycle N')
        if j==0: ax.set_ylabel('% of the picture'); 
        else: ax.set_yticklabels([])
        ax.set_title(f'low-stretch share: cycle 13 {100*sh[13,0]:.1f}%,  cycle 60 {100*sh[60,0]:.2f}%',fontsize=11)
    h,l=ax.get_legend_handles_labels(); fig.legend(h[::-1],[f'S {x}' if not x.startswith('S') else x for x in l[::-1]],loc='lower center',ncol=5,fontsize=12,frameon=False)
    fig.suptitle('Figure 4.  How much of the picture is stretched by how much, cycle by cycle',fontsize=17,fontweight='bold',y=0.975)
    fig.text(0.06,0.925,'Top: the picture at cycle 13. Bottom: at every cycle, the share of the picture in each stretch class (darkest = low-stretch, palest = most stretched). The dashed line is cycle 13.',fontsize=12.5)
    fig.savefig('F4_stretch_shares.png',dpi=72)
def fig5():
    modes=['annealed','quenched','cosine']; C={'annealed':'#1f77b4','quenched':'#ff7f0e','cosine':'#2ca02c'}
    fig,ax=plt.subplots(figsize=(12,8)); N=np.arange(61)
    for m in modes:
        f=D[f'{m}_share'][:,0]; ax.semilogy(N[1:],100*f[1:],'o-',ms=3,color=C[m],label=NAME[m])
    sel=(N>=20)
    a=np.polyfit(N[sel],D['annealed_mean'][sel],1)[0]; v=np.polyfit(N[sel],D['annealed_var'][sel],1)[0]
    pr=a*a/(2*v); f=D['annealed_share'][:,0]; meas=-np.polyfit(N[sel],np.log(f[sel]),1)[0]
    xx=np.arange(20,61); ax.semilogy(xx,100*f[20]*np.exp(-pr*(xx-20)),'k--',lw=1.5,label=f'log-normal prediction for annealed (falls {100*pr:.1f}% per cycle; measured {100*meas:.1f}%)')
    ax.set_xlabel('cycle N'); ax.set_ylabel('% of the picture that is low-stretch (S ≤ 2)')
    ax.set_title('Figure 5.  How fast the low-stretch part of the picture disappears',fontsize=15,fontweight='bold',loc='left')
    ax.grid(alpha=.3,which='both'); ax.legend(fontsize=11,loc='lower left'); ax.set_xlim(0,60)
    fig.tight_layout(); fig.savefig('F5_lowstretch_decay.png',dpi=90)
    print('annealed alpha v pred meas',a,v,pr,meas)
    for m in modes: print(m,[round(100*D[f'{m}_share'][n,0],3) for n in (5,13,20,40,60)])
for f in sys.argv[1:]: eval(f)
