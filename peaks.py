import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from core import *
plt.rcParams.update({'font.size':12})
D=np.load('stats.npz'); r=D['annealed_red_rho'][3]; xc=np.arange(480)*4+2
pic,S,orow=picture_and_S('annealed',3,1); lab=band_of(orow)
im=overlay(pic,S)
fade=np.where((lab==1)[...,None],im,0.75*im+0.25)   # red band stays vivid, rest faded
fig=plt.figure(figsize=(15,12.5))
g=fig.add_gridspec(2,1,height_ratios=[1.45,0.75],hspace=0.12,left=0.08,right=0.98,top=0.90,bottom=0.07)
a=fig.add_subplot(g[0]); a.imshow(fade[600:1300],extent=[0,R,1300,600],interpolation='nearest',aspect='auto')
a.set_ylabel('row of the picture\n(= direction the light travels)'); plt.setp(a.get_xticklabels(),visible=False)
cols=[(162,'peak: band stands\nUPRIGHT here','w'),(718,'peak','w'),(1534,'peak','w'),(1134,'dip: band lies\nFLAT and thin','k')]
for x,t,c in cols:
    a.axvspan(x-12,x+12,color='cyan' if 'dip' not in t else 'orange',alpha=0.35)
    a.text(x,620,t,ha='center',va='top',fontsize=11,fontweight='bold',color='k',bbox=dict(fc='white',ec='none',alpha=.8))
a.set_title('Cycle 3, annealed: your picture. Only the RED band is being counted. Shaded strips = the x-slots of the peaks and the dip.',loc='left',fontsize=12.5)
b=fig.add_subplot(g[1],sharex=a); b.semilogy(xc,r,color='#d62728',lw=1.7); b.axhline(1,color='k',ls=':')
for x,t,c in cols: b.axvspan(x-12,x+12,color='cyan' if 'dip' not in t else 'orange',alpha=0.35)
b.set_ylim(0.05,60); b.set_xlim(0,R); b.set_xlabel('x (px)'); b.set_ylabel('amount of red band\n(1 = evenly spread)')
b.text(162,6.5,'4x',ha='center',fontweight='bold'); b.text(718,6.5,'4x',ha='center',fontweight='bold'); b.text(1534,6.5,'4x',ha='center',fontweight='bold'); b.text(1134,0.09,'0.2x',ha='center',fontweight='bold')
fig.suptitle('What a peak is: count how much red band is stacked up above each slot of x',fontsize=16,fontweight='bold')
fig.savefig('peaks_explained.png',dpi=80)
