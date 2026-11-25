"""Shared machinery: the user's starting picture, the kick-drift map, the
picture at cycle N (exact pullback), and S at the current location."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, kd_engine as E
from matplotlib.colors import ListedColormap
R=1920
def params(mode): return E.Params(cycles=60,mode=mode,seed=0,K_STD=0.5)
_K={}
def kicks(mode):
    if mode not in _K: _K[mode]=E.make_kicks(params(mode))
    return _K[mode]
P=params('annealed'); Tau=P.Tau; row0=P.row0
# ---- the starting picture: white, 42x42 black dots, red/green/blue bands ----
u=R/41.0
def row_of(yd): return (41-yd)*u
BANDS={'red':(row_of(23.5),row_of(22.5)),'green':(row_of(20.5),row_of(19.5)),'blue':(row_of(17.5),row_of(16.5))}
BANDRGB={'red':(1,0,0),'green':(0,1,0),'blue':(0,0,1)}
def start_image():
    img=np.ones((R,R,3))
    yy,xx=np.mgrid[0:R,0:R].astype(float)
    for k,(a,b) in BANDS.items():
        img[int(round(a)):int(round(b))]=BANDRGB[k]
    rad=4.4
    for i in range(42):
        for j in range(42):
            cx=i*u; cy=row_of(j)
            x0,x1=int(max(cx-6,0)),int(min(cx+7,R)); y0,y1=int(max(cy-6,0)),int(min(cy+7,R))
            if x1<=x0 or y1<=y0: continue
            m=(xx[y0:y1,x0:x1]-cx)**2+(yy[y0:y1,x0:x1]-cy)**2<=rad**2
            img[y0:y1,x0:x1][m]=0
    return img
IMG0=start_image()
def inv_step(x,y,dp):
    x=x-Tau*(y-row0); xf=np.floor(x); i0=np.mod(xf.astype(np.int64),R); d0=dp[i0]
    return x, y-(d0+(dp[(i0+1)%R]-d0)*(x-xf))
def picture_and_S(mode,N,step=1):
    """picture at cycle N (each pixel shows the colour of the material now there)
    and S of that material, on the pixel grid (step px)."""
    ys,xs=np.mgrid[0:R:step,0:R:step].astype(float)
    x=xs.ravel(); y=ys.ravel(); kk=kicks(mode)
    for k in range(N-1,-1,-1): x,y=inv_step(x,y,kk[k])
    oi=np.mod(np.rint(y).astype(int),R); oj=np.mod(np.rint(x).astype(int),R)
    pic=IMG0[oi,oj].reshape(ys.shape+(3,))
    st=E.new_state(x,y)
    for k in range(N): E.step(st,kk[k],Tau,row0)
    S=E.sigma_max(st['a'],st['b'],st['c'],st['d']).reshape(ys.shape)
    origin_row=y.reshape(ys.shape)
    return pic,S,origin_row
t=np.linspace(0,1,256)[:,None]
CMAP=ListedColormap(np.hstack([0.0*(1-t)**2+0.129*2*t*(1-t)+0.993*t**2,0.25*(1-t)**2+0.566*2*t*(1-t)+0.906*t**2,0.35*(1-t)**2+0.692*2*t*(1-t)+0.144*t**2]))
def overlay(pic,S,alpha=0.28,smax=15):
    tt=np.clip(np.log(np.maximum(S,1))/np.log(smax),0,1)
    col=CMAP(tt)[...,:3]
    return (1-alpha)*pic+alpha*col
def band_of(origin_row):
    lab=np.zeros(origin_row.shape,int)
    for n,(k,(a,b)) in enumerate(BANDS.items(),1): lab[(origin_row>=a)&(origin_row<b)]=n
    return lab
