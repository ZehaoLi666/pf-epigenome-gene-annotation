import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
mpl.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.titlesize':7,'axes.labelsize':6.5,'xtick.labelsize':5.5,'ytick.labelsize':5.5,
    'axes.linewidth':0.6,'xtick.major.width':0.5,'ytick.major.width':0.5,'xtick.major.size':2,'ytick.major.size':2,
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
INK='#222222'
COL=dict(agree='#4a78c2', ol_unit='#a6400a', ol_lnc='#e08214', ol_anti='#f5b041', ol_none='#fbd9a6', short='#7b4fa0', noutr='#2a9d8f',
         between='#c4a46b', na='#e3e3e3', ref='#1f3b73', cds='#1f3b73', utr='#8fb3e0', lnc='#e08214', cage='#c23b6e', daft='#4a78c2',
         tss='#7b4fa0', tts='#2a9d8f', body='#4a78c2', inter='#b5b5b5', trunc='#9a9a9a', art='#d9822b')
def plabel(ax,l,dx=-0.02,dy=1.02):
    p=ax.get_position(); ax.figure.text(max(0.012,p.x0+dx*p.width),p.y0+dy*p.height,l,fontsize=9,fontweight='bold',va='bottom',ha='right')
def title(ax,t,y=1.03): ax.text(0,y,t,transform=ax.transAxes,fontsize=7,va='bottom',ha='left')
def box(ax,x,y,w,h,fc='white',ec=INK,lw=0.8,ls='-'): ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0,rounding_size=1.0",fc=fc,ec=ec,lw=lw,ls=ls))
def arrow(ax,a,b,c=INK): ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=7,color=c,lw=0.8))
def overlaps(fig):
    r_=fig.canvas.get_renderer(); tx=[(t,t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    return [(a.get_text()[:25],b.get_text()[:25]) for i,(a,ba) in enumerate(tx) for b,bb in tx[i+1:] if ba.overlaps(bb)]
