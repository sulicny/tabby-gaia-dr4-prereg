#!/usr/bin/env python
"""Independent offline DR4-format generator; orbital functions copied verbatim from sim.py."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
G=6.67430e-11; MS=1.98847e30; MJ=1.89813e27; AU=1.495978707e11; PLX=2.254545
WINDOWS=[(1029.7,1031.5),(1034.5,1035.3),(1301,1303.6),(1345.640,1345.656),(1367.3,1368.3)]
TABBY=2081900940499099136
def alpha(K,P,e,si): return K*P*86400*np.sqrt(1-e*e)/(2*np.pi*AU)*PLX/si
# Exact two-body RV relation (star mass fixed).
def kval(mp,P,e,si):
 q=mp*MJ; return (2*np.pi*G/(P*86400))**(1/3)*q*si/(1.19*MS+q)**(2/3)/np.sqrt(1-e*e)
def mass_from_a(a,P):
 low=np.zeros_like(np.asarray(a)); high=np.full_like(low,1000.)
 for _ in range(45):
  mid=(low+high)/2; aa=alpha(kval(mid,P,.09,1),P,.09,1); low=np.where(aa<a,mid,low); high=np.where(aa>=a,mid,high)
 return (low+high)/2

def orbit(P,e,w,ci,node,times,opposite=False):
 wp=w+(0 if opposite else np.pi); nu=np.pi/2-wp
 E0=2*np.arctan2(np.sqrt(1-e)*np.sin(nu/2),np.sqrt(1+e)*np.cos(nu/2)); M0=E0-e*np.sin(E0)
 M=M0[:,None]+2*np.pi*(times[None,:]-2458730.164)/P[:,None]; E=M.copy()
 for _ in range(10): E-=(E-e[:,None]*np.sin(E)-M)/(1-e[:,None]*np.cos(E))
 X=np.cos(E)-e[:,None]; Y=np.sqrt(1-e[:,None]**2)*np.sin(E)
 A=np.cos(node)*np.cos(w)-np.sin(node)*np.sin(w)*ci; B=np.sin(node)*np.cos(w)+np.cos(node)*np.sin(w)*ci
 F=-np.cos(node)*np.sin(w)-np.sin(node)*np.cos(w)*ci; H=-np.sin(node)*np.sin(w)+np.cos(node)*np.cos(w)*ci
 return (A[:,None]*X+F[:,None]*Y),(B[:,None]*X+H[:,None]*Y)

def make_dataset(outdir,case='injected',seed=1,sigma_scale=.75,fmt='gaia_array',window=4,mass_mjup=9.4,n_controls=300):
 if case not in {'injected','null'} or fmt not in {'gaia_array','standard'}: raise ValueError('invalid case/format')
 out=Path(outdir); out.mkdir(parents=True,exist_ok=True); (out/'SYNTHETIC_MARKER').write_text('Synthetic; never real DR4 data.\n')
 rng=np.random.default_rng(seed)
 g=pd.read_csv(ROOT/'data/gost_tabby_dr4window.csv',skipinitialspace=True)
 jd=g['ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]'].to_numpy(); psi=g['scanAngle[rad]'].to_numpy(); pf=g.parallaxFactorAlongScan.to_numpy()
 raw=pd.read_csv(ROOT/'data/controls_raw.csv'); tab=pd.read_csv(ROOT/'data/tabby_dr3.csv').iloc[0]
 kt=float(tab.parallax_error)*np.sqrt(float(tab.astrometric_n_good_obs_al))
 k=raw.parallax_error.to_numpy()*np.sqrt(raw.astrometric_n_good_obs_al.to_numpy()); good=np.isfinite(k)&(k>0); raw=raw[good].reset_index(drop=True); k=k[good]
 truth={}
 def epochs(sid,scale,inject):
  keep=rng.random(len(jd))<.9; times=jd[keep]; ang=psi[keep]; factors=pf[keep]; n=len(times)
  cc=times[:,None]+np.arange(9)[None,:]*4.4/86400; angles=ang[:,None]+rng.normal(0,1e-7,(n,9))
  t=(cc-jd.mean())/365.25; model=2.254545*factors[:,None]+(.2+3*t)*np.sin(angles)+(-.1-2*t)*np.cos(angles)
  if inject:
   lo,hi=WINDOWS[window-1]; p=rng.uniform(lo,hi); e=rng.normal(.09,.05)
   while not 0<=e<.5: e=rng.normal(.09,.05)
   z=rng.normal(); w=np.deg2rad(141+z*(89 if z<0 else 84)); z=rng.normal(); b=.76+z*(.04 if z<0 else .05); z=rng.normal(); rho=max(1,485+z*(96 if z<0 else 100))
   ar=(G*rho*(p*86400)**2/(3*np.pi))**(1/3); ci=np.clip(b/ar*(1+e*np.sin(w+np.pi))/(1-e*e),-.999,.999)*rng.choice([-1,1]); si=np.sqrt(1-ci*ci); node=rng.uniform(0,2*np.pi)
   dec,ra=orbit(*[np.array([v]) for v in (p,e,w,ci,node)],cc.ravel()); a=alpha(kval(mass_mjup,p,e,si),p,e,si)
   model+=a*(ra.reshape(n,9)*np.sin(angles)+dec.reshape(n,9)*np.cos(angles)); truth.update(P=p,e=e,omega=w,cos_i=ci,node=node,alpha=a)
  sig=sigma_scale*.0417*scale*3; corr=rng.uniform(-500,500,n)*1e9
  values=model+rng.normal(0,sig,(n,9)); used=np.ones((n,9),bool)
  for i in range(n):
   if rng.random()<.1: used[i,rng.integers(9)]=False
  rows=[]
  for i in range(n):
   arr=dict(obs_time_tcb=(cc[i]-2455197.5)*86400e9-corr[i],scan_pos_angle=np.rad2deg(angles[i]),centroid_pos_al=values[i],centroid_pos_error_al=np.full(9,.9*sig),used_by_agis_al=used[i])
   base=dict(source_id=sid,transit_id=i,parallax_factor_al=factors[i],agis_source_excess_noise=.436*sig,obs_time_bary_corr=corr[i],ra0=301.5643,dec0=44.4568)
   if fmt=='gaia_array': rows.append(dict(base,**{key:json.dumps(v.tolist()) for key,v in arr.items()}))
   else:
    for j in range(9): rows.append(dict(base,**{key:v[j] for key,v in arr.items()}))
  return pd.DataFrame(rows)
 epochs(TABBY,1,case=='injected').to_csv(out/'tabby_epoch_astrometry.csv',index=False)
 source=dict(source_id=TABBY,parallax=PLX+rng.normal(0,.003),parallax_error=.006,ruwe=.9,ipd_frac_multi_peak=0,ipd_gof_harmonic_amplitude=.001)
 pd.DataFrame([source]).to_csv(out/'tabby_dr4_source.csv',index=False)
 pd.DataFrame([dict(source_id=TABBY,has_epoch_astrometry=True)]).to_csv(out/'tabby_dr4_flags.csv',index=False)
 controls=[]
 for i in range(n_controls):
  j=int(rng.integers(len(raw))); sid=100000+i; row=raw.iloc[j].to_dict(); row['source_id']=sid; controls.append(row)
  epochs(sid,k[j]/kt,False).to_csv(out/f'ctl_{sid}.csv',index=False)
 pd.DataFrame(controls).to_csv(out/'controls_dr4.csv',index=False)
 (out/'simulation_truth.json').write_text(json.dumps(dict(case=case,seed=seed,fmt=fmt,**truth),indent=2))
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser(); p.add_argument('outdir'); p.add_argument('--case',choices=['injected','null'],default='injected'); p.add_argument('--seed',type=int,default=1); p.add_argument('--fmt',choices=['gaia_array','standard'],default='gaia_array'); p.add_argument('--n-controls',type=int,default=300)
 a=p.parse_args(); make_dataset(**vars(a))
