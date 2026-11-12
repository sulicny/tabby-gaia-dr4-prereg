#!/usr/bin/env python
"""Offline Gaia AL injection recovery; omega_p = omega_s + pi at transit.
Opposite-convention sensitivity uses omega_p = omega_s. Gaussian split errors
are represented by equal-probability half normals (quoted centre is median).
Independent dropout masks are drawn for all 1,500 simulation indices;
matched noise and masks are shared across truth families for comparisons.
All fits use the scenario's known nominal sigma, including empirical nulls.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MPLCONFIGDIR']='/tmp/tabby-matplotlib'
import csv,json,time
from pathlib import Path
import numpy as np
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
rng=np.random.default_rng(260720643); START=time.monotonic()
G=6.67430e-11; MS=1.98847e30; MJ=1.89813e27; AU=1.495978707e11; PLX=2.254545
N=1500; GROUPS=1500; BATCH=N//GROUPS; masses=np.array([2,3,5,7.5,9.4,12,15.])
windows=[(1029.7,1031.5),(1034.5,1035.3),(1301,1303.6),(1345.648,1345.648),(1367.3,1368.3)]
def read(name):
 with open(ROOT/'data'/name) as f: return list(csv.DictReader(f,skipinitialspace=True))
g=read('gost_tabby_dr4window.csv'); tab=read('tabby_dr3.csv')[0]; controls=read('controls_raw.csv')
jd=np.array([float(x['ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]']) for x in g]); th=np.array([float(x['scanAngle[rad]']) for x in g]); pf=np.array([float(x['parallaxFactorAlongScan']) for x in g]); s=np.sin(th); c=np.cos(th)
def design(mask):
 t=(jd[mask]-jd[mask].mean())/365.25
 return np.column_stack([s[mask],c[mask],pf[mask],s[mask]*t,c[mask]*t]),t
m3=np.array([x['ObservationTimeAtGaia[UTC]'][:10]<='2017-05-28' for x in g]); x3,_=design(m3); sigma3=float(tab['parallax_error'])/np.sqrt(np.linalg.inv(x3.T@x3)[2,2])
k=np.array([float(x['parallax_error'])*np.sqrt(float(x['astrometric_n_good_obs_al'])) for x in controls]); k=k[np.isfinite(k)&(k>0)]; kt=float(tab['parallax_error'])*np.sqrt(float(tab['astrometric_n_good_obs_al']))
cal={'sigma_FoV_DR3_eff_mas':sigma3,'DR3_FoV':int(m3.sum()),'DR3_good_CCD':int(tab['astrometric_n_good_obs_al']),'CCD_per_FoV':float(tab['astrometric_n_good_obs_al'])/m3.sum(),'control_count':len(k),'control_log_k_sd':float(np.std(np.log(k),ddof=1)),'tabby_control_percentile':float(100*np.mean(k<=kt))}
def split(mid,lo,hi,n):
 z=rng.normal(size=n); return mid+z*np.where(z<0,lo,hi)
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
def template(P,times,opp=False): return orbit(np.array([P]),np.array([.09]),np.deg2rad(np.array([141.])),np.array([0.]),np.array([0.]),times,opp)[0][0]
def truth(lo,hi,kind,opp):
 P=rng.uniform(lo,hi,N); e=np.full(N,.09); w=np.full(N,np.deg2rad(141.)); b=np.full(N,.76); rho=np.full(N,485.)
 if kind=='posterior':
  e=rng.normal(.09,.05,N)
  while np.any((e<0)|(e>=.5)):
   bad=(e<0)|(e>=.5); e[bad]=rng.normal(.09,.05,bad.sum())
  w=np.deg2rad(rng.normal(141,87,N)); b=split(.76,.04,.05,N); rho=np.maximum(split(485,96,100,N),1)
 ar=(G*rho*(P*86400)**2/(3*np.pi))**(1/3); ci=b/ar*(1+e*np.sin(w+(0 if opp else np.pi)))/(1-e*e); ci=np.clip(ci,-.999,.999)*rng.choice([-1,1],N); si=np.sqrt(1-ci*ci)
 K=np.full(N,162.)
 if kind=='posterior':
  K=split(162,69,76,N)
  while np.any(K<0):
   bad=K<0; K[bad]=split(162,69,76,bad.sum())
 elif isinstance(kind,float): K=kval(kind,P,e,si)
 a=alpha(K,P,e,si); dec,ra=orbit(P,e,w,ci,rng.uniform(0,2*np.pi,N),jd,opp)
 return a[:,None]*(ra*s+dec*c),a,float(np.max(abs(ci))),P
def basis(z,base):
 r=z-base@(base.T@z); u,sv,_=np.linalg.svd(r,full_matrices=False); return u[:,sv>1e-10]
def upper(coeff,cov):
 # Flat amplitude, uniform orientation: integrate likelihood in da dphi.
 vals=np.linalg.eigvalsh(cov); amax=np.linalg.norm(coeff,axis=1).max()+8*np.sqrt(vals.max()); aa=np.linspace(0,amax,240); phi=np.arange(64)*2*np.pi/64; dirs=np.column_stack([np.cos(phi),np.sin(phi)]); inv=np.linalg.inv(cov)
 v=dirs@inv@dirs.T; quad=np.diag(v); lin=coeff@inv@dirs.T
 log=-.5*aa[None,:,None]**2*quad[None,None,:]+aa[None,:,None]*lin[:,None,:]; log-=log.max(axis=(1,2),keepdims=True)
 prob=np.exp(log).mean(axis=2); prob[:,0]*=.5; prob[:,-1]*=.5; cum=np.cumsum(prob,axis=1); idx=(cum<.95*cum[:,-1,None]).sum(axis=1)
 return aa[np.minimum(idx,len(aa)-1)]
if __name__ == '__main__' and '--self-test' in __import__('sys').argv:
 from scipy.special import erfinv
 measured=upper(np.zeros((1,2)),np.eye(2))[0]
 expected=np.sqrt(2)*erfinv(.95)
 assert abs(measured/expected-1)<.02, (measured,expected)
 for opp in [False,True]:
  d,_=orbit(np.array([1345.648]),np.array([.09]),np.deg2rad(np.array([141.])),np.array([0.]),np.array([0.]),np.array([2458730.164]),opp)
  assert abs(d[0,0])<1e-10
 mass=np.array([2.,9.4,15.]); aa=alpha(kval(mass,1345.648,.09,1),1345.648,.09,1)
 assert np.allclose(mass_from_a(aa,1345.648),mass,rtol=1e-8)
 print('PASS: analytic flat-amplitude upper limit (<2% quadrature error), both transit conventions, exact mass inversion.')
 raise SystemExit(0)
periods=1/np.linspace(1/2500,1/500,int(np.ceil(2010*(1/500-1/2500)/.05))+1)
scenarios=[(name,r,floor,keep,opp,emp) for name,r,floor,keep,opp,emp in [('pessimistic',1,0,.9,False,False),('central',.75,0,.9,False,False),('optimistic',.5,0,.9,False,False),('pessimistic_floor',1,.03,.9,False,False),('central_floor',.75,.03,.9,False,False),('optimistic_floor',.5,.03,.9,False,False),('empirical',.75,0,.9,False,True),('dropout_80',.75,0,.8,False,False),('opposite_phase',.75,0,.9,True,False)]]
records=[]; example=None
for wi,(lo,hi) in enumerate(windows,1):
 P0=(lo+hi)/2; nominal=alpha(162,P0,.09,1)
 for name,r,floor,keep,opp,emp in scenarios:
  sig=np.hypot(r*sigma3,floor); kinds=['null','posterior','face_value']+list(masses)+['blend_%s_%s'%(sw,e) for sw in ('0.3','0.5','0.8') for e in ('1e-4','1e-3','1e-2','1e-1')]; ns=len(kinds)
  signals=[]; amps=[]; true_periods=[]; maxci=0
  for kind in kinds:
   if kind in ['null'] or isinstance(kind,str) and kind.startswith('blend'):
    y=np.zeros((N,len(jd))); a=np.zeros(N); truth_P=np.full(N,P0)
    if str(kind).startswith('blend'):
     sw=float(kind.split('_')[1]); eps=float(kind.split('_')[2]); proj=np.cos(th-np.deg2rad(95.967)); y[:]=eps*10**(-.4*(17.603-11.759))*1949.5*proj*np.exp(-.5*(1.9495*proj/sw)**2)
   else: y,a,mc,truth_P=truth(lo,hi,kind,opp); maxci=max(maxci,mc)
   signals.append(y); amps.append(a); true_periods.append(truth_P)
  stats={key:np.zeros((ns,N)) for key in ['chi0','dacc','dc','dfree','pbest','afit','ratio','ul','free_penalty','two_penalty']}; static_error=0
  for group in range(GROUPS):
   ix=slice(group*BATCH,(group+1)*BATCH); mask=rng.random(len(jd))<keep; x,t=design(mask); n=mask.sum(); q=np.linalg.qr(x)[0]; tc=template(P0,jd[mask],False); z=np.column_stack([s[mask]*tc,c[mask]*tc]); zr=z-q@(q.T@z); cov=np.linalg.inv(zr.T@zr)*sig**2; pin=np.linalg.pinv(zr); qc=basis(z,q); qa=basis(np.column_stack([s[mask]*t*t/2,c[mask]*t*t/2]),q)
   y=np.stack([v[ix][:,mask] for v in signals]); scale=np.ones((ns,BATCH))
   if emp: scale[:]=rng.choice(k/np.median(k),(ns,BATCH))
   pass
   noise=rng.normal(size=y.shape)*sig*scale[:,:,None]; noise[-3:]=noise[0]; y+=noise
   flat=y.reshape(-1,n).T; residual=flat-q@(q.T@flat); ch0=np.sum(residual**2,axis=0)/sig**2; dc=np.sum((qc.T@residual)**2,axis=0)/sig**2; da=np.sum((qa.T@residual)**2,axis=0)/sig**2; coeff=(pin@residual).T; af=np.linalg.norm(coeff,axis=1)
   best=np.full(len(ch0),-np.inf); bp=np.zeros(len(ch0))
   for pp in periods:
    ph=2*np.pi*(jd[mask]-jd[mask].mean())/pp; zf=np.column_stack([s[mask]*np.cos(ph),s[mask]*np.sin(ph),c[mask]*np.cos(ph),c[mask]*np.sin(ph)]); qf=basis(zf,q); dd=np.sum((qf.T@residual)**2,axis=0)/sig**2; hit=dd>best; best[hit]=dd[hit]; bp[hit]=pp
   stats['free_penalty'][:,ix]=5*np.log(n); stats['two_penalty'][:,ix]=2*np.log(n)
   for key,v in [('chi0',ch0),('dc',dc),('dacc',da),('dfree',best),('pbest',bp),('afit',af)]: stats[key][:,ix]=v.reshape(ns,BATCH)
   for j in [0,1]: stats['ul'][j,ix]=upper(coeff.reshape(ns,BATCH,2)[j],cov)
   # A static projected separation must lie exactly in the position columns.
   static=1949.5*np.cos(th[mask]-np.deg2rad(95.967)); static_error=max(static_error,float(np.max(abs(static-q@(q.T@static)))))
   if wi==4 and name=='central' and group==0:
    j=2*BATCH; fit=zr@coeff[j]; example=(jd[mask],residual[:,j],fit,sig)
  from scipy.stats import chi2 as _c2
  thresholds={key:{str(p):(float(_c2.isf(p,2)) if key in ('dc','dacc') else float(np.quantile(stats[key][0],1-p))) for p in [.0027,.001]} for key in ['dc','dfree','dacc']}
  out={'window':'W'+str(wi),'period_days':P0,'scenario':name,'sigma_mas':sig,'keep_probability':keep,'N_sim':N,'alpha_nominal_mas':nominal,'max_abs_cos_i':maxci,'thresholds':thresholds,'null_distributions':{key:stats[key][0].tolist() for key in ['dc','dfree','dacc']},'static_blend_absorption_max_mas':static_error,'truths':{}}
  for j,kind in enumerate(kinds):
   det=stats['dc'][j]>thresholds['dc']['0.001']; free=(stats['dfree'][j]>thresholds['dfree']['0.001'])&(abs(stats['pbest'][j]/true_periods[j]-1)<.03)&(stats['dfree'][j]-stats['free_penalty'][j]>10)
   a=amps[j]; valid=a>0
   d={'true_period_days':true_periods[j].tolist(),'fit_distributions':{key:stats[key][j].tolist() for key in ['chi0','dc','dacc','dfree','pbest','afit','two_penalty','free_penalty']},'P_detect_Mc':float(det.mean()),'P_detect_free_recover_BIC10':float(free.mean()),'median_delta_chi2_Mc':float(np.median(stats['dc'][j])),'median_amplitude_mas':float(np.median(stats['afit'][j])),'median_recovered_true_amplitude_ratio':float(np.median(stats['afit'][j][valid]/a[valid])) if valid.any() else None,'median_chi2_M0':float(np.median(stats['chi0'][j])),'median_chi2_Mc':float(np.median(stats['chi0'][j]-stats['dc'][j])),'median_delta_chi2_M1':float(np.median(stats['dacc'][j])),'median_delta_chi2_Mfree':float(np.median(stats['dfree'][j])),'median_BIC_improvement_Mc':float(np.median(stats['dc'][j]-stats['two_penalty'][j])),'median_BIC_improvement_M1':float(np.median(stats['dacc'][j]-stats['two_penalty'][j])),'median_BIC_improvement_Mfree':float(np.median(stats['dfree'][j]-stats['free_penalty'][j])),'median_recovered_period_days':float(np.median(stats['pbest'][j]))}
   if j<2:
    ul=mass_from_a(stats['ul'][j],P0); detect=det&(stats['afit'][j]/nominal>.4)&(stats['afit'][j]/nominal<2); exclude=(ul<5.3)&~detect
    d.update(amplitude_UL95_mas_quantiles=np.quantile(stats['ul'][j],[.16,.5,.84]).tolist(),mass_UL95_MJup_quantiles=np.quantile(ul,[.16,.5,.84]).tolist(),claim_outcomes={'detect':float(detect.mean()),'exclude':float(exclude.mean()),'inconclusive':float((~(detect|exclude)).mean())})
   out['truths'][str(kind)]=d
  power=np.maximum.accumulate([out['truths'][str(v)]['P_detect_Mc'] for v in masses]); out['mass_detection_MJup']={str(level):float(np.interp(level,power,masses)) if power[0]<=level<=power[-1] else None for level in [.5,.9]}
  records.append(out)
 print('W%d complete; elapsed %.1fs'%(wi,time.monotonic()-START),flush=True)
result={'seed':260720643,'calibration':cal,'chi2_2_thresholds':{str(p):float(chi2.isf(p,2)) for p in [.0027,.001]},'period_grid_days':periods.tolist(),'method_notes':[__doc__,'UL numerical quadrature: 240 amplitude points, 64 orientations; tail extends 8 covariance sigmas beyond maximum fitted amplitude.','Exclude is evaluated only if not detected, making claim outcomes exclusive.','Mass-power interpolation uses monotone envelope; null means crossing outside supplied grid.','BJD TCB input used directly with supplied BJD transit anchor; no network/time-standard conversion.'], 'cases':records,'runtime_seconds':time.monotonic()-START}
with open(OUT/'forecast.json','w') as f: json.dump(result,f,indent=2,allow_nan=False)
rows=[]
for rec in records:
 for kind,d in rec['truths'].items():
  row=dict(window=rec['window'],scenario=rec['scenario'],truth=kind,sigma_mas=rec['sigma_mas'],**{k:v for k,v in d.items() if not isinstance(v,(dict,list))})
  for model,threshold in rec['thresholds'].items():
   for p,value in threshold.items(): row['threshold_'+model+'_'+p]=value
  for field in ['amplitude_UL95_mas_quantiles','mass_UL95_MJup_quantiles']:
   for label,value in zip(['p16','median','p84'],d.get(field,[None]*3)): row[field+'_'+label]=value
  for outcome in ['detect','exclude','inconclusive']: row['claim_'+outcome]=d.get('claim_outcomes',{}).get(outcome)
  for level,value in rec['mass_detection_MJup'].items(): row['mass_power_'+level+'_MJup']=value
  rows.append(row)
with open(OUT/'forecast.csv','w') as f:
 writer=csv.DictWriter(f,fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
for rec in records:
 if rec['window'] in ['W3','W4'] and rec['scenario'] in ['pessimistic','central','optimistic']:
  ax[0].plot(masses,[rec['truths'][str(v)]['P_detect_Mc'] for v in masses],ls='-' if rec['window']=='W4' else '--',color={'pessimistic':'C0','central':'C1','optimistic':'C2'}[rec['scenario']],label=rec['window']+' '+rec['scenario'])
ax[0].axvline(9.4,color='gray',ls=':'); ax[0].set(xlabel='Injected companion mass (MJup)',ylabel='Detection probability',ylim=(0,1)); ax[0].legend(fontsize=8)
times,res,fit,sig=example; ax[1].errorbar((times-times.mean())/365.25,res,yerr=sig,fmt='.',alpha=.65,label='AL residual after M0'); ax[1].plot((times-times.mean())/365.25,fit,'.',label='Mc fitted residual'); ax[1].set(xlabel='Years from mean kept epoch',ylabel='AL displacement (mas)'); ax[1].legend(fontsize=8)
central=next(v for v in records if v['window']=='W4' and v['scenario']=='central'); fig.suptitle('W4 central: nominal K power %.1f%%; null median mass UL %.1f MJup'%(100*central['truths']['face_value']['P_detect_Mc'],central['truths']['null']['mass_UL95_MJup_quantiles'][1])); fig.tight_layout(); fig.savefig(OUT/'forecast.png',dpi=170); plt.close(fig)
print('Calibration:',json.dumps(cal)); print('Window scenario       posterior power  nominal power   null mass UL95 (MJup)')
for rec in records: print('%s %-18s %6.3f %6.3f %7.2f'%(rec['window'],rec['scenario'],rec['truths']['posterior']['P_detect_Mc'],rec['truths']['face_value']['P_detect_Mc'],rec['truths']['null']['mass_UL95_MJup_quantiles'][1]))
print('Runtime %.1fs; outputs forecast.json, forecast.csv, forecast.png'%(time.monotonic()-START))
