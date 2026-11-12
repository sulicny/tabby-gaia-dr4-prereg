#!/usr/bin/env python
"""Offline frozen-preregistration Gaia DR4 analysis."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MPLCONFIGDIR']='/tmp/tabby-phaseb-mpl'
import argparse,csv,hashlib,json,tempfile,re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
from gen_sim_dr4 import make_dataset
ROOT=Path(__file__).resolve().parents[1]
BARY_SIGN=1  # VERIFY on Dec 2: barycentric TCB = Gaia TCB + correction.
PERIODS=np.array([1030.6,1034.9,1302.3,1345.648,1367.8])
NOMINAL=np.array([.0345,.0346,.0435,.0450,.0457])
TABBY='2081900940499099136'

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hash_status(path):
 if not Path(path).is_file(): return False
 for line in Path(path).read_text().splitlines():
  fields=line.split()
  if len(fields)>1 and Path(fields[1].lstrip('*')).name=='prereg.md': return fields[0]==sha(ROOT/'prereg.md')
 return False

def fitter_hash_status(path):
 if Path(path).is_file():
  for line in Path(path).read_text().splitlines():
   fields=line.split()
   if len(fields)>1 and fields[1].lstrip('*')=='phaseB/02_fit.py': return fields[0]==sha(__file__)
 return None

class DataValidationError(RuntimeError): pass

def cell(value):
 text=str(value).strip()
 text=text.strip('[]{}()').strip()
 return re.split(r'[,\s]+',text.lower()) if text else ['']
def number(x):
 try: return float(x)
 except (TypeError,ValueError): return np.nan

def collapse(path,sign=BARY_SIGN,real=False,source_id=None):
 groups={}; total=failed=0
 keys=['obs_time_tcb','scan_pos_angle','centroid_pos_al','centroid_pos_error_al','used_by_agis_al','parallax_factor_al','agis_source_excess_noise','obs_time_bary_corr']
 with open(path,newline='') as f:
  for row in csv.DictReader(f):
   if source_id is not None and str(row.get('source_id',''))!=str(source_id): continue
   arrays={k:cell(row[k]) for k in keys}; n=max(map(len,arrays.values()))
   if any(len(v) not in (1,n) for v in arrays.values()):
    total+=n; failed+=n; continue
   for j in range(n):
    v={k:a[j if len(a)>1 else 0] for k,a in arrays.items()}
    total+=1
    vals=np.array([number(v[k]) for k in keys if k!='used_by_agis_al'])
    valid_bool=v['used_by_agis_al'] in ('true','t','1','1.0','false','f','0','0.0')
    time,angle,y,err,pf,excess,corr=vals
    jd=2455197.5+(time+sign*corr)/86400e9
    if real and np.isfinite(time) and np.isfinite(corr) and not (2456860<=jd<=2458870 and abs(corr/1e9)<=600):
     raise DataValidationError(f'ABORT: {path}: barycentric JD must be within 2456860..2458870 and |obs_time_bary_corr| <= 600 s (input ns); got JD={jd}, correction={corr/1e9} s')
    if not np.isfinite(vals).all() or not valid_bool:
     failed+=1; continue
    if v['used_by_agis_al'] not in ('true','t','1','1.0'): continue
    if err<0 or excess<0 or err*err+excess*excess<=0: continue
    groups.setdefault(row['transit_id'],[]).append([2455197.5+(time+sign*corr)/86400e9,np.deg2rad(angle),y,pf,1/(err*err+excess*excess)])
 if real and total and failed/total>.10:
  raise DataValidationError(f'ABORT: {path}: {failed}/{total} samples failed to parse (>10%)')
 result=[]
 for rows in groups.values():
  a=np.array(rows); w=a[:,4]; sw=w.sum(); mean=lambda x: float(w@x/sw)
  result.append([mean(a[:,0]),np.arctan2(mean(np.sin(a[:,1])),mean(np.cos(a[:,1]))),mean(a[:,2]),mean(a[:,3]),sw**-.5])
 result=np.array(result)
 if len(result)<=7: raise ValueError('insufficient valid transits for orbital fit')
 return result[np.argsort(result[:,0])]

def template(p,jd,opposite=False):
 e=.09; w=np.deg2rad(141); nu=np.pi/2-w-(0 if opposite else np.pi)
 E0=2*np.arctan2(np.sqrt(1-e)*np.sin(nu/2),np.sqrt(1+e)*np.cos(nu/2)); m=E0-e*np.sin(E0)+2*np.pi*(jd-2458730.164)/p; E=m.copy()
 for _ in range(12): E-=(E-e*np.sin(E)-m)/(1-e*np.cos(E))
 return np.cos(w)*(np.cos(E)-e)-np.sin(w)*np.sqrt(1-e*e)*np.sin(E)

def amplitude(m,p,plx):
 G=6.67430e-11; MS=1.98847e30; MJ=1.89813e27; AU=1.495978707e11
 q=m*MJ; K=(2*np.pi*G/(p*86400))**(1/3)*q/(1.19*MS+q)**(2/3)/np.sqrt(1-.09**2)
 return K*p*86400*np.sqrt(1-.09**2)/(2*np.pi*AU)*plx

def mass(a,p,plx):
 if not np.isfinite(plx) or plx<=0: return float('inf')
 high=1000.
 while amplitude(high,p,plx)<a: high*=2
 return brentq(lambda m:amplitude(m,p,plx)-a,0,high)

def upper(c,cov):
 # Flat da dphi prior; no radial Jacobian (which would impose a different prior).
 radius=np.linspace(0,np.linalg.norm(c)+9*np.sqrt(np.linalg.eigvalsh(cov).max()),512)
 phi=np.arange(128)*2*np.pi/128; dirs=np.array([np.cos(phi),np.sin(phi)]).T
 d=radius[:,None,None]*dirs[None,:,:]-c; log=-.5*np.einsum('aki,ij,akj->ak',d,np.linalg.inv(cov),d)
 lp=logsumexp(log,axis=1); density=np.exp(lp-lp.max()); cdf=cumulative_trapezoid(density,radius,initial=0)
 return float(np.interp(.95*cdf[-1],cdf,radius))

def fit(data,limits=True,explore=False):
 jd,psi,y,pf,sig=data.T; t=(jd-jd.mean())/365.25; s=np.sin(psi); c=np.cos(psi)
 x=np.column_stack([s,c,pf,t*s,t*c]); xw=x/sig[:,None]; yw=y/sig
 beta=np.linalg.lstsq(xw,yw,rcond=None)[0]; q=np.linalg.qr(xw)[0]; r=yw-xw@beta; chi=float(r@r); s2=max(1,chi/(len(y)-5)); cov0=np.linalg.inv(xw.T@xw)*s2
 rows=[]
 for p,nom in zip(PERIODS,NOMINAL):
  stats=[]
  for opp in (False,True):
   z=np.column_stack([s,c])*template(p,jd,opp)[:,None]/sig[:,None]; zp=z-q@(q.T@z); cv=np.linalg.inv(zp.T@zp); coeff=cv@zp.T@r; delta=float(coeff@np.linalg.solve(cv,coeff)); stats.append((delta/s2,coeff,cv*s2,zp))
  S,coeff,cv,zp=stats[0]; a=float(np.linalg.norm(coeff)); ratio=a/nom
  rows.append(dict(S=S,S_opp=stats[1][0],a=a,ratio=ratio,node=float(np.degrees(np.arctan2(coeff[0],coeff[1]))%360),M95=mass(upper(coeff,cv),p,beta[2]) if limits else None,BIC_change=2*np.log(len(y))-S,coeff=coeff,orbit_residual=zp@coeff*sig))
 out=dict(n=len(y),sigma=float(np.median(sig)),chi2=chi,s2=s2,parallax=float(beta[2]),parallax_error=float(np.sqrt(cov0[2,2])),rows=rows,residual=r*sig,data=data)
 if explore:
  def improvement(z):
   zw=z/sig[:,None]; zp=zw-q@(q.T@zw); b=np.linalg.lstsq(zp,r,rcond=None)[0]; return float(np.sum((zp@b)**2)/s2)
  out['M1_BIC_change']=2*np.log(len(y))-improvement(np.column_stack([t*t*s/2,t*t*c/2]))
  best=(float('inf'),0)
  for p in 1/np.linspace(1/2500,1/500,129):
   ph=2*np.pi*(jd-jd.mean())/p; z=np.column_stack([s*np.cos(ph),s*np.sin(ph),c*np.cos(ph),c*np.sin(ph)]); bic=5*np.log(len(y))-improvement(z)
   best=min(best,(bic,float(p)))
  out['Mfree']=best
 return out

def detect(r,T): return bool(r['S']>=max(T,19.1) and .4<=r['ratio']<=2 and r['S_opp']<=r['S']+4)
def outcome(rows,good=True):
 if not good: return 'INCONCLUSIVE (data quality)'
 if any(r['DETECT'] for r in rows): return 'CLAIM SUPPORTED'
 if all(r['EXCLUDE'] for r in rows): return 'CLAIM EXCLUDED (strong)'
 if all(r['WEAK_EXCLUDE'] for r in rows): return 'face value (9.4 MJup) excluded'
 return 'INCONCLUSIVE'

def pipeline(directory,min_controls=200,explore=True,injections=True,real=False):
 directory=Path(directory); source=pd.read_csv(directory/'tabby_dr4_source.csv',dtype={'source_id':str}); matches=source[source.source_id==TABBY]
 if len(matches)!=1: raise DataValidationError('ABORT: expected exactly one Tabby source_id matching the DR3 id; check all_source_match')
 source=matches.iloc[0]; tabby_id=source.source_id
 pd.read_csv(directory/'tabby_dr4_flags.csv',dtype={'source_id':str}) # Require the fetched flags product.
 catalog=pd.read_csv(directory/'controls_dr4.csv',dtype={'source_id':str}); catalog=catalog.drop_duplicates('source_id'); catalog=catalog[catalog.source_id.map(lambda sid:(directory/f'ctl_{sid}.csv').is_file())]; catalog=catalog.iloc[np.argsort([int(s) for s in catalog.source_id])][:500]
 epoch_path=directory/f'tabby_epoch_astrometry_{tabby_id}.csv'
 if not epoch_path.is_file(): epoch_path=directory/'tabby_epoch_astrometry.csv'
 tab=fit(collapse(epoch_path,real=real,source_id=tabby_id),explore=explore)
 controls=[]; excluded=[]; retained=[]
 for sid in catalog.source_id:
  try:
   ct=fit(collapse(directory/f'ctl_{sid}.csv',real=real,source_id=sid),limits=False)
   if not all(np.isfinite(r['S']) and np.isfinite(r['S_opp']) for r in ct['rows']): raise ValueError('nonfinite scores')
  except (ValueError,np.linalg.LinAlgError) as exc:
   excluded.append(f'{sid} ({exc})'); continue
  controls.append(ct); retained.append(sid)
 catalog=catalog[catalog.source_id.isin(retained)]
 scores=np.array([[r['S'] for r in ct['rows']] for ct in controls]).reshape(-1,5); fpr=np.mean(scores>13.82,axis=0) if len(controls) else np.full(5,np.nan); thresholds=np.maximum(13.82,np.quantile(scores,.99,axis=0)) if len(controls) else np.full(5,13.82)
 thresholds=np.where(fpr>.015,thresholds,13.82)
 g2=abs(tab['parallax']-source.parallax)<=3*np.hypot(tab['parallax_error'],source.parallax_error)
 gates={'G1':bool(tab['n']>=60 and source.ruwe<1.4),'G2':bool(g2),'G3':bool(len(controls)>=min_controls and np.all(fpr<=.015)),'G4':bool(len(controls) and all(source[k]<=catalog[k].quantile(.95) for k in ['ipd_frac_multi_peak','ipd_gof_harmonic_amplitude']))}
 diagnostic=''
 if not g2:
  alt=fit(collapse(epoch_path,-BARY_SIGN,real=real,source_id=tabby_id),limits=False)
  diagnostic=f"G2 diagnostic: +sign parallax {tab['parallax']:.7f}; -sign {alt['parallax']:.7f}; catalogue {source.parallax:.7f} mas. Weak diagnostic: a sign flip need not fail G2."
  print(diagnostic)
 for r,T in zip(tab['rows'],thresholds):
  r.update(T=float(T),DETECT=detect(r,T)); r['EXCLUDE']=not r['DETECT'] and r['M95']<5.3; r['WEAK_EXCLUDE']=r['M95']<9.4
 recovery=[]; rng=np.random.default_rng(20261202)
 if injections:
  for j,p in enumerate(PERIODS):
   detected=[]; ratios=[]
   for ct in controls:
    d=ct['data'].copy(); node=rng.uniform(0,2*np.pi); a=amplitude(9.4,p,ct['parallax']); psi=d[:,1]; d[:,2]+=a*template(p,d[:,0],True)*(np.sin(node)*np.sin(psi)+np.cos(node)*np.cos(psi))
    try: r=fit(d,limits=False)['rows'][j]
    except (ValueError,np.linalg.LinAlgError): continue
    if not (a>0 and np.isfinite(a) and np.isfinite(r['S']) and np.isfinite(r['S_opp'])): continue
    r['ratio']=r['a']/a
    detected.append(detect(r,thresholds[j])); ratios.append(r['ratio'])
   recovery.append((float(np.mean(detected)) if detected else np.nan,float(np.median(ratios)) if ratios else np.nan))
 nss=[]
 for path in sorted(directory.glob('tabby_dr4_nss_*.csv')):
  df=pd.read_csv(path,dtype={'source_id':str})
  if 'source_id' in df: df=df[df.source_id==tabby_id]
  for _,row in df.iterrows():
   p=number(row.get('period',np.nan)); matches=[f'W{i+1}' for i,w in enumerate(PERIODS) if abs(p/w-1)<=.05]
   nss.append(f"{path.stem}: solution_type={row.get('solution_type','unspecified')}, period={p}, matches={','.join(matches) or 'none'}")
 return dict(excluded=excluded,alternative=outcome(tab['rows'],all(v for k,v in gates.items() if k!='G3')),tab=tab,source=source,controls=controls,scores=scores,fpr=fpr,thresholds=thresholds,gates=gates,outcome=outcome(tab['rows'],all(gates.values())),recovery=recovery,nss=nss,diagnostic=diagnostic)

def report(result,out,plot=None,drop=None,citation=None,dry=False,hash_ok=True):
 r=result; t=r['tab']; rows=t['rows']; src=r['source']; lines=[f"# {r['outcome']}",f"Synthetic dry run: {dry}. Prereg hash match: {hash_ok}.",f"prereg.md SHA256: {sha(ROOT/'prereg.md')}",f"02_fit.py SHA256: {sha(__file__)}",f"N FoV={t['n']}; median sigma={t['sigma']:.6f} mas; RUWE={src.ruwe}; chi2/dof={t['chi2']:.3f}/{t['n']-5}; s2={t['s2']:.4f}.",f"Parallax fit={t['parallax']:.7f} +/- {t['parallax_error']:.7f}; catalogue={src.parallax:.7f} +/- {src.parallax_error:.7f} mas.",'','|Window|S|S opposite|T|a mas|a/alpha_nom|node PA deg|M95 MJup|DETECT|EXCLUDE|WEAK-EXCLUDE|delta BIC|','|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|']
 for i,x in enumerate(rows): lines.append(f"|W{i+1}|{x['S']:.4f}|{x['S_opp']:.4f}|{x['T']:.4f}|{x['a']:.6f}|{x['ratio']:.4f}|{x['node']:.2f}|{x['M95']:.4f}|{x['DETECT']}|{x['EXCLUDE']}|{x['WEAK_EXCLUDE']}|{x['BIC_change']:.4f}|")
 lines += ['',f"alternative reading (thresholds raised, G3 not vetoing): {r['alternative']}",'Excluded control IDs: '+('; '.join(r['excluded']) or 'none')+'.', 'Barycentric +/- sign diagnostic is weak: a sign flip need not fail G2.']
 lines+=['','|Gate|Pass|Cause / criterion|','|---|---|---|']
 descriptions={'G1':'>=60 FoV and RUWE<1.4','G2':'parallax agreement within 3 combined sigma','G3':'>=200 controls and every FPR<=1.5%','G4':'both IPD values <= control 95th percentiles'}
 lines += [f'|{k}|{v}|{descriptions[k]}|' for k,v in r['gates'].items()]
 lines+=['',f"Controls n={len(r['controls'])}; median s2={np.median([c['s2'] for c in r['controls']]) if r['controls'] else np.nan:.4f}.",'|Window|FPR S>13.82|9.4 MJup recovery|median a_fit/alpha|','|---|---:|---:|---:|']
 for j in range(5):
  rec=r['recovery'][j] if r['recovery'] else (np.nan,np.nan); lines.append(f'|W{j+1}|{r["fpr"][j]:.4f}|{rec[0]:.4f}|{rec[1]:.4f}|')
 lines += ['',f"Exploratory M1 delta BIC={t.get('M1_BIC_change')}; Mfree (delta BIC, period d)={t.get('Mfree')}. BIC uses rescaled chi2; Mfree counts the searched period as a fifth added parameter.",*(r['nss'] or ['NSS: no Tabby rows available.']),r['diagnostic']]
 if drop:
  retained=[x for j,x in enumerate(rows) if f'W{j+1}' not in drop]; lines += ['',f"Restricted outcome: {outcome(retained,all(r['gates'].values()))}. Dropped: {','.join(drop)}. Citation supplied by operator: {citation}. Publication date and transit non-detection eligibility must be verified before interpreting this restricted result."]
 lines+=['','Caveats: draft DR4 serialization/units and barycentric sign require Dec 2 verification; BJD TCB is used directly with the supplied transit anchor. Controls injection uses opposite-convention truth and is reported only. Bayesian limits use flat amplitude and uniform orientation; other parameters profiled. Null simulations are Gaussian, not a validation of real DR4 systematics. No post hoc changes to prereg.md.']
 out=Path(out); out.parent.mkdir(parents=True,exist_ok=True); out.write_text('\n'.join(lines)+'\n')
 if plot:
  import matplotlib; matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  fig,axes=plt.subplots(1,3,figsize=(15,4)); ax=axes[0]
  for j in range(5):
   ax.hist(r['scores'][:,j],bins=20,histtype='step',label=f'W{j+1}'); ax.axvline(rows[j]['S'],color=f'C{j}',ls='--')
  ax.set(xlabel='S (histograms controls; dashed Tabby)',ylabel='Count'); ax.legend(fontsize=7)
  best=int(np.argmax([x['S'] for x in rows])); d=t['data']; phase=((d[:,0]-2458730.164)/PERIODS[best])%1; ix=np.argsort(phase)
  axes[1].errorbar(phase,t['residual'],yerr=d[:,4],fmt='.',alpha=.5); axes[1].plot(phase[ix],rows[best]['orbit_residual'][ix],'.',label='Mc projected AL fit'); axes[1].set(xlabel=f'W{best+1} template phase',ylabel='AL residual after M0 (mas)'); axes[1].legend(fontsize=7)
  axes[2].plot(range(1,6),[x['M95'] for x in rows],'o-'); axes[2].axhline(5.3,ls='--'); axes[2].axhline(9.4,ls=':'); axes[2].set(xlabel='Window',ylabel='95% mass limit (MJup)',xticks=range(1,6)); fig.tight_layout(); Path(plot).parent.mkdir(parents=True,exist_ok=True); fig.savefig(plot); plt.close(fig)

def main():
 p=argparse.ArgumentParser(); p.add_argument('--data-dir',default='data_dr4'); p.add_argument('--out',default='memo/phaseB_result.md'); p.add_argument('--plot',default='results/phaseB_fit.png'); p.add_argument('--dry-run',action='store_true'); p.add_argument('--case',choices=['injected','null'],default='injected'); p.add_argument('--seed',type=int,default=1); p.add_argument('--drop-windows'); p.add_argument('--drop-citation'); p.add_argument('--hashes-file',default=str(ROOT/'results/hashes.txt')); p.add_argument('--fmt',choices=['gaia_array','standard'],default='gaia_array'); p.add_argument('--skip-plots',action='store_true')
 a=p.parse_args(); drop=a.drop_windows.split(',') if a.drop_windows else None
 if drop and (not a.drop_citation or any(x not in [f'W{i}' for i in range(1,6)] for x in drop) or len(set(drop))==5): p.error('restricted outcome requires a citation, valid windows, and at least one retained window')
 ok=hash_status(a.hashes_file); print(f'Prereg hash match: {ok}')
 fit_ok=fitter_hash_status(a.hashes_file); print(f'02_fit.py hash match: {fit_ok}' if fit_ok is not None else 'WARNING: phaseB/02_fit.py entry missing from hashes file')
 if not a.dry_run:
  if not ok: p.exit(2,'ABORT: prereg hash mismatch\n')
  if (Path(a.data_dir)/'SYNTHETIC_MARKER').exists(): p.exit(2,'ABORT: SYNTHETIC_MARKER in real-mode data\n')
  if fit_ok is False: p.exit(2,'ABORT: phaseB/02_fit.py hash mismatch\n')
  try: r=pipeline(a.data_dir,real=True)
  except DataValidationError as exc: p.exit(2,str(exc)+'\n')
  report(r,a.out,None if a.skip_plots else a.plot,drop,a.drop_citation,hash_ok=ok)
 else:
  with tempfile.TemporaryDirectory(prefix='phaseb-') as td:
   make_dataset(td,a.case,a.seed,fmt=a.fmt); r=pipeline(td); report(r,ROOT/f'memo/phaseB_dryrun_{a.case}.md',None if a.skip_plots else ROOT/f'results/phaseB_dryrun_{a.case}.png',drop,a.drop_citation,True,ok)
 print(r['outcome']); print('S:',', '.join(f"W{i+1}={x['S']:.3f}" for i,x in enumerate(r['tab']['rows']))); print('Gates:',r['gates'])
if __name__=='__main__': main()
