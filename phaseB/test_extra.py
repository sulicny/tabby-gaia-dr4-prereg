import sys,tempfile,json
sys.path.insert(0,'phaseB')
from test_fit import f
import pandas as pd
import numpy as np
with tempfile.TemporaryDirectory(prefix='phaseb-extra-') as td:
 d=f.Path(td); f.make_dataset(d,'null',7,n_controls=200)
 for sid in range(100000,100010):
  path=d/f'ctl_{sid}.csv'; df=pd.read_csv(path)
  for i,row in df.iterrows():
   times=np.array(json.loads(row.obs_time_tcb)); jd=2455197.5+(times+row.obs_time_bary_corr)/86400e9; psi=np.deg2rad(json.loads(row.scan_pos_angle)); values=np.array(json.loads(row.centroid_pos_al)); values+=.15*f.template(f.PERIODS[3],jd)*np.sin(psi); df.at[i,'centroid_pos_al']=json.dumps(values.tolist())
  df.to_csv(path,index=False)
 r=f.pipeline(d,explore=False,injections=False)
 assert r['gates']['G3'] and all(r['gates'].values()) and not r['outcome'].startswith('INCONCLUSIVE')
 assert np.any(r['fpr']>.015) and np.all(r['thresholds']>=13.82)
 assert np.all(r['thresholds'][r['fpr']>.015]>13.82)
 print('Contaminated controls: FPR',r['fpr'],'thresholds',r['thresholds'],'G3 passes without quality veto PASS')
 # Serialization variants, padding and boolean masks exercised through actual collapse.
 path=d/'tabby_epoch_astrometry.csv'; df=pd.read_csv(path); baseline=f.collapse(path)
 for key in ['obs_time_tcb','scan_pos_angle','centroid_pos_al','centroid_pos_error_al','used_by_agis_al']:
  def convert(text):
   values=json.loads(text)
   if key=='used_by_agis_al': values=[('t' if x else 'f') for x in values]+['null']
   else: values=[str(x) for x in values]+['NaN']
   return '{'+','.join(values)+'}'
  df[key]=df[key].map(convert)
 df.to_csv(path,index=False); assert np.allclose(f.collapse(path),baseline)
 print('Brace arrays, t/f booleans, NaN/null padding through collapse PASS')
 # Restricted result preserves unrestricted headline.
 out=d/'restricted.md'; f.report(r,out,drop=['W1'],citation='test-only citation')
 assert out.read_text().startswith('# '+r['outcome']) and 'Restricted outcome:' in out.read_text()
 print('Restricted-result headline preservation PASS')
 # Independently solve the full seven-parameter model for both phases.
 data=r['tab']['data']; jd,psi,y,pf,sig=data.T; t=(jd-jd.mean())/365.25
 x=np.column_stack([np.sin(psi),np.cos(psi),pf,t*np.sin(psi),t*np.cos(psi)])
 for j,p in enumerate(f.PERIODS):
  row=r['tab']['rows'][j]
  for opposite in [False,True]:
   z=np.column_stack([np.sin(psi),np.cos(psi)])*f.template(p,jd,opposite)[:,None]
   design=np.column_stack([x,z])/sig[:,None]; coeff=np.linalg.lstsq(design,y/sig,rcond=None)[0]
   chi=float(np.sum((y/sig-design@coeff)**2)); scale=max(1,chi/(len(y)-7))
   assert np.isclose(row['chi2_Mc_opp' if opposite else 'chi2_Mc'],chi)
   assert np.isclose(row['s2_opp' if opposite else 's2_W'],scale)
   assert np.isclose(row['S_opp' if opposite else 'S'],(r['tab']['chi2']-chi)/scale)
   if not opposite:
    cov=np.linalg.inv(design.T@design)[-2:,-2:]*scale
    expected=f.mass(f.upper(coeff[-2:],cov),p,r['tab']['parallax'])
    assert np.isclose(row['M95'],expected)
 print('Full-model residual, phase-specific rescale, statistic and covariance UL PASS')
