#!/usr/bin/env python
"""Reproducible offline T1--T4; full run writes results/phaseB_tests_v2.md."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import importlib.util,sys,tempfile,subprocess,json,time
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
spec=importlib.util.spec_from_file_location('fit',Path(__file__).with_name('02_fit.py')); f=importlib.util.module_from_spec(spec); spec.loader.exec_module(f)
PYTHON=sys.executable
def run():
 lines=['# Phase B offline tests','Run: '+pd.Timestamp.now(tz='America/Detroit').strftime('%Y-%m-%d %H:%M ET'), 'Commands: `python phaseB/test_fit.py`; `python phaseB/02_fit.py --dry-run --case injected|null`. No network; prereg.md unchanged.']
 start=time.monotonic(); evidence={}
 with tempfile.TemporaryDirectory(prefix='phaseb-tests-') as tmp:
  root=Path(tmp)
  lines+=['## T1 — end-to-end, both serializations']
  for case in ['injected','null']:
   outcomes=[]
   for fmt in ['gaia_array','standard']:
    d=root/f'{case}-{fmt}'; f.make_dataset(d,case,1,fmt=fmt); r=f.pipeline(d); outcomes.append(r)
    f.report(r,f.ROOT/f'memo/phaseB_dryrun_{case}.md',None,dry=True)
    lines.append(f"- {case}, {fmt}: **{r['outcome']}**; S="+str([round(x['S'],4) for x in r['tab']['rows']])+f"; gates={r['gates']}.")
    assert all(r['gates'].values()),r['gates']
   assert np.allclose([x['S'] for x in outcomes[0]['tab']['rows']],[x['S'] for x in outcomes[1]['tab']['rows']],rtol=1e-8)
   evidence[case]=d
  lines+=['## T2 — 40 independent seeds per truth']
  # Reduced simulations relax only the valid-control count to 100 for both truths.
  records={}
  for case in ['null','injected']:
   det=[]; w4=[]; limits=[]; gatepass=[]; scales=[]
   for seed in range(40):
    d=root/'repeat'; f.make_dataset(d,case,1000+seed,n_controls=100)
    # Remove stale higher-ID files/catalogue is authoritative for selection.
    r=f.pipeline(d,min_controls=100,explore=False,injections=False)
    rows=r['tab']['rows']; det.append(any(x['DETECT'] for x in rows)); w4.append(rows[3]['DETECT']); limits.append(rows[3]['M95']); gatepass.append(all(r['gates'].values())); scales.append(rows[3]['s2_W'])
    print(case,seed+1,'/40',flush=True)
   records[case]=dict(detections=sum(det),fraction=float(np.mean(det)),W4_fraction=float(np.mean(w4)),median_W4_M95=float(np.median(limits)),all_gates_pass=sum(gatepass),all_gates_fraction=float(np.mean(gatepass)),median_s2_W4=float(np.median(scales)))
   lines.append(f'- {case}: {records[case]}. DETECT is the per-window rule, before the overall quality veto; all-gates counts are also shown.')
  lines+=['## Claim-true posterior — 60 seeds', 'Seeds 2000–2059; independent RNG 20261112 selects uniform W1–W5 and K=162+z*(69 if z<0 else 76), rejecting K<0. Mass maps K at the window midpoint, e=0.09, sin i=1; generator draws the remaining truth parameters. Central noise=0.75; 100 controls and relaxed G3.']
  rng=np.random.default_rng(20261112); posterior=[]
  from gen_sim_dr4 import alpha,mass_from_a
  for seed in range(2000,2060):
   window=int(rng.integers(1,6)); K=-1
   while K<0:
    z=rng.normal(); K=162+z*(69 if z<0 else 76)
   m=float(mass_from_a(alpha(K,f.PERIODS[window-1],.09,1),f.PERIODS[window-1]))
   d=root/'posterior'; f.make_dataset(d,'injected',seed,n_controls=100,window=window,mass_mjup=m)
   r=f.pipeline(d,min_controls=100,explore=False,injections=False); rows=r['tab']['rows']
   posterior.append(dict(window=window,mass=m,any=any(x['DETECT'] for x in rows),own=rows[window-1]['DETECT'],raw=rows[window-1]['S']*rows[window-1]['s2_W']>=13.82,gates=all(r['gates'].values())))
   print('posterior',seed-1999,'/60',flush=True)
  records['posterior']=posterior
  lines+=['|Truth|n|Any DETECT|Own DETECT|Own raw >=13.82|Central raw forecast|','|---|---:|---:|---:|---:|---:|']
  for w,forecast in [(0,'—'),(1,'~0.70'),(2,'~0.70'),(3,'0.80'),(4,'0.81'),(5,'0.84')]:
   subset=[x for x in posterior if w==0 or x['window']==w]
   lines.append(f"|{'Overall' if w==0 else 'W'+str(w)}|{len(subset)}|{np.mean([x['any'] for x in subset]):.3f}|{np.mean([x['own'] for x in subset]):.3f}|{np.mean([x['raw'] for x in subset]):.3f}|{forecast}|")
  lines.append('Forecast uses a raw statistic and nominal central assumptions; full DETECT adds amplitude, opposite-phase and 19.1 criteria. Small per-window samples limit comparison precision.')
  lines+=['|Fixed truth|Any DETECT|W4 DETECT|Median W4 UL|Median s2 W4|All gates fraction|','|---|---:|---:|---:|---:|---:|']
  for case in ['null','injected']:
   v=records[case]; lines.append(f"|{case}|{v['fraction']:.3f}|{v['W4_fraction']:.3f}|{v['median_W4_M95']:.3f}|{v['median_s2_W4']:.3f}|{v['all_gates_fraction']:.3f}|")
  lines+=['## T3 — empirical DR3 scales and gate failures']
  d=evidence['null']; r=f.pipeline(d,injections=False,explore=False)
  assert r['gates']['G3']; scales=np.array([c['s2'] for c in r['controls']]); lines.append(f"- DR3-derived per-star scales: FPR={r['fpr'].tolist()}, G3 PASS; control s2 min/median/max={np.min(scales):.4f}/{np.median(scales):.4f}/{np.max(scales):.4f}. Errors track the per-star scale, so s2 is near 1, clipped below at 1.")
  srcpath=d/'tabby_dr4_source.csv'; src=pd.read_csv(srcpath); original=src.copy(); src['ruwe']=1.6; src.to_csv(srcpath,index=False)
  bad=f.pipeline(d,injections=False,explore=False); assert not bad['gates']['G1'] and bad['outcome']=='INCONCLUSIVE (data quality)'; lines.append('- RUWE 1.6: G1 failed and headline INCONCLUSIVE (data quality), PASS.')
  original.to_csv(srcpath,index=False); ep=d/'tabby_epoch_astrometry.csv'; raw=pd.read_csv(ep); saved=raw.copy(); raw=raw[raw.transit_id.isin(raw.transit_id.unique()[:40])]; raw.to_csv(ep,index=False)
  bad=f.pipeline(d,injections=False,explore=False); assert bad['tab']['n']==40 and not bad['gates']['G1'] and bad['outcome']=='INCONCLUSIVE (data quality)'; lines.append('- 40 FoV: G1 failed and quality veto applied, PASS.'); saved.to_csv(ep,index=False)
  src=original.copy(); src['parallax']+=1; src.to_csv(srcpath,index=False)
  bad=f.pipeline(d,injections=False,explore=False); assert not bad['gates']['G2'] and bad['outcome']=='INCONCLUSIVE (data quality)'; lines.append('- Catalogue parallax offset +1 mas: G2 failed, both barycentric signs diagnosed, quality veto applied, PASS.'); original.to_csv(srcpath,index=False)
  src=original.copy(); src['ipd_gof_harmonic_amplitude']=10; src.to_csv(srcpath,index=False)
  bad=f.pipeline(d,injections=False,explore=False); assert not bad['gates']['G4'] and bad['outcome']=='INCONCLUSIVE (data quality)'; lines.append('- IPD harmonic amplitude 10: G4 failed and quality veto applied, PASS.'); original.to_csv(srcpath,index=False)
  sparse=d/'controls_dr4.csv'; ctl=pd.read_csv(sparse); ctl.iloc[:199].to_csv(sparse,index=False); bad=f.pipeline(d,injections=False,explore=False); assert not bad['gates']['G3']; ctl.to_csv(sparse,index=False); lines.append('- Only 199 controls: G3 failed, PASS.')
  # Deliberately underreported control errors: rescaling must contain the false positives.
  c=r['controls'][0]['data'].copy(); baseline=f.fit(c,limits=False); c[:,4]/=4; scaled=f.fit(c,limits=False)
  assert np.isclose(scaled['s2'],16*baseline['chi2']/(len(c)-5))
  assert np.allclose([x['S'] for x in scaled['rows']],[x['S'] for x in baseline['rows']],rtol=.2)
  lines.append(f"- Errors divided by 4: s2={scaled['s2']:.4f}; S uses each Mc residual rescale (M0 diagnostic shown), PASS.")
  lines+=['## T4 — real-mode safety']
  # Isolate rule files in a temporary ROOT; never create the operator's prereg_v2.md.
  fixture=root/'hash-fixture'; fixture.mkdir(); (fixture/'results').mkdir(); (fixture/'phaseB').mkdir()
  import shutil
  for name in ['02_fit.py','gen_sim_dr4.py']: shutil.copy(f.ROOT/'phaseB'/name,fixture/'phaseB'/name)
  shutil.copy(f.ROOT/'prereg.md',fixture/'prereg.md'); (fixture/'prereg_v2.md').write_text('test-only v2 fixture\n')
  entries=['prereg.md','prereg_v2.md','phaseB/02_fit.py']
  for bad,label in [('prereg.md','prereg hash mismatch'),('prereg_v2.md','prereg_v2.md hash mismatch'),('phaseB/02_fit.py','phaseB/02_fit.py hash mismatch'),(None,'SYNTHETIC_MARKER'),('missing-v2','prereg_v2.md hash mismatch'),('missing-fit','phaseB/02_fit.py hash mismatch')]:
   hashes=fixture/'results/hashes.txt'
   text='0'*64+'  phaseB/02_fit.py\n' # stale first entry must be ignored
   for entry in entries:
    if (bad=='missing-v2' and entry=='prereg_v2.md') or (bad=='missing-fit' and entry=='phaseB/02_fit.py'): continue
    text+=('0'*64 if entry==bad else f.sha(fixture/entry))+'  '+entry+'\n'
   if bad=='missing-fit': text=''.join(line+'\n' for line in text.splitlines() if 'phaseB/02_fit.py' not in line)
   hashes.write_text(text); output=root/'must-not-exist.md'
   proc=subprocess.run([PYTHON,str(fixture/'phaseB/02_fit.py'),'--data-dir',str(d),'--out',str(output)],capture_output=True,text=True)
   assert proc.returncode==2 and label in proc.stderr and not output.exists(),proc
   lines.append(f'- {bad or label}: exit 2; no result written, PASS.')
  # Mathematical, serialization and independence checks.
  from scipy.special import erfinv
  assert abs(f.upper(np.zeros(2),np.eye(2))/(np.sqrt(2)*erfinv(.95))-1)<.002
  for p in f.PERIODS:
   for opp in (False,True): assert abs(f.template(p,np.array([2458730.164]),opp)[0])<1e-12
   assert abs(f.mass(f.amplitude(9.4,p,2.254545),p,2.254545)-9.4)<1e-8
  for text in ['{1,null,NaN}','[1, null, NaN]','(1,null,NaN)']: assert len(f.cell(text))==3 and np.isnan(f.number(f.cell(text)[1]))
  import ast
  ref=ast.parse((f.ROOT/'scripts/sim.py').read_text()); gen=ast.parse((f.ROOT/'phaseB/gen_sim_dr4.py').read_text())
  for name in ['alpha','kval','mass_from_a','orbit']:
   find=lambda tree:next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
   assert ast.dump(find(ref))==ast.dump(find(gen))
  lines+=['## Additional checks','- Analytic flat-amplitude Gaussian 95% limit, exact mass inversion, zero displacement at transit for both conventions, null/NaN parsing, and AST equality of copied generator orbital functions: PASS.']
  (f.ROOT/'results/phaseB_test_metrics_v2.json').write_text(json.dumps(records,indent=2))
  lines += ['## Fixed-mass comparison',f"- Null W4 median UL {records['null']['median_W4_M95']:.3f} MJup versus forecast ~4.1. Injected DETECT fraction {records['injected']['fraction']:.3f} versus ~0.98 Phase A raw-statistic power. The posterior comparison above uses a different estimand.",f'Elapsed {time.monotonic()-start:.1f} seconds.']
 lines += ['## v1 comparison', '|Rule|Injected any|Injected W4|Null any FPR|Null median W4 UL|', '|---|---:|---:|---:|---:|', '|v1|0.675|0.600|0/40|3.89|', f"|v2|{records['injected']['fraction']:.3f}|{records['injected']['W4_fraction']:.3f}|{records['null']['detections']}/40|{records['null']['median_W4_M95']:.3f}|", 'Reduced runs explicitly relax G3 to 100 controls for both truths; production requires 200. DETECT fractions precede the quality veto. No network. Real DR4 performance remains unverified.']
 extra=subprocess.run([PYTHON,str(f.ROOT/'phaseB/test_extra.py')],check=True,capture_output=True,text=True)
 lines+=['## G3 raised-threshold and independent fit checks','```text',extra.stdout.strip(),'```']
 (f.ROOT/'results/phaseB_tests_v2.md').write_text('\n'.join('\n'+x if x.startswith('## ') else x for x in lines)+'\n')
 print('\n'.join(lines),flush=True)
if __name__=='__main__': run()
