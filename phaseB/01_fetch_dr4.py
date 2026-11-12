#!/usr/bin/env python
"""Phase B step 1: fetch Gaia DR4 catalogue rows + epoch astrometry for Tabby's star and ~controls.
Anonymous access only (no login). Run on/after 2026-12-02. Usage: python 01_fetch_dr4.py [--dry]
TAP schema/table names follow the DRAFT DR4 data model (2026-06-26): gaiadr4.gaia_source,
all_source_astrometry, all_source_flags, nss_two_body_orbit, nss_acceleration_astro, nss_masses.
VERIFY names via TAP_SCHEMA once live (step 0 below prints them).
Lessons from Phase A (DR3 archive, anonymous): async jobs are killed at 60 s; cone searches on
ra/dec time out on gaia_source; use source_id ranges (HEALPix level-12 id >> 35, indexed) instead.
"""
import requests, time, io, sys, os, pandas as pd, numpy as np, healpy as hp
TAP="https://gea.esac.esa.int/tap-server/tap"; DL="https://gea.esac.esa.int/data-server/data"
TABBY=2081900940499099136; NEIGH=2081900944807842560     # DR3 ids; DR4 ids normally identical (verify via all_source_match / DR3 neighbours table)
OUT=os.path.join(os.path.dirname(__file__),'..','data_dr4'); os.makedirs(OUT,exist_ok=True)
S=requests.Session()
def adql(q,name,timeout=900):
    r=S.post(TAP+"/async",data=dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",PHASE="RUN",QUERY=q),allow_redirects=False,timeout=60)
    url=r.headers['Location']; t=time.time()
    while True:
        ph=S.get(url+"/phase",timeout=60).text.strip()
        if ph in('COMPLETED','ERROR','ABORTED'): break
        if time.time()-t>timeout: raise TimeoutError(name)
        time.sleep(5)
    if ph!='COMPLETED': raise RuntimeError(S.get(url,timeout=60).text[:3000])
    d=pd.read_csv(io.StringIO(S.get(url+"/results/result",timeout=300).text)); d.to_csv(f'{OUT}/{name}.csv',index=False); return d
def sid_range(ra,dec,radius_deg,level=6):
    pix=sorted(hp.query_disc(2**level,hp.ang2vec(ra,dec,lonlat=True),np.radians(radius_deg),nest=True,inclusive=True))
    return " OR ".join(f"(source_id>={p*2**35*4**(12-level)} AND source_id<{(p+1)*2**35*4**(12-level)})" for p in pix)
def epoch_astrometry(sid,fmt="csv",structure="DATAMODEL_GAIA"):
    r=S.get(DL,params=dict(ID=f"Gaia DR4 {sid}",RETRIEVAL_TYPE="EPOCH_ASTROMETRY",DATA_STRUCTURE=structure,FORMAT=fmt),timeout=300)
    r.raise_for_status(); return r.content
if __name__=="__main__":
    if '--dry' in sys.argv: print(sid_range(301.5643,44.4568,2.5)[:200]); sys.exit()
    # step 0: confirm schema + tables exist
    try: print(adql("SELECT table_name FROM TAP_SCHEMA.tables WHERE schema_name='gaiadr4'","dr4_tables").table_name.tolist())
    except Exception as e: print("schema listing failed",str(e)[:300])
    cols="source_id,ra,dec,phot_g_mean_mag,bp_rp,parallax,parallax_error,pmra,pmdec,ruwe,astrometric_params_solved,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_n_obs_al,astrometric_n_good_obs_al,astrometric_chi2_al,visibility_periods_used,ipd_gof_harmonic_amplitude,ipd_frac_multi_peak,ipd_frac_odd_win,duplicated_source,nu_eff_used_in_astrometry"
    adql(f"SELECT {cols} FROM gaiadr4.gaia_source WHERE source_id IN ({TABBY},{NEIGH})","tabby_dr4_source")
    # all_source_flags for has_epoch_astrometry; NSS tables
    adql(f"SELECT source_id,has_epoch_astrometry FROM gaiadr4.all_source_flags WHERE source_id IN ({TABBY},{NEIGH})","tabby_dr4_flags")
    for t in ("nss_two_body_orbit","nss_acceleration_astro","nss_masses"):
        try: adql(f"SELECT * FROM gaiadr4.{t} WHERE source_id IN ({TABBY},{NEIGH})",f"tabby_dr4_{t}")
        except Exception as e: print(t,"failed",str(e)[:200])
    # epoch astrometry (DataLink, one source per call; anonymous; <=5000 sources per call per ESA docs)
    open(f'{OUT}/tabby_epoch_astrometry.csv','wb').write(epoch_astrometry(TABBY))
    # controls: same recipe as Phase A (healpix-range selection, G 11.2-12.4, |bp_rp-0.776|<0.25, 5p, RUWE<1.4), 
    ctl=adql(f"SELECT TOP 2500 {cols} FROM gaiadr4.gaia_source WHERE ({sid_range(301.5643,44.4568,2.5)}) AND phot_g_mean_mag BETWEEN 11.2 AND 12.4 AND bp_rp BETWEEN 0.5257 AND 1.0257 AND astrometric_params_solved=31 AND ruwe<1.4 AND source_id<>{TABBY}","controls_dr4")
    # pick controls deterministically (prereg: lowest 600 source_id mod hash -> first 500 with epoch astrometry), fetch their epoch astrometry sequentially w/ resume
    ctl=ctl.sort_values('source_id').iloc[:500]
    for sid in ctl.source_id:
        f=f'{OUT}/ctl_{sid}.csv'
        if os.path.exists(f): continue
        for k in range(4):
            try: open(f,'wb').write(epoch_astrometry(int(sid))); break
            except Exception as e: time.sleep(10*(k+1))
