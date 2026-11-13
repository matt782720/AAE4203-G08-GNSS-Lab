"""One-factor weighting experiment; no hidden reference trajectory is used."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'course'))
import WLSE_spp as course
from process_g08 import rot
def solve(sp,corrected,seed,weighted):
    x=np.r_[seed,0.]
    for _ in range(30):
        d=sp-x[:3];r=np.linalg.norm(d,axis=1);H=np.column_stack((-d/r[:,None],np.ones(len(r))));sag=7.2921151467e-5/299792458*(sp[:,0]*x[1]-sp[:,1]*x[0]);w=np.ones(len(r))
        if weighted and np.linalg.norm(x[:3])>6e6:w=np.maximum(((d@rot(x[:3]).T)[:,2]/r),.05)**2
        scale=np.sqrt(w);dx,_,rank,_=np.linalg.lstsq(H*scale[:,None],(corrected-r-sag-x[3])*scale,rcond=None)
        if rank<4:raise ValueError('rank')
        x+=dx
        if np.linalg.norm(dx)<1e-4:
            if 6.3e6<np.linalg.norm(x[:3])<6.42e6:return x[:3]
            break
    raise ValueError('convergence')
def run(folder,weighted):
    a,epochs=course.read_inputs(folder);seed=np.zeros(3);rows=[]
    for i,(week,tow) in enumerate(epochs):
        sp=a['satellite_positions'][:,i*3:i*3+3];v=[a[n][:,i] for n in ['pseudoranges_meas','satellite_clock_bias','ionospheric_delay','tropospheric_delay','group_delay']];mask=np.isfinite(sp).all(1)
        for z in v:mask&=np.isfinite(z)
        if mask.sum()<4:continue
        pr,cl,I,T,g=[z[mask] for z in v]
        try:p=solve(sp[mask],pr+cl-I-T-g,seed,weighted)
        except ValueError:continue
        seed=p;rows.append([week,tow,*p])
    return np.array(rows)
def main():
    results=[];ref=pd.read_csv(ROOT/'data/NAV-PVT.csv');lat=np.radians(ref.lat.to_numpy());lon=np.radians(ref.lon.to_numpy());h=ref.height.to_numpy()*.001;n=6378137/np.sqrt(1-.00669437999014*np.sin(lat)**2);xyz=np.array([(n+h)*np.cos(lat)*np.cos(lon),(n+h)*np.cos(lat)*np.sin(lon),(n*(1-.00669437999014)+h)*np.sin(lat)]).T
    # Raw pyubx2 iTOW is numeric milliseconds in the exported dictionary.
    refmap={round(float(t)*.001,1):x for t,x in zip(ref.iTOW,xyz)};R=rot(xyz.mean(0))
    for dataset,folder in [('own',ROOT/'data'),('benchmark',ROOT/'course/benchmark/gps_csv')]:
        variants=[]
        for weighted in [False,True]:
            a=run(folder,weighted);variant='elevation_WLS' if weighted else 'uniform_LS';pd.DataFrame(a,columns=course.HEADER).to_csv(ROOT/f'data/{dataset}_{variant}.csv',index=False);record=dict(dataset=dataset,variant=variant,epochs=len(a))
            if dataset=='own':
                d=np.array([(p[2:]-refmap[round(p[1],1)])@R.T for p in a if round(p[1],1) in refmap]);record.update(receiver_matched=len(d),receiver_difference_rmse2d_m=float(np.sqrt(np.mean(np.sum(d[:,:2]**2,1)))),receiver_difference_rmse3d_m=float(np.sqrt(np.mean(np.sum(d**2,1)))))
            else:record.update(accuracy='Not available: instructor reference withheld')
            results.append(record);variants.append(a)
        assert np.array_equal(variants[0][:,:2],variants[1][:,:2]);delta=variants[1][:,2:]-variants[0][:,2:];results[-1]['change_from_uniform_3d_rms_m']=float(np.sqrt(np.mean(np.sum(delta**2,1))))
    pd.DataFrame(results).to_csv(ROOT/'data/ablation.csv',index=False)
    own=pd.read_csv(ROOT/'data/own_solution_full.csv').to_numpy();matrix=pd.read_csv(ROOT/'matrix_check/matrix_check.csv').to_numpy();assert np.allclose(own[:,:2],matrix[:,:2],atol=.01,rtol=0)
    difference=float(np.max(np.abs(own[:,2:]-matrix[:,2:])));assert difference<.01,difference
    (ROOT/'data/verification.json').write_text(json.dumps(dict(matrix_solver_max_coordinate_difference_m=difference,ablation=results),indent=2));print(json.dumps(results,indent=2))
if __name__=='__main__':main()
