"""GPS L1 SPP from checksum-validated UBX RAWX and broadcast RINEX navigation.
Run from the G08_lab directory: python code/process_g08.py
Dependencies: numpy, pandas, matplotlib, pyubx2. No receiver coordinates enter SPP.
"""
from pathlib import Path
import struct, csv, json, math
from collections import Counter
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyubx2 import UBXReader

ROOT=Path(__file__).resolve().parents[1]
C=299792458.; OM=7.2921151467e-5; MU=3.986005e14
def llh(x):
    lon=np.arctan2(x[1],x[0]); p=np.hypot(*x[:2]); lat=np.arctan2(x[2],p*(1-0.00669437999014))
    for _ in range(10):
        n=6378137/np.sqrt(1-0.00669437999014*np.sin(lat)**2)
        lat=np.arctan2(x[2]+0.00669437999014*n*np.sin(lat),p)
    return lat,lon,p/np.cos(lat)-n
def rot(x):
    a,b,_=llh(x);s,c=np.sin,np.cos
    return np.array([[-s(b),c(b),0],[-s(a)*c(b),-s(a)*s(b),c(a)],[c(a)*c(b),c(a)*s(b),s(a)]])
def frames(path):
    b=path.read_bytes(); i=0
    while True:
        i=b.find(b'\xb5\x62',i)
        if i<0 or i+8>len(b):break
        n=int.from_bytes(b[i+4:i+6],'little');raw=b[i:i+n+8];a=z=0
        for v in raw[2:-2]:a=(a+v)&255;z=(z+a)&255
        if len(raw)==n+8 and raw[-2:]==bytes([a,z]):yield raw
        i+=max(n+8,2)
def read_ubx():
    epochs=[]; receiver=[]; sky=[];dop=[];counts=Counter(); rows={}
    for raw in frames(ROOT/'raw/G08_owndata_20261030.ubx'):
        cl,mid=raw[2:4];p=raw[6:-2];counts[f'{cl:02X}-{mid:02X}']+=1
        if (cl,mid)==(2,21):
            tow,week=struct.unpack_from('<dH',p);n=p[11];ms=[]
            for j in range(n):
                q=p[16+32*j:48+32*j];pr,cp,do=struct.unpack_from('<ddf',q)
                gn,sv,sg,fr=struct.unpack_from('4B',q,20)
                if q[30]&1:ms.append(dict(gn=gn,sv=sv,sig=sg,pr=pr,cno=q[26]))
            epochs.append((week,tow,ms))
        if (cl,mid) in [(1,1),(1,19),(1,7),(1,53),(1,4)]:
            try:
                m=UBXReader.parse(raw);d={k:v for k,v in vars(m).items() if not k.startswith('_')};rows.setdefault(m.identity,[]).append(d)
            except Exception:pass
        if (cl,mid)==(1,19):
            t=struct.unpack_from('<I',p,4)[0]/1000;x=np.array(struct.unpack_from('<iii',p,8))*.01+np.array(struct.unpack_from('<bbb',p,20))*.0001
            if not p[3]&1:receiver.append([t,*x])
        if (cl,mid)==(1,1):
            t,x,y,z=struct.unpack_from('<Iiii',p);receiver.append([t/1000,x*.01,y*.01,z*.01])
        if (cl,mid)==(1,7) and not counts['01-01'] and not counts['01-13']:
            lon,lat,h=struct.unpack_from('<iii',p,24);a,b=np.radians([lat*1e-7,lon*1e-7]);h*=.001;n=6378137/np.sqrt(1-.00669437999014*np.sin(a)**2)
            receiver.append([struct.unpack_from('<I',p)[0]/1000,(n+h)*np.cos(a)*np.cos(b),(n+h)*np.cos(a)*np.sin(b),(n*(1-.00669437999014)+h)*np.sin(a)])
        if (cl,mid)==(1,53):
            t=struct.unpack_from('<I',p)[0]/1000
            for j in range(p[5]):
                gn,sv,cn,el,az=struct.unpack_from('<BBBbh',p,8+j*12);sky.append([t,gn,sv,cn,el,az])
        if (cl,mid)==(1,4):dop.append([struct.unpack_from('<I',p)[0]/1000,*[v*.01 for v in struct.unpack_from('<7H',p,4)]])
        if (cl,mid)==(1,7):dop.append([struct.unpack_from('<I',p)[0]/1000,np.nan,struct.unpack_from('<H',p,76)[0]*.01,np.nan,np.nan,np.nan,np.nan,np.nan])
    for name,rr in rows.items():pd.DataFrame(rr).to_csv(ROOT/f'data/{name}.csv',index=False)
    pd.DataFrame(sky,columns=['tow','gnss','sv','cno','elevation','azimuth']).to_csv(ROOT/'data/sky.csv',index=False)
    (ROOT/'data/message_counts.json').write_text(json.dumps(counts,indent=2))
    return epochs,np.unique(np.array(receiver),axis=0),np.array(sky),np.array(dop)
def navread():
    lines=(ROOT/'data/G08.nav').read_text().splitlines();ion=[]; eph={};i=0
    for l in lines:
        if l.startswith(('GPSA','GPSB')):ion.extend(float(l[j:j+12].replace('D','E')) for j in [5,17,29,41])
        i+=1
        if 'END OF HEADER' in l:break
    while i<len(lines):
        l=lines[i];n=4 if l[0] in 'RS' else 8
        if l[0]=='G':
            a=[float(l[j:j+19].replace('D','E')) for j in [23,42,61]]
            for k in range(1,8):a.extend(float(lines[i+k][j:j+19].replace('D','E') or 0) for j in [4,23,42,61])
            date=datetime(*map(int,l[4:23].split()));toc=(date-datetime(1980,1,6)).total_seconds()%604800
            eph.setdefault(int(l[1:3]),[]).append((toc,np.array(a)))
        i+=n
    return eph,np.array(ion)
def wrap(t):return (t+302400)%604800-302400
def orbit(t,e):
    toc,a=e;af0,af1,af2=a[:3];crs,dn,m0=a[4:7];cuc,ec,cus,sq=a[7:11];toe,cic,O,cis=a[11:15];inc,crc,w,Od=a[15:19];idot=a[19]
    A=sq*sq;dt=wrap(t-toe);M=m0+(np.sqrt(MU/A**3)+dn)*dt;E=M
    for _ in range(12):E=M+ec*np.sin(E)
    v=np.arctan2(np.sqrt(1-ec*ec)*np.sin(E),np.cos(E)-ec);u=v+w;s=np.sin(2*u);c=np.cos(2*u)
    r=A*(1-ec*np.cos(E))+crs*s+crc*c;ii=inc+idot*dt+cis*s+cic*c;u+=cus*s+cuc*c;OO=O+(Od-OM)*dt-OM*toe
    xx=r*np.cos(u);yy=r*np.sin(u);pos=np.array([xx*np.cos(OO)-yy*np.cos(ii)*np.sin(OO),xx*np.sin(OO)+yy*np.cos(ii)*np.cos(OO),yy*np.sin(ii)])
    tc=wrap(t-toc);clk=af0+af1*tc+af2*tc*tc-4.442807633e-10*ec*sq*np.sin(E)
    return pos,clk,a[25]
def atmos(x,sp,t,ion):
    en=(sp-x)@rot(x).T;r=np.linalg.norm(en,axis=1);el=np.arcsin(en[:,2]/r);az=np.arctan2(en[:,0],en[:,1]);lat,lon,h=llh(x);es=np.maximum(el,0.01)
    psi=.0137/(es/np.pi+.11)-.022;phi=np.clip(lat/np.pi+psi*np.cos(az),-.416,.416);lam=lon/np.pi+psi*np.sin(az)/np.cos(phi*np.pi);pm=phi+.064*np.cos((lam-1.617)*np.pi)
    tt=(43200*lam+t)%86400;pol=np.array([np.ones(len(pm)),pm,pm**2,pm**3]);amp=np.maximum(ion[:4]@pol,0);per=np.maximum(ion[4:]@pol,72000);phase=2*np.pi*(tt-50400)/per
    I=C*(1+16*(.53-es/np.pi)**3)*(5e-9+np.where(abs(phase)<1.57,amp*(1-phase**2/2+phase**4/24),0))
    h=np.clip(h,0,10000);pres=1013.25*(1-2.2557e-5*h)**5.2568;temp=288.16-.0065*h;vap=6.108*.7*np.exp((17.15*temp-4684)/(temp-38.45));T=.002277/np.sin(es)*(pres+(1255/temp+.05)*vap)
    return I,T,el
def solve_epoch(t,ms,eph,ion,seed,weighted=True):
    sp=[];pr=[];clk=[];tgd=[];cn=[];ids=[]
    for m in ms:
        if m['gn']!=0 or m['sig']!=0 or m['cno']<15 or m['sv'] not in eph:continue
        e=min(eph[m['sv']],key=lambda e:abs(wrap(t-e[1][11])))
        if abs(wrap(t-e[1][11]))>14400 or e[1][24]!=0:continue
        tx=t-m['pr']/C;_,d,_=orbit(tx,e);p,d,g=orbit(tx-d,e)
        sp.append(p);pr.append(m['pr']);clk.append(d*C);tgd.append(g*C);cn.append(m['cno']);ids.append(m['sv'])
    sp=np.array(sp);pr=np.array(pr);clk=np.array(clk);tgd=np.array(tgd);x=np.r_[seed,0.]
    if len(pr)<4:raise ValueError('fewer than four GPS L1 observations')
    for j in range(30):
        delta=sp-x[:3];r=np.linalg.norm(delta,axis=1);H=np.column_stack([-delta/r[:,None],np.ones(len(r))]);sag=OM/C*(sp[:,0]*x[1]-sp[:,1]*x[0])
        if np.linalg.norm(x[:3])>6e6:I,T,el=atmos(x[:3],sp,t,ion);use=el>=np.radians(15)
        else:I=T=np.zeros(len(r));el=np.ones(len(r));use=np.ones(len(r),bool)
        if use.sum()<4:raise ValueError('fewer than four above mask')
        weights=np.maximum(np.sin(el[use]),.05)**2 if weighted else np.ones(use.sum());v=pr+clk-tgd-I-T-r-sag-x[3];scale=np.sqrt(weights)
        dx,_,rank,_=np.linalg.lstsq(H[use]*scale[:,None],v[use]*scale,rcond=None)
        if rank<4:raise ValueError('rank deficient')
        x+=dx
        if np.linalg.norm(dx)<.001:
            if not 6.30e6<np.linalg.norm(x[:3])<6.42e6:raise ValueError('invalid radius')
            return x, (sp,pr,clk,I,T,tgd,el,use,ids),v[use],H[use]
    raise ValueError('not converged')
def main():
    epochs,ref,sky,dop=read_ubx();eph,ion=navread();seed=np.zeros(3);sol=[];res=[];omitted=[];inputs=[]
    for week,t,ms in epochs:
        try:x,inp,v,H=solve_epoch(t,ms,eph,ion,seed);seed=x[:3];sol.append([week,t,*seed]);inputs.append((t,inp));res.extend([[t,float(r)] for r in v])
        except (ValueError,np.linalg.LinAlgError) as e:omitted.append([week,t,str(e)])
    out=np.array(sol);header=['gps_week','gps_tow_s','x_ecef_m','y_ecef_m','z_ecef_m']
    with (ROOT/'G08_owndata.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(header)
        for a in out:w.writerow([int(a[0]),f'{a[1]:.3f}',*[f'{z:.3f}' for z in a[2:]]])
    pd.DataFrame(out,columns=header).to_csv(ROOT/'data/own_solution_full.csv',index=False);pd.DataFrame(omitted,columns=['week','tow','reason']).to_csv(ROOT/'data/omitted_epochs.csv',index=False)
    pd.DataFrame([[w,t] for w,t,_ in epochs],columns=header[:2]).to_csv(ROOT/'data/epochs.csv',index=False)
    # Preserve the lab matrix layout, GPS PRN rows and epoch columns.
    names=['pseudoranges_meas','satellite_clock_bias','ionospheric_delay','tropospheric_delay','group_delay','elevation_deg'];mats={n:np.full((32,len(epochs)),np.nan) for n in names};S=np.full((32,3*len(epochs)),np.nan);index={round(t,6):k for k,(_,t,_) in enumerate(epochs)}
    for t,(sp,pr,cl,I,T,g,el,use,ids) in inputs:
        k=index[round(t,6)]
        for j,sv in enumerate(ids):
            if not use[j]:continue
            for n,val in zip(names,[pr[j],cl[j],I[j],T[j],g[j],np.degrees(el[j])]):mats[n][sv-1,k]=val
            S[sv-1,3*k:3*k+3]=sp[j]
    for n,a in mats.items():np.savetxt(ROOT/f'data/{n}.csv',a,delimiter=',')
    np.savetxt(ROOT/'data/satellite_positions.csv',S,delimiter=',');pd.DataFrame({'satellite_id':[f'G{i:02}' for i in range(1,33)]}).to_csv(ROOT/'data/satellites.csv',index=False)
    origin=np.mean(ref[:,1:],axis=0);R=rot(origin);enu=(out[:,2:]-origin)@R.T;reference={round(r[0],1):r[1:] for r in ref};rtk={}
    for l in (ROOT/'data/RTKLIB_GPS.pos').read_text().splitlines():
        if l.startswith('%') or not l.strip():continue
        a=l.split();rtk[round(float(a[1]),1)]=np.array(list(map(float,a[2:5])))
    stats=[]
    def compare(name,left,right):
        d=np.array([(p-right[t])@R.T for t,p in left.items() if t in right]);
        if not len(d):return
        stats.append(dict(comparison=name,n=len(d),**{f'mean_{a}_m':float(v) for a,v in zip('ENU',d.mean(0))},**{f'std_{a}_m':float(v) for a,v in zip('ENU',d.std(0,ddof=1))},rmse2d_m=float(np.sqrt(np.mean(np.sum(d[:,:2]**2,1)))),rmse3d_m=float(np.sqrt(np.mean(np.sum(d**2,1))))))
    own={round(a[1],1):a[2:] for a in out};compare('Custom GPS WLS - receiver',own,reference);compare('Custom GPS WLS - RTKLIB GPS',own,rtk);compare('RTKLIB GPS - receiver',rtk,reference);pd.DataFrame(stats).to_csv(ROOT/'data/comparison_statistics.csv',index=False)
    fig,ax=plt.subplots(figsize=(6.5,5));
    for name,vals in [('Custom GPS WLS',out[:,2:]),('Receiver',ref[:,1:]),('RTKLIB GPS',np.array(list(rtk.values())))]:
        en=(vals-origin)@R.T;ax.scatter(en[:,0],en[:,1],s=5,label=name,alpha=.5)
    ax.set(xlabel='East relative to receiver mean (m)',ylabel='North relative to receiver mean (m)',title='Own-data solution comparison');ax.axis('equal');ax.grid(alpha=.3);ax.legend();fig.tight_layout();fig.savefig(ROOT/'figures/scatter.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(8,7),sharex=True);tt=out[:,1]-out[0,1]
    for j,a in enumerate(axes):a.plot(tt,enu[:,j]);a.set_ylabel('ENU'[j]+' (m)');a.grid(alpha=.3)
    axes[-1].set_xlabel('Seconds from first observation');fig.tight_layout();fig.savefig(ROOT/'figures/enu_time.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(9,6));df=pd.DataFrame(sky,columns=['tow','gnss','sv','cno','el','az']);labels={0:'GPS',1:'SBAS',2:'Galileo',3:'BeiDou',5:'QZSS',6:'GLONASS'}
    for gn,g in df[df.cno>0].groupby('gnss'):
        counts=g.groupby('tow').sv.nunique();axes[0,0].plot(counts.index-out[0,1],counts.values,label=labels.get(gn,str(gn)));axes[0,1].hist(g.cno,bins=25,alpha=.4,label=labels.get(gn,str(gn)))
    axes[0,0].set(title='Tracked satellites per system',ylabel='Satellite count',xlabel='Elapsed seconds');axes[0,0].legend(fontsize=7,ncol=2);axes[0,1].set(title='Tracked-signal C/N0',xlabel='dB-Hz',ylabel='Count');
    if len(dop):axes[1,0].plot(dop[:,0]-out[0,1],dop[:,2],label='Receiver PDOP (NAV-PVT)');axes[1,0].legend()
    axes[1,0].set(title='Receiver DOP',xlabel='Elapsed seconds');axes[1,1].hist(np.array(res)[:,1],bins=50);axes[1,1].set(title='Custom WLS post-fit residuals',xlabel='Metres');fig.tight_layout();fig.savefig(ROOT/'figures/diagnostics.png',dpi=180);plt.close(fig)
    fig=plt.figure(figsize=(6,6));ax=fig.add_subplot(projection='polar');g=df[(df.cno>0)&(df.el>=0)];s=ax.scatter(np.radians(g.az),90-g.el,c=g.cno,s=2,cmap='viridis');ax.set_theta_zero_location('N');ax.set_theta_direction(-1);ax.set_ylim(0,90);ax.set_yticks([0,30,60,90],['90','60','30','0']);ax.set_title('Receiver sky plot - elevation rings (deg)');fig.colorbar(s,ax=ax,label='C/N0 (dB-Hz)',shrink=.7);fig.savefig(ROOT/'figures/skyplot.png',dpi=180);plt.close(fig)
    centered=enu-enu.mean(0);summary=dict(observation_epochs=len(epochs),solved_epochs=len(out),omitted=len(omitted),week=int(out[0,0]),first_tow=float(epochs[0][1]),last_tow=float(epochs[-1][1]),duration_s=float(epochs[-1][1]-epochs[0][1]),receiver_epochs=len(ref),receiver_message='NAV-HPPOSECEF' if (ROOT/'data/NAV-HPPOSECEF.csv').exists() else ('NAV-POSECEF' if (ROOT/'data/NAV-POSECEF.csv').exists() else 'NAV-PVT LLH converted to ECEF'),mean_llh=[float(np.degrees(v)) for v in llh(origin)[:2]]+[float(llh(origin)[2])],trajectory_std_enu_m=list(centered.std(0,ddof=1)),trajectory_horizontal_p95_about_mean_m=float(np.percentile(np.linalg.norm(centered[:,:2],axis=1),95)),comparisons=stats)
    (ROOT/'data/summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
