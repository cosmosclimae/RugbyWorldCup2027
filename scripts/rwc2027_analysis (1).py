#!/usr/bin/env python3
# =====================================================================
#  RWC 2027 — ENSO/IOD-conditioned match-weather typology
#  FULL REPRODUCIBLE ANALYSIS (classification + typology + 7 figures + CIs)
#
#  INPUT  (place in ./ or edit UP):  rugbywc2027_event_samples_master.csv
#  OUTPUT (written to OUT):
#     season_classification_1996_2025.csv
#     fig1_map_schedule.png ... fig7_sensitivity.png
#     typology_final_by_phase.csv
#
#  Run:  python rwc2027_analysis.py
#  Deps: numpy, pandas, matplotlib
#
#  Definitions (documented for verification)
#   - ENSO phase : ONI on OND season ; El Nino>=+0.5, La Nina<=-0.5, else neutral
#   - IOD  phase : DMI(SON) ; Positive>=+0.4, Negative<=-0.4, else neutral
#   - Families   : Thermal WBGT_max>=26C ; Wet (rain_window>=1 OR rain_24h>=5
#                  OR rain_5d>=25 mm) ; Wind gust_max>=12 m/s
#   - Types      : Compound if >=2 families active, else single family,
#                  else Fair (mutually exclusive)
#   - Estimate   : per-day frequency over all +-2d realisations of a phase
#   - Uncertainty: block-bootstrap resampling whole YEARS (proper unit)
# =====================================================================
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch

UP="/mnt/user-data/uploads"; OUT="/mnt/user-data/outputs"
rng=np.random.default_rng(0)

# ---------------------------------------------------------------- 0. INDICES
# NOAA CPC ONI (Year, DJF..NDJ) ; NOAA PSL DMI HadISST long (Year, Jan..Dec)
ONI="""1996 -0.9 -0.8 -0.6 -0.4 -0.3 -0.3 -0.3 -0.3 -0.4 -0.4 -0.4 -0.5
1997 -0.5 -0.4 -0.1 0.3 0.8 1.2 1.6 1.9 2.1 2.3 2.4 2.4
1998 2.2 1.9 1.4 1.0 0.5 -0.1 -0.8 -1.1 -1.3 -1.4 -1.5 -1.6
1999 -1.5 -1.3 -1.1 -1.0 -1.0 -1.0 -1.1 -1.1 -1.2 -1.3 -1.5 -1.7
2000 -1.7 -1.4 -1.1 -0.8 -0.7 -0.6 -0.6 -0.5 -0.5 -0.6 -0.7 -0.7
2001 -0.7 -0.5 -0.4 -0.3 -0.3 -0.1 -0.1 -0.1 -0.2 -0.3 -0.3 -0.3
2002 -0.1 0.0 0.1 0.2 0.4 0.7 0.8 0.9 1.0 1.2 1.3 1.1
2003 0.9 0.6 0.4 0.0 -0.3 -0.2 0.1 0.2 0.3 0.3 0.4 0.4
2004 0.4 0.3 0.2 0.2 0.2 0.3 0.5 0.6 0.7 0.7 0.7 0.7
2005 0.6 0.6 0.4 0.4 0.3 0.1 -0.1 -0.1 -0.1 -0.3 -0.6 -0.8
2006 -0.9 -0.8 -0.6 -0.4 -0.1 0.0 0.1 0.3 0.5 0.8 0.9 0.9
2007 0.7 0.2 -0.1 -0.3 -0.4 -0.5 -0.6 -0.8 -1.1 -1.3 -1.5 -1.6
2008 -1.6 -1.5 -1.3 -1.0 -0.8 -0.6 -0.4 -0.2 -0.2 -0.4 -0.6 -0.7
2009 -0.8 -0.8 -0.6 -0.3 0.0 0.3 0.5 0.6 0.7 1.0 1.4 1.6
2010 1.5 1.2 0.8 0.4 -0.2 -0.7 -1.0 -1.3 -1.6 -1.6 -1.6 -1.5
2011 -1.3 -1.0 -0.8 -0.6 -0.5 -0.4 -0.4 -0.6 -0.8 -1.0 -1.0 -0.9
2012 -0.7 -0.6 -0.5 -0.4 -0.2 0.1 0.3 0.4 0.4 0.3 0.1 -0.1
2013 -0.3 -0.3 -0.2 -0.2 -0.3 -0.3 -0.4 -0.3 -0.2 -0.1 -0.1 -0.2
2014 -0.3 -0.3 -0.1 0.2 0.3 0.2 0.1 0.1 0.3 0.5 0.7 0.8
2015 0.7 0.6 0.7 0.8 1.0 1.3 1.6 1.9 2.2 2.5 2.6 2.8
2016 2.6 2.3 1.7 1.0 0.5 0.0 -0.3 -0.5 -0.6 -0.6 -0.6 -0.5
2017 -0.2 0.0 0.2 0.3 0.4 0.4 0.2 -0.1 -0.3 -0.6 -0.8 -0.9
2018 -0.8 -0.7 -0.6 -0.4 -0.1 0.1 0.1 0.3 0.5 0.8 1.0 0.9
2019 0.9 0.9 0.8 0.8 0.6 0.5 0.3 0.2 0.2 0.4 0.6 0.7
2020 0.6 0.6 0.5 0.3 0.0 -0.2 -0.4 -0.5 -0.8 -1.1 -1.2 -1.1
2021 -0.9 -0.8 -0.7 -0.5 -0.4 -0.3 -0.3 -0.4 -0.6 -0.8 -0.9 -0.9
2022 -0.8 -0.8 -0.9 -1.0 -0.9 -0.8 -0.8 -0.9 -1.0 -0.9 -0.8 -0.7
2023 -0.5 -0.3 0.0 0.3 0.6 0.8 1.1 1.4 1.6 1.8 2.0 2.1
2024 1.9 1.6 1.3 0.8 0.5 0.2 0.1 -0.1 -0.2 -0.2 -0.3 -0.4
2025 -0.4 -0.2 -0.1 0.0 0.0 0.0 -0.1 -0.3 -0.4 -0.5 -0.6 -0.5"""
DMI="""1996 -0.021 -0.033 -0.085 -0.369 -0.266 -0.394 -0.643 -0.681 -0.712 -1.108 -0.797 -0.413
1997 -0.110 0.079 0.043 0.054 0.025 0.082 0.447 0.634 0.771 0.873 1.279 0.863
1998 0.525 0.422 -0.055 0.047 0.140 0.147 -0.385 -0.580 -0.496 -0.743 -0.653 -0.336
1999 -0.130 -0.038 0.102 -0.002 -0.186 -0.145 0.112 0.023 -0.050 -0.190 -0.116 -0.148
2000 -0.125 -0.009 0.140 0.141 0.113 -0.022 0.077 0.131 -0.096 -0.142 -0.306 -0.248
2001 -0.431 -0.017 -0.010 0.136 0.137 0.127 -0.140 -0.301 -0.223 -0.451 -0.245 -0.051
2002 -0.142 -0.098 0.008 -0.352 -0.329 -0.190 -0.260 -0.208 0.286 0.405 0.096 -0.158
2003 -0.239 0.017 -0.039 -0.099 -0.164 0.139 0.131 0.108 -0.061 -0.243 -0.159 0.189
2004 0.044 0.114 0.077 -0.080 -0.565 -0.384 -0.301 -0.132 -0.106 -0.004 -0.153 -0.132
2005 -0.124 -0.560 -0.432 0.084 -0.033 -0.176 -0.384 -0.345 -0.534 -0.436 -0.272 -0.300
2006 -0.135 -0.305 -0.226 -0.015 -0.207 -0.059 0.040 0.229 0.428 0.577 0.501 0.172
2007 0.224 0.150 0.116 0.100 0.270 0.030 0.045 0.236 0.235 0.068 -0.023 -0.227
2008 0.115 -0.072 0.097 -0.153 0.178 0.215 0.246 0.124 0.086 0.029 -0.128 -0.053
2009 0.031 0.163 0.100 0.126 0.256 0.100 -0.192 -0.104 -0.103 -0.013 -0.067 0.160
2010 0.294 0.023 0.458 0.370 -0.030 -0.140 -0.001 -0.062 -0.268 -0.437 -0.495 -0.212
2011 0.192 0.242 0.367 0.154 -0.089 0.055 0.223 0.344 0.202 0.356 0.336 -0.128
2012 0.046 -0.078 0.026 -0.271 -0.370 0.001 0.546 0.652 0.453 0.110 -0.100 0.268
2013 -0.065 0.189 0.083 -0.297 -0.506 -0.497 -0.180 -0.194 -0.310 -0.168 0.199 0.141
2014 -0.101 -0.089 -0.151 -0.058 -0.092 -0.028 -0.363 -0.372 -0.145 0.141 0.010 0.046
2015 -0.100 -0.345 -0.241 -0.008 0.240 0.296 0.225 0.567 0.294 0.483 0.347 0.272
2016 0.266 -0.110 -0.009 0.146 -0.113 -0.443 -0.758 -0.444 -0.437 -0.372 -0.382 -0.310
2017 -0.086 0.101 0.357 0.499 0.536 0.424 0.520 0.349 0.034 0.016 0.289 0.109
2018 -0.200 0.215 -0.120 -0.083 0.122 0.155 0.053 0.122 0.604 0.685 0.500 0.309
2019 0.387 0.416 0.224 0.258 0.539 0.605 0.597 0.436 0.893 0.964 0.835 0.243
2020 0.173 0.054 0.019 -0.011 0.298 0.454 0.320 -0.183 -0.190 0.074 0.020 0.030
2021 0.051 0.243 0.266 0.250 0.009 -0.002 -0.228 -0.099 -0.058 -0.091 0.069 -0.120
2022 -0.056 -0.083 -0.093 -0.068 -0.122 -0.335 -0.195 -0.246 -0.322 -0.691 -0.269 -0.092
2023 0.109 0.157 0.415 0.560 0.443 0.665 0.498 0.825 0.946 0.804 0.920 0.851
2024 0.765 0.328 0.421 0.440 0.297 0.197 0.033 0.267 0.115 -0.196 -0.383 -0.331"""
def _parse(b):
    return pd.DataFrame([[int(p[0])]+[float(x) for x in p[1:13]] for p in
                         (ln.split() for ln in b.strip().splitlines())])
oni=_parse(ONI); oni.columns=["year","DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ"]
dmi=_parse(DMI); dmi.columns=["year"]+list(range(1,13))
cls=pd.DataFrame({"year":oni.year,"ONI_OND":oni.OND.values})
cls["DMI_SON"]=dmi.set_index("year").reindex(cls.year)[[9,10,11]].mean(axis=1).values
cls["ENSO"]=cls.ONI_OND.map(lambda v:"El Nino" if v>=.5 else "La Nina" if v<=-.5 else "Neutral")
cls["IOD"] =cls.DMI_SON.map(lambda v:np.nan if pd.isna(v) else "Positive" if v>=.4 else "Negative" if v<=-.4 else "Neutral")
cls.to_csv(f"{OUT}/season_classification_1996_2025.csv",index=False)

# ---------------------------------------------------------------- 1. TYPOLOGY
m=pd.read_csv(f"{UP}/rugbywc2027_event_samples_master.csv").merge(cls[["year","ENSO","IOD"]],on="year",how="left")
m["match_date"]=pd.to_datetime(m["match_date"])
th=m.WBGT_max>=26
wet=(m.rain_window>=1)|(m.rain_24h_before>=5)|(m.rain_5d_before>=25)
wd=m.gust_max>=12
K=th.astype(int)+wet+wd
m["type"]=np.select([K>=2,th,wet,wd],["Compound","Hot","Wet","Windy"],default="Fair")
m["moist"]=m.type.isin(["Wet","Compound"])
m["wet_trig"]=wet   # H_W: at least one rainfall trigger active
# ---- Wet-family decomposition (descriptive supplement, Table S1) ----
def wet_decomp(mode,order):
    r=[]
    for ph in order:
        d=m[m[mode]==ph]
        if len(d)==0: continue
        r.append({"Mode":mode,"Phase":ph,"n_years":int(d.year.nunique()),"n_realisations":len(d),
            "pct_Pmatch_ge1":(d.rain_window>=1).mean()*100,
            "pct_P24h_ge5":(d.rain_24h_before>=5).mean()*100,
            "pct_P5d_ge25":(d.rain_5d_before>=25).mean()*100,
            "pct_any_wet":d.wet_trig.mean()*100,
            "pct_Wet_only":(d.type=="Wet").mean()*100,
            "pct_Compound_wet":((d.type=="Compound")&(d.wet_trig)).mean()*100,
            "pct_Compound_nowet":((d.type=="Compound")&(~d.wet_trig)).mean()*100})
    return pd.DataFrame(r)
wet_decomposition=pd.concat([wet_decomp("ENSO",["La Nina","Neutral","El Nino"]),
                             wet_decomp("IOD",["Negative","Neutral","Positive"])],ignore_index=True)
wet_decomposition.round(2).to_csv(f"{OUT}/wet_family_decomposition_by_phase.csv",index=False)
TYPES=["Fair","Windy","Wet","Hot","Compound"]
COL={"Fair":"#bcd9b0","Windy":"#7b5ea7","Wet":"#2b6cb0","Hot":"#e07a3f","Compound":"#2d2d2d"}
legend=[Patch(facecolor=COL[t],label=t) for t in TYPES]
def comp(df,by,order): return (df.groupby(by).type.value_counts(normalize=True).unstack(fill_value=0)*100).reindex(index=order,columns=TYPES,fill_value=0)
def stacked(ax,tab,rot=0):
    b=np.zeros(len(tab))
    for t in TYPES:
        ax.bar(tab.index,tab[t],bottom=b,color=COL[t],edgecolor="white",lw=.5)
        for i,(bb,v) in enumerate(zip(b,tab[t])):
            if v>=3.5: ax.text(i,bb+v/2,f"{v:.0f}",ha="center",va="center",fontsize=6.5,color="white" if t in("Wet","Compound","Hot") else "#333")
        b+=tab[t].values
    ax.set_ylim(0,100); ax.tick_params(axis="x",rotation=rot)
mt=m.groupby("stadium").agg(city=("city","first"),lat=("lat","first"),lon=("lon","first"),n=("match_id","nunique")).reset_index()
enso=comp(m.dropna(subset=["ENSO"]),"ENSO",["La Nina","Neutral","El Nino"])
iod =comp(m.dropna(subset=["IOD"]),"IOD",["Negative","Neutral","Positive"])
pd.concat({"ENSO":enso,"IOD":iod}).round(2).to_csv(f"{OUT}/typology_final_by_phase.csv")

# ---------------------------------------------------------------- FIG 1 (map)
AUS=[(113.4,-22.0),(113.8,-26.1),(114.6,-28.8),(115.0,-30.0),(115.7,-31.6),(115.0,-33.6),(115.1,-34.4),(118.0,-35.1),(120.5,-33.9),(123.6,-33.9),(126.1,-32.3),(129.0,-31.7),(131.5,-31.5),(134.2,-32.6),(135.2,-34.6),(136.8,-35.3),(137.6,-35.0),(138.0,-35.6),(139.8,-37.0),(140.9,-38.1),(143.6,-38.8),(146.3,-38.7),(146.4,-39.1),(148.4,-37.8),(149.9,-37.5),(150.8,-35.1),(151.3,-33.9),(152.0,-32.4),(153.0,-30.3),(153.6,-28.6),(153.1,-25.9),(151.0,-24.0),(149.1,-21.0),(146.3,-19.0),(145.5,-16.8),(143.9,-14.5),(142.6,-10.9),(141.7,-12.7),(140.8,-17.4),(139.4,-17.4),(137.0,-16.0),(136.7,-12.2),(135.0,-14.9),(132.6,-11.5),(130.6,-12.5),(129.0,-14.8),(126.5,-14.0),(122.2,-18.1),(120.0,-19.8),(116.0,-20.6),(113.4,-22.0)]
TAS=[(145.5,-40.7),(148.3,-40.9),(148.3,-42.0),(147.9,-43.6),(146.0,-43.5),(145.2,-42.2),(144.7,-40.9),(145.5,-40.7)]
off={"stadium_australia":(2.0,0.6),"sydney_football_stadium":(2.4,-1.2),"newcastle_stadium":(2.2,0.9)}
fig,ax=plt.subplots(1,3,figsize=(14,4.6),gridspec_kw={"width_ratios":[1.5,1,1]})
ax[0].add_patch(Polygon(AUS,closed=True,facecolor="#eef2e9",edgecolor="#9bb069",lw=1))
ax[0].add_patch(Polygon(TAS,closed=True,facecolor="#eef2e9",edgecolor="#9bb069",lw=1))
ax[0].scatter(mt.lon,mt.lat,s=40+90*mt.n,c="#2b6cb0",alpha=.85,edgecolor="white",zorder=4)
for _,r in mt.iterrows():
    dx,dy=off.get(r.stadium,(0,1.1))
    ax[0].annotate(f"{r.city} (n={r.n})",(r.lon,r.lat),(r.lon+dx,r.lat+dy),fontsize=7,ha="left" if dx else "center",va="center",zorder=5,arrowprops=dict(arrowstyle="-",lw=.5,color="#777") if dx else None)
ax[0].set_xlim(110,156); ax[0].set_ylim(-45,-9); ax[0].set_aspect(1.15)
ax[0].set_xlabel("Longitude"); ax[0].set_ylabel("Latitude"); ax[0].set_title("Host venues (bubble = n matches)",fontsize=9)
ax[1].hist(m.groupby("match_id").kickoff_hour_local.first(),bins=np.arange(11,23,1),color="#7b5ea7",edgecolor="white")
ax[1].set_title("Local kick-off hour",fontsize=9); ax[1].set_xlabel("hour"); ax[1].set_ylabel("matches")
md=m.groupby("match_id").agg(d=("match_date","first"),ph=("phase","first")); pc={"pool":"#2b6cb0","round_of_16":"#5aa0d0","quarter_final":"#e07a3f","semi_final":"#c0392b","bronze_final":"#888","final":"#2d2d2d"}
for _,r in md.iterrows(): ax[2].scatter(r.d,0,color=pc.get(r.ph,"#999"),s=16)
ax[2].set_yticks([]); ax[2].set_title("Tournament calendar",fontsize=9); ax[2].tick_params(axis="x",rotation=30)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_map_schedule.png",dpi=150)

# ---------------------------------------------------------------- FIG 2 (classification)
fig,ax=plt.subplots(2,1,figsize=(11,5),sharex=True)
ce={"La Nina":"#2b6cb0","Neutral":"#bbb","El Nino":"#c0392b"}; ci={"Negative":"#2b6cb0","Neutral":"#bbb","Positive":"#c0392b"}
ax[0].bar(cls.year,cls.ONI_OND,color=[ce[e] for e in cls.ENSO]); ax[0].axhline(.5,ls="--",c="#c0392b",lw=.8); ax[0].axhline(-.5,ls="--",c="#2b6cb0",lw=.8)
ax[0].set_ylabel("ONI (OND)"); ax[0].set_title("ENSO phase of each Oct-Nov season",fontsize=9)
ax[1].bar(cls.year,cls.DMI_SON,color=[ci.get(i,"#ddd") for i in cls.IOD.fillna("Neutral")]); ax[1].axhline(.4,ls="--",c="#c0392b",lw=.8); ax[1].axhline(-.4,ls="--",c="#2b6cb0",lw=.8)
ax[1].set_ylabel("DMI (SON)"); ax[1].set_title("IOD phase of each Oct-Nov season",fontsize=9); ax[1].set_xlabel("year")
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_phase_classification.png",dpi=150)

# ---------------------------------------------------------------- FIG 3 (baseline)
m["kcat"]=pd.cut(m["kickoff_hour_local"],[0,14,17.5,24],labels=["Midday","Afternoon","Evening"])
overall=(m.type.value_counts(normalize=True)*100).reindex(TYPES).fillna(0)
fig,ax=plt.subplots(1,3,figsize=(14,4.4),gridspec_kw={"width_ratios":[.6,1.6,1]})
b=0
for t in TYPES:
    ax[0].bar(["All"],overall[t],bottom=b,color=COL[t])
    if overall[t]>=3.5: ax[0].text(0,b+overall[t]/2,f"{overall[t]:.0f}",ha="center",va="center",fontsize=7,color="white" if t in("Wet","Compound","Hot") else "#333")
    b+=overall[t]
ax[0].set_ylim(0,100); ax[0].set_title("All 52 matches",fontsize=9); ax[0].set_ylabel("% realizations")
stacked(ax[1],comp(m,"city",sorted(m.city.unique())),rot=35); ax[1].set_title("By host city",fontsize=9)
stacked(ax[2],comp(m,"kcat",["Midday","Afternoon","Evening"])); ax[2].set_title("By kick-off period",fontsize=9)
fig.legend(handles=legend,ncol=5,fontsize=8,loc="lower center",bbox_to_anchor=(.5,-.02))
fig.suptitle("Baseline match-type composition (all phases)",fontsize=10); fig.tight_layout(rect=[0,0.05,1,0.96]); fig.savefig(f"{OUT}/fig3_baseline_typology.png",dpi=150)

# ---------------------------------------------------------------- FIG 4 (headline)
fig,ax=plt.subplots(1,2,figsize=(10,4.4))
stacked(ax[0],enso); ax[0].set_title("by ENSO phase",fontsize=9); ax[0].set_ylabel("% of match realizations")
stacked(ax[1],iod);  ax[1].set_title("by IOD phase",fontsize=9)
fig.legend(handles=legend,ncol=5,fontsize=8,loc="lower center",bbox_to_anchor=(.5,-.01))
fig.suptitle("RWC 2027 — match-type composition by climate-mode phase",fontsize=10); fig.tight_layout(rect=[0,0.05,1,0.95]); fig.savefig(f"{OUT}/fig4_typology_by_phase.png",dpi=150)

# ---------------------------------------------------------------- FIG 5 (north/south)
def share(mode,ph): s=m[m[mode]==ph]; return s.groupby("stadium").moist.mean()*100
lat=m.groupby("stadium").lat.first(); city=m.groupby("stadium").city.first()
g=pd.DataFrame({"lat":lat,"city":city,"ENSO":share("ENSO","La Nina")-share("ENSO","El Nino"),"IOD":share("IOD","Negative")-share("IOD","Positive")}).sort_values("lat")
fig,ax=plt.subplots(figsize=(8,4.6)); y=np.arange(len(g))
ax.hlines(y,np.minimum(g.ENSO,g.IOD),np.maximum(g.ENSO,g.IOD),color="#ddd",lw=2)
ax.scatter(g.ENSO,y,color="#c0392b",label="ENSO (La Nina - El Nino)",zorder=3); ax.scatter(g.IOD,y,color="#2b6cb0",label="IOD (Negative - Positive)",zorder=3)
ax.axvline(0,color="#999",lw=.8); ax.set_yticks(y); ax.set_yticklabels([f"{c} ({l:.0f})" for c,l in zip(g.city,g.lat)],fontsize=8)
ax.set_xlabel("Wet+Compound sensitivity (pp)"); ax.set_title("Moisture sensitivity to ENSO vs IOD, venues S->N",fontsize=9); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/fig5_north_south.png",dpi=150)

# ---------------------------------------------------------------- FIG 6 (operational + bootstrap CI)
yr_rate=m.groupby("year").moist.mean()*100; pe=dict(zip(cls.year,cls.ENSO)); pi=dict(zip(cls.year,cls.IOD))
def boot(years,n=4000):
    v=yr_rate.loc[[y for y in years if y in yr_rate.index]].values
    s=[rng.choice(v,len(v),True).mean() for _ in range(n)]; return v.mean(),np.percentile(s,2.5),np.percentile(s,97.5)
rows=[("La Nina",[y for y,p in pe.items() if p=="La Nina"],"#2b6cb0"),("Neutral",[y for y,p in pe.items() if p=="Neutral"],"#bbb"),("El Nino",[y for y,p in pe.items() if p=="El Nino"],"#c0392b"),
      ("IOD-",[y for y,p in pi.items() if p=="Negative"],"#1f77b4"),("IOD0",[y for y,p in pi.items() if p=="Neutral"],"#999"),("IOD+",[y for y,p in pi.items() if p=="Positive"],"#d62728")]
labs,pts,lo,hi,cc=[],[],[],[],[]
for lab,ys,c in rows: p,a,b=boot(ys); labs.append(lab);pts.append(p);lo.append(p-a);hi.append(b-p);cc.append(c)
fig,ax=plt.subplots(1,2,figsize=(13,4.6),gridspec_kw={"width_ratios":[1,1.3]})
ax[0].bar(labs,pts,yerr=[lo,hi],color=cc,capsize=4,edgecolor="white"); ax[0].set_ylabel("Wet+Compound share (%)"); ax[0].set_title("Tournament moisture by phase (95% CI, year-bootstrap)",fontsize=9); ax[0].tick_params(axis="x",rotation=20)
hv=(m.dropna(subset=["ENSO"]).groupby(["stadium","ENSO"]).moist.mean()*100).unstack().reindex(columns=["La Nina","Neutral","El Nino"]).reindex(lat.sort_values().index)
im=ax[1].imshow(hv.values,cmap="Blues",aspect="auto",vmin=0); ax[1].set_xticks(range(3)); ax[1].set_xticklabels(["La Nina","Neutral","El Nino"]); ax[1].set_yticks(range(len(hv))); ax[1].set_yticklabels([city[s] for s in hv.index],fontsize=8)
for i in range(hv.shape[0]):
    for j in range(hv.shape[1]): ax[1].text(j,i,f"{hv.values[i,j]:.0f}",ha="center",va="center",fontsize=7,color="white" if hv.values[i,j]>hv.values.max()*.6 else "#333")
ax[1].set_title("Wet+Compound % by venue x ENSO phase",fontsize=9); fig.colorbar(im,ax=ax[1],fraction=.046)
fig.tight_layout(); fig.savefig(f"{OUT}/fig6_operational.png",dpi=150)

# ---------------------------------------------------------------- FIG 7 (sensitivity)
def contrast(wbgt=26,rm=1,r24=5,r5=25,gst=12):
    t=m.WBGT_max>=wbgt; w=(m.rain_window>=rm)|(m.rain_24h_before>=r24)|(m.rain_5d_before>=r5); g2=m.gust_max>=gst
    kk=t.astype(int)+w+g2; ty=np.select([kk>=2,t,w,g2],["Compound","Hot","Wet","Windy"],default="Fair")
    d=m.assign(mo=np.isin(ty,["Wet","Compound"]))
    def sh(col,a,b): gg=d.dropna(subset=[col]).groupby(col).mo.mean()*100; return gg.get(a,np.nan)-gg.get(b,np.nan)
    return sh("ENSO","La Nina","El Nino"), sh("IOD","Negative","Positive")
grids={"rain_match (mm)":("rm",[0.5,1,2,5]),"rain_24h (mm)":("r24",[3,5,10,15]),"rain_5d (mm)":("r5",[15,25,35,50]),"gust (m/s)":("gst",[10,12,15,18]),"WBGT (C)":("wbgt",[24,26,28,30])}
base=dict(wbgt=26,rm=1,r24=5,r5=25,gst=12)
fig,ax=plt.subplots(1,5,figsize=(16,3.4),sharey=True)
for a,(title,(key,vals)) in zip(ax,grids.items()):
    e=[];ii=[]
    for v in vals: kw=dict(base);kw[key]=v;x,z=contrast(**kw);e.append(x);ii.append(z)
    a.plot(vals,e,"o-",color="#c0392b",label="ENSO"); a.plot(vals,ii,"s-",color="#2b6cb0",label="IOD")
    a.axhline(0,color="#999",lw=.8); a.axvline(base[key],color="#bbb",ls=":",lw=1); a.set_title(title,fontsize=9); a.set_xlabel("threshold")
ax[0].set_ylabel("Wet+Compound contrast (pp)"); ax[0].legend(fontsize=7,loc="lower right")
fig.suptitle("Sensitivity of the moisture contrast to thresholds (baseline dotted)",fontsize=10); fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig(f"{OUT}/fig7_sensitivity.png",dpi=150)

# ---------------------------------------------------------------- print key numbers
print("ENSO season counts:",cls.ENSO.value_counts().to_dict())
print("\nTypology by ENSO (%):\n",enso.round(1).to_string())
print("\nTypology by IOD (%):\n",iod.round(1).to_string())
print("\nMoisture (Wet+Compound) by phase, point [95% CI year-bootstrap]:")
for lab,ys,_ in rows: p,a,b=boot(ys); print(f"  {lab:8s} n={len([y for y in ys if y in yr_rate.index]):2d}  {p:5.1f}  [{a:4.1f},{b:4.1f}]")
def boot_contrast(ya,yb,n=10000):
    va=yr_rate.loc[[y for y in ya if y in yr_rate.index]].values
    vb=yr_rate.loc[[y for y in yb if y in yr_rate.index]].values
    s=np.array([rng.choice(va,len(va),True).mean()-rng.choice(vb,len(vb),True).mean() for _ in range(n)])
    return va.mean()-vb.mean(), np.percentile(s,2.5), np.percentile(s,97.5)
ln=[y for y,p in pe.items() if p=="La Nina"]; en=[y for y,p in pe.items() if p=="El Nino"]
ng=[y for y,p in pi.items() if p=="Negative"]; ps=[y for y,p in pi.items() if p=="Positive"]
ec=boot_contrast(ln,en); ic=boot_contrast(ng,ps)
print("\nMoisture (Wet+Compound) CONTRAST -- single source of truth [95% CI, year-bootstrap]:")
print(f"  ENSO  La Nina - El Nino  : {ec[0]:+.1f} pp  [{ec[1]:+.1f}, {ec[2]:+.1f}]")
print(f"  IOD   Negative - Positive: {ic[0]:+.1f} pp  [{ic[1]:+.1f}, {ic[2]:+.1f}]")
print("baseline contrast (ENSO,IOD):",[round(x,1) for x in contrast(**base)])
print("\nAll figures + CSVs written to",OUT)
