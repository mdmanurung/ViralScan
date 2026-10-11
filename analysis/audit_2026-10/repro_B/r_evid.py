from viralscan import evidence as e
from viralscan import molecule_verdict as mv
print("gini([100])",e._gini([100]),"gini([50,50])",e._gini([50,50]),"gini([1,99])",e._gini([1,99]))
print("target_assignment same-virus overlapping genes:",e._target_assignment([1,2],{1,2},{1,2},"em-global"))
rows=[{"reference_class":"virus","max_50bp_window_fraction":1.0,"breadth_1x":0.01,"reads":1,"host_competitive_fraction":0.0}]
for r in e.interpretation_flags(rows,[],[]): 
    if r["flag"] in("coverage_hotspot","eve_or_integration_like_hotspot"): print(r["flag"],r["status"],r["metric"])
lin=[{"ambiguity_class":"virus_virus_ambiguous"}]*3+[{"ambiguity_class":"unique_target"}]
print([ (r["flag"],r["status"]) for r in e.interpretation_flags(rows,[],lin) if r["flag"]=="sibling_or_host_ambiguity"])
# blank bitscore counted as host-preferred
print([ (r["flag"],r["status"],r["metric"]) for r in e.interpretation_flags(rows,[{"viral_minus_host_bitscore":""}],[]) if r["flag"]=="host_homology"])
# window grid split
