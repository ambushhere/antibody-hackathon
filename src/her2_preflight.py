"""Conservative local screen, not an abguard or Adaptyv database replacement."""
import argparse
import json
import pandas as pd
from immunum import Annotator
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
from her2_candidates import VH_PREFIX, VH_SUFFIX, WT_CDRH3


def annotate(vh, vl):
    result = {'numbering': {}, 'cdrh3_by_scheme': {}}
    for scheme in ['imgt','kabat','martin','chothia']:
        a = Annotator(chains=['H','K','L'], scheme=scheme, min_confidence=0.1)
        n = a.number(vh)
        s = a.segment(vh)
        wt = a.segment(VH_PREFIX + WT_CDRH3 + VH_SUFFIX)
        cdr = s.cdr3 if not s.error else None
        ref = wt.cdr3
        # Only aligned equal-length comparisons are claimed here.
        identity = sum(x == y for x,y in zip(cdr,ref))/len(ref) if cdr and len(cdr)==len(ref) else None
        result['cdrh3_by_scheme'][scheme] = {'sequence':cdr,'wt':ref,'identity':identity}
        if scheme == 'imgt':
            light=a.number(vl)
            result['numbering']={'vh_confidence':n.confidence,'vl_confidence':light.confidence,
                'vh_chain':n.chain,'vl_chain':light.chain}
            result['numbering_pass']=not n.error and not light.error and n.chain=='H' and light.chain in ['K','L'] and min(n.confidence,light.confidence)>=0.1
            result['canonical_cysteines_only']=(vh.count('C')==2 and vl.count('C')==2 and all(z.numbering.get('23')=='C' and z.numbering.get('104')=='C' for z in [n,light]))
    result['all_schemes_identity_below_70']=all(x['identity'] is not None and x['identity']<0.7 for x in result['cdrh3_by_scheme'].values())
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates',required=True)
    p.add_argument('--facs',required=True)
    p.add_argument('--spr',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    known=sorted(set(pd.read_csv(a.facs,usecols=['cdrh3']).cdrh3)|set(pd.read_csv(a.spr,usecols=['cdrh3']).cdrh3))
    rows=[]
    for row in pd.read_csv(a.candidates).fillna('').to_dict('records'):
        seq=row['cdrh3'];vh=row['vh_sequence'];vl=row['vl_sequence']
        nearest=process.extractOne(seq,known,scorer=Levenshtein.distance)
        row.update(annotate(vh,vl))
        row['nearest_provided_cdrh3']=nearest[0];row['nearest_provided_edit_distance']=int(nearest[1])
        row['paired_length']=len(vh)+len(vl)
        row['scfv_length_with_15aa_linker']=len(vh)+len(vl)+15
        row['strict_local_gate']=bool(row['numbering_pass'] and row['canonical_cysteines_only'] and row['all_schemes_identity_below_70'] and nearest[1]>=3 and row['paired_length']<=250)
        row['official_validation']='pending'
        rows.append(row)
    with open(a.output,'w') as f: json.dump(rows,f,indent=2,allow_nan=False)

if __name__=='__main__': main()
