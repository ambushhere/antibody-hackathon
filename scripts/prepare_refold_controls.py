"""Add WT and a prespecified FACS-low control to the optimisation pilot queries."""
import argparse, copy, csv, json
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--prepared',type=Path,required=True)
p.add_argument('--facs',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
queries=json.loads(a.prepared.read_text())['queries']
selected=queries[:2]
wt=copy.deepcopy(queries[0]);wt['query_id']='control_wt'
heavy=next(c for c in wt['chains'] if c['chain_id']=='B')
# Library window is fixed between CSR and YWGQGTLVTVSS, independent of CDR scheme.
prefix,suffix=heavy['sequence'].split('CSR',1)
assert suffix.endswith('YWGQGTLVTVSS')
heavy['sequence']=prefix+'CSR'+'WGGDGFYAMD'+'YWGQGTLVTVSS'
with a.facs.open() as f:
 low=next(r for r in csv.DictReader(f) if r['binding_class']=='low' and 'C' not in r['cdrh3'])
negative=copy.deepcopy(wt);negative['query_id']='control_facs_low'
next(c for c in negative['chains'] if c['chain_id']=='B')['sequence']=low['vh_sequence']
a.output.write_text(json.dumps({'queries':[wt,negative]+selected},indent=2)+'\n')
