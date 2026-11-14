"""Use the same target/format as WT for a bounded shortlist refolding batch."""
import argparse,copy,csv,json
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--controls',type=Path,required=True)
p.add_argument('--candidates',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--count',type=int,default=2)
a=p.parse_args()
base=next(x for x in json.loads(a.controls.read_text())['queries'] if x['query_id']=='control_wt')
with a.candidates.open() as f: rows=list(csv.DictReader(f))[:a.count]
queries=[]
for row in rows:
 q=copy.deepcopy(base);q['query_id']='stat_'+row['cdrh3']
 for chain in q['chains']:
  if chain['chain_id']=='A':chain['sequence']=row['vl_sequence']
  if chain['chain_id']=='B':chain['sequence']=row['vh_sequence']
 queries.append(q)
a.output.write_text(json.dumps({'queries':queries},indent=2)+'\n')
