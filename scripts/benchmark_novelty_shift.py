"""Diagnostic: does FACS discrimination survive distance from training data?"""
import argparse,json
from pathlib import Path
import pandas as pd
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
from sklearn.metrics import roc_auc_score,average_precision_score
p=argparse.ArgumentParser()
p.add_argument('--train',required=True)
p.add_argument('--predictions',required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--sample-size',type=int,default=2000)
a=p.parse_args()
train=pd.read_csv(a.train).seq.tolist()
test=pd.read_csv(a.predictions).query("split == 'test'")
x=test.sample(n=min(a.sample_size,len(test)),random_state=20261114).copy()
d=[]
for s in x.cdrh3:
 match=process.extractOne(s,train,scorer=Levenshtein.distance,score_cutoff=2)
 d.append(int(match[1]) if match else 3)
x['distance_to_train_capped_at_3']=d
report={'sample_seed':20261114,'sample_size':len(x),'interpretation':'Diagnostic only. Distance 3 means at least 3. No model selection on these results. Train-only distance, not full submission novelty corpus.','cohorts':{}}
for threshold in [0,2,3]:
 y=x[x.distance_to_train_capped_at_3>=threshold]
 report['cohorts'][f'train_distance_at_least_{threshold}']={'n':len(y),'high_prevalence':float(y.label.mean()) if len(y) else None,'roc_auc':float(roc_auc_score(y.label,y.p_facs_high)) if y.label.nunique()==2 else None,'average_precision':float(average_precision_score(y.label,y.p_facs_high)) if y.label.nunique()==2 else None}
a.output.write_text(json.dumps(report,indent=2)+'\n')
x.to_csv(a.output.with_suffix('.csv'),index=False)
