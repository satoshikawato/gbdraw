from pathlib import Path
import sys,json,gzip,hashlib
root=Path.cwd();sys.path[:0]=[str(root),str(root/'tools')]
import benchmark_protein_comparison as b
from protein_comparison_browser import run_path_browser
pc,cc=b.load_source(root)
inputs=b.path_browser_inputs(root,pc,cc,names=('gallery-collinear',))
result={'purpose':'One separate memory/counter diagnostic per viewport, no timing comparison','source':b.source_info(root),'result':run_path_browser(root,inputs,1,'memory')}
result['result']['artifact']='S08 source SPA and final generated local wheel'
p=root/'docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data/s08-browser-memory.json.gz'
p.write_bytes(gzip.compress((json.dumps(result,indent=2)+'\n').encode(),mtime=0))
