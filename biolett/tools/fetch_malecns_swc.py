"""Fetch MaleCNS v1.0 SWC skeletons listed in a text file (one bodyId per line) into $FLYBRAIN_DATA/swc/ (default data/swc/)."""
import os, sys, requests
from concurrent.futures import ThreadPoolExecutor, as_completed
SWC='https://storage.googleapis.com/flyem-male-cns/v1.0/segmentation/skeletons-malecns/skeletons-swc'
D=os.path.join(os.environ.get('FLYBRAIN_DATA',os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'data')),'swc')
os.makedirs(D,exist_ok=True)
ids=[l.strip() for l in open(sys.argv[1]) if l.strip()]
def fetch(i):
    dst=f'{D}/{i}.swc'
    if os.path.exists(dst) and os.path.getsize(dst)>0: return 'skip'
    r=requests.get(f'{SWC}/{i}.swc',timeout=120)
    if r.status_code!=200: return f'http{r.status_code}'
    open(dst+'.part','wb').write(r.content); os.replace(dst+'.part',dst); return 'ok'
res={}
with ThreadPoolExecutor(24) as ex:
    futs={ex.submit(fetch,i):i for i in ids}
    for k,f in enumerate(as_completed(futs)):
        try: s=f.result()
        except Exception as e: s='err'
        res[s]=res.get(s,0)+1
        if (k+1)%1000==0: print(k+1,res,flush=True)
print('final',res)
