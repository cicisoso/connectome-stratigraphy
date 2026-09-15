"""Download the public MaleCNS v1.0 files used by this repository into $FLYBRAIN_DATA (default ../data).

Files: the body annotation table (14 MB), the body statistics table (778 MB, needed only by
birth_consequences.py and dev_fig4.py) and the SWC skeletons of the 17,225 neurons listed in
results/swc_needed.txt (about 2 GB). Existing files are skipped, so the script can be re-run.
Usage: python download_data.py [--no-stats] [--workers 16]
"""
import os, sys, argparse, requests
from concurrent.futures import ThreadPoolExecutor, as_completed
BASE='https://storage.googleapis.com/flyem-male-cns/v1.0'
FLAT=f'{BASE}/connectome-data/flat-connectome'
SWC=f'{BASE}/segmentation/skeletons-malecns/skeletons-swc'
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.environ.get('FLYBRAIN_DATA', os.path.join(HERE,'..','data'))
def fetch(url,dst):
    if os.path.exists(dst) and os.path.getsize(dst)>0: return 'skip'
    r=requests.get(url,stream=True,timeout=120); r.raise_for_status()
    tmp=dst+'.part'
    with open(tmp,'wb') as f:
        for chunk in r.iter_content(1<<20): f.write(chunk)
    os.replace(tmp,dst); return 'ok'
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--no-stats',action='store_true',help='skip the 778 MB body-stats table'); ap.add_argument('--workers',type=int,default=16); a=ap.parse_args()
    os.makedirs(os.path.join(D,'swc'),exist_ok=True)
    tables=['body-annotations-male-cns-v1.0-minconf-0.5.feather']+([] if a.no_stats else ['body-stats-male-cns-v1.0-minconf-0.5.feather'])
    for t in tables: print(t,fetch(f'{FLAT}/{t}',os.path.join(D,t)))
    ids=[l.strip() for l in open(os.path.join(HERE,'..','results','swc_needed.txt')) if l.strip()]
    done=0; fail=[]
    with ThreadPoolExecutor(a.workers) as ex:
        futs={ex.submit(fetch,f'{SWC}/{i}.swc',os.path.join(D,'swc',f'{i}.swc')):i for i in ids}
        for f in as_completed(futs):
            try: f.result(); done+=1
            except Exception as e: fail.append((futs[f],str(e)))
            if done%1000==0: print(f'{done}/{len(ids)} skeletons')
    print(f'skeletons: {done} ok, {len(fail)} failed'); 
    if fail: print(fail[:10])
