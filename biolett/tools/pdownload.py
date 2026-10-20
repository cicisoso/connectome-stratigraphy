"""Parallel HTTP range downloader: python pdownload.py URL DEST TOTAL_BYTES [workers]. Resumes from DEST.part.* chunk files."""
import os, sys, requests, threading, time
from concurrent.futures import ThreadPoolExecutor
url, dest, total = sys.argv[1], sys.argv[2], int(sys.argv[3]); W = int(sys.argv[4]) if len(sys.argv) > 4 else 16
prefix = dest + '.prefix'                       # bytes already downloaded by a previous single-stream run
start0 = os.path.getsize(prefix) if os.path.exists(prefix) else 0
CH = 64 << 20                                   # 64 MB chunks
chunks = [(i, s, min(s + CH, total) - 1) for i, s in enumerate(range(start0, total, CH))]
def get(c):
    i, s, e = c; part = f'{dest}.part{i:05d}'
    need = e - s + 1
    for attempt in range(30):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have >= need: return i
        try:
            r = requests.get(url, headers={'Range': f'bytes={s+have}-{e}'}, stream=True, timeout=120)
            if r.status_code not in (200, 206): raise RuntimeError(r.status_code)
            with open(part, 'ab') as f:
                for b in r.iter_content(1 << 20): f.write(b)
        except Exception as ex:
            time.sleep(5 + attempt)
    raise RuntimeError(f'chunk {i} failed')
t0 = time.time(); done = 0
with ThreadPoolExecutor(W) as ex:
    for i in ex.map(get, chunks):
        done += 1
        if done % 5 == 0: print(f'{done}/{len(chunks)} chunks, {time.time()-t0:.0f}s', flush=True)
with open(dest, 'wb') as out:
    if start0:
        with open(prefix, 'rb') as f:
            while (b := f.read(1 << 24)): out.write(b)
    for i, s, e in chunks:
        p = f'{dest}.part{i:05d}'
        with open(p, 'rb') as f:
            while (b := f.read(1 << 24)): out.write(b)
assert os.path.getsize(dest) == total, (os.path.getsize(dest), total)
for i, s, e in chunks: os.remove(f'{dest}.part{i:05d}')
if start0: os.remove(prefix)
print('done', dest, os.path.getsize(dest), f'{time.time()-t0:.0f}s')
