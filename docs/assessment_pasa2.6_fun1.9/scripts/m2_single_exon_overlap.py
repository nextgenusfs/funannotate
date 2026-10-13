import sys,re,collections,bisect
asm=collections.defaultdict(list)
for l in open(sys.argv[1]):
    if l.startswith('#'): continue
    c=l.rstrip('\n').split('\t')
    if len(c)<9: continue
    m=re.search(r'Target=(\S+)',c[8]); 
    asm[m.group(1)].append((c[0],int(c[3]),int(c[4]),c[6]))
single=[];spliced=collections.defaultdict(list)
for a,ex in asm.items():
    ch=ex[0][0]; s=min(e[1] for e in ex); e=max(e[2] for e in ex)
    if len(ex)==1: single.append((ch,s,e,ex[0][3]))
    else: spliced[ch].append((s,e,ex[0][3]))
for ch in spliced: spliced[ch].sort()
ov=cont=0
for ch,s,e,st in single:
    L=spliced.get(ch,[]); hit=False; inside=False
    for (ss,ee,st2) in L:
        if ss>e: break
        if ee>=s: hit=True; inside = inside or (ss<=s and ee>=e)
    ov+=hit; cont+=inside
n=len(asm); print(f"assemblies={n} single={len(single)} ({100*len(single)/n:.1f}%) single_overlapping_spliced={ov} ({100*ov/max(1,len(single)):.1f}%) single_contained_in_spliced_span={cont} ({100*cont/max(1,len(single)):.1f}%)")
