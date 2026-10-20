# ESM-2 pseudo-perplexity by sequence age -- self-contained, for a Colab GPU runtime.
# Inputs are fetched from public sources only: LiveProteinBench (pinned commit) and UniProt.
# Groups: benchmark proteins first seen before / after ESM-2's training data (UniRef 2021_04),
# and NovelFam (sequence and 50% cluster new after 2023-01-01, no Swiss-Prot relative).
import json, math, subprocess, sys, time, urllib.parse, urllib.request

subprocess.run([sys.executable, "-m", "pip", "-q", "install", "fair-esm"], check=True)
import esm
import numpy as np
import torch

COMMIT = "8b9710932511f816651fc1a7e8a1bc85c3c77767"
TASKS = ["EC_number", "GO_biological_process", "GO_cellular_component", "GO_molecular_function", "active_site",
         "catalytic_activity", "cofactor", "motif_position", "pathway", "ph", "temperature", "transmembrane"]
NOVELFAM = ("P0DXT9 P0DY40 P0DXX0 P0DRC3 P0DY43 C0HMF8 P0DXX2 P0DRC0 P0DRI7 A0A410N6V5 P0DRM7 C0HML8 C0HMG8 P0DXN9 "
            "P9WEJ5 P0DY41 A0AAN2L0I3 P0DRI2 P0DXO1 P0DRM2 C0HME6 P0DXE4 P9WEI0 P0DRH7 P9WEL0 C0HMC8 A0A1D0BN92 P9WEL4 "
            "A0AAU7VAB8 A0AAN2H888 P9WEK7 P0DXW4 P0DXW7 P0DY21 C0HME5 C0HMD7 P0DXW8 P0DXA1 A0AAN2H6H5 P0DXX6 P9WEJ2 "
            "P0DXW3 C0HMD6 P0DRI5 A0ACI8P6N6 A0AAA9SUZ2 A0AA86QTQ0 A0AAF6ZA81 A0AA36NFT9 A0ACI8S133 A0AB32TTZ9 A0AC58IPS6 "
            "A0ABD2QNE0 A0AAQ5BHX1 A0AA89BZV0 A0AAN9ZIB5 A0AAN9T4W0 A0AA86TG18 A0AAN8JV67 A0ABD0LYH4 A0AAA9YHZ6 A0AA88I4V7 "
            "A0AAF6Z103 A0AA40LFM2").split()
ESM2_DATA = "2021-04-01"
MODELS = ["esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D", "esm2_t30_150M_UR50D"]
SAMPLE = 64  # masked positions per protein, seeded


def get(url, data=None):
    for i in range(7):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers={"User-Agent": "protein-cutoff-audit"}), timeout=180) as r:
                return r.read()
        except Exception:
            if i == 6:
                raise
            time.sleep(2 ** i)


def tsv(text):
    lines = text.strip("\n").split("\n")
    keys = lines[0].split("\t")
    return [dict(zip(keys, l.split("\t"))) for l in lines[1:]]


def idmap_first_seen(ids):
    out = {}
    for i in range(0, len(ids), 5000):
        data = urllib.parse.urlencode({"ids": ",".join(ids[i:i + 5000]), "from": "UniProtKB_AC-ID", "to": "UniParc"}).encode()
        job = json.loads(get("https://rest.uniprot.org/idmapping/run", data))["jobId"]
        while True:
            st = json.loads(get(f"https://rest.uniprot.org/idmapping/status/{job}"))
            if st.get("jobStatus") in (None, "FINISHED") or "results" in st:
                break
            time.sleep(3)
        for r in tsv(get(f"https://rest.uniprot.org/idmapping/uniparc/results/stream/{job}?format=tsv&fields=upi,first_seen").decode()):
            out[r["From"]] = r["First seen"]
    return out


seqs = {}
for t in TASKS:
    for q in json.loads(get(f"https://raw.githubusercontent.com/Rongdingyi/LiveProteinBench/{COMMIT}/dataset/QA/{t}.json")):
        seqs.setdefault(q["Protein Image"][0].split("/")[1], q["Protein Sequence"])
nf = {}
for i in range(0, len(NOVELFAM), 50):
    q = " OR ".join(f"accession:{a}" for a in NOVELFAM[i:i + 50])
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": q, "format": "tsv", "fields": "accession,sequence", "size": 500})
    for r in tsv(get(url).decode()):
        nf[r["Entry"]] = r["Sequence"]
fs = idmap_first_seen(sorted(seqs) + sorted(nf))
rows = [(a, s, fs.get(a, ""), "before_esm2" if fs.get(a, "9") < ESM2_DATA else "after_esm2") for a, s in seqs.items()]
rows += [(a, s, fs.get(a, ""), "novelfam") for a, s in nf.items()]
print(f"proteins: benchmark {len(seqs)}, novelfam {len(nf)}; groups:",
      {g: sum(r[3] == g for r in rows) for g in ("before_esm2", "after_esm2", "novelfam")}, flush=True)

dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("device:", dev, torch.cuda.get_device_name(0) if dev.type == "cuda" else "", flush=True)


@torch.no_grad()
def ppl(model, alphabet, seq, sample=SAMPLE, batch=64, seed=0):
    seq = seq[:1022]
    _, _, toks = alphabet.get_batch_converter()([("p", seq)])
    toks = toks.to(dev)
    L = len(seq)
    pos_all = np.sort(np.random.default_rng(seed).choice(L, size=min(L, sample), replace=False))
    nll = 0.0
    for s in range(0, len(pos_all), batch):
        pos = torch.from_numpy(pos_all[s:s + batch]).to(dev)
        x = toks.repeat(len(pos), 1)
        x[torch.arange(len(pos), device=dev), pos + 1] = alphabet.mask_idx
        with torch.autocast("cuda", dtype=torch.float16, enabled=dev.type == "cuda"):
            logits = model(x)["logits"][torch.arange(len(pos), device=dev), pos + 1]
        nll -= torch.log_softmax(logits.float(), -1).gather(1, toks[0, pos + 1][:, None])[:, 0].sum().item()
    return math.exp(nll / len(pos_all)), L


def summarise(res):
    rng = np.random.default_rng(0)
    g = {k: np.array([r["ppl"] for r in res if r["group"] == k]) for k in ("before_esm2", "after_esm2", "novelfam")}
    out = {k: {"n": int(len(v)), "median": round(float(np.median(v)), 3),
               "q25": round(float(np.percentile(v, 25)), 3), "q75": round(float(np.percentile(v, 75)), 3)} for k, v in g.items()}
    for a, b in (("before_esm2", "novelfam"), ("after_esm2", "novelfam"), ("before_esm2", "after_esm2")):
        d = [np.median(rng.choice(g[a], len(g[a]))) - np.median(rng.choice(g[b], len(g[b]))) for _ in range(4000)]
        out[f"median_{a}_minus_{b}"] = [round(float(np.median(g[a]) - np.median(g[b])), 3),
                                        round(float(np.percentile(d, 2.5)), 3), round(float(np.percentile(d, 97.5)), 3)]
    # length-adjusted: log ppl ~ 1 + log L + [before] + [after]   (NovelFam is the reference)
    y = np.log([r["ppl"] for r in res])
    X = np.column_stack([np.ones(len(res)), np.log([r["L"] for r in res]),
                         [r["group"] == "before_esm2" for r in res], [r["group"] == "after_esm2" for r in res]]).astype(float)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    boots = []
    for _ in range(2000):
        i = rng.integers(0, len(y), len(y))
        boots.append(np.linalg.lstsq(X[i], y[i], rcond=None)[0])
    boots = np.array(boots)
    for j, name in ((2, "before_esm2"), (3, "after_esm2")):
        out[f"length_adjusted_ratio_{name}_vs_novelfam"] = [round(float(np.exp(beta[j])), 3),
                                                            round(float(np.exp(np.percentile(boots[:, j], 2.5))), 3),
                                                            round(float(np.exp(np.percentile(boots[:, j], 97.5))), 3)]
    return out


ALL = {"benchmark_commit": COMMIT, "sample_positions": SAMPLE, "esm2_training_data": "UniRef 2021_04", "models": {}}
PER = []
for name in MODELS:
    t0 = time.time()
    model, alphabet = getattr(esm.pretrained, name)()
    model = model.eval().to(dev)
    res = []
    for i, (a, s, f, g) in enumerate(rows):
        p, L = ppl(model, alphabet, s)
        res.append({"acc": a, "group": g, "first_seen": f, "L": L, "ppl": p})
        if i % 200 == 0:
            print(f"{name}: {i}/{len(rows)}  {time.time() - t0:.0f}s", flush=True)
    ALL["models"][name] = summarise(res)
    PER += [{**r, "model": name} for r in res]
    print(name, json.dumps(ALL["models"][name]), flush=True)
    del model
    torch.cuda.empty_cache()

with open("esm2_per_protein.tsv", "w") as fh:
    fh.write("model\taccession\tgroup\tfirst_seen\tlength\tpseudo_perplexity\n")
    for r in PER:
        fh.write(f"{r['model']}\t{r['acc']}\t{r['group']}\t{r['first_seen']}\t{r['L']}\t{r['ppl']:.4f}\n")
print("===RESULT_JSON===")
print(json.dumps(ALL))
print("===END===")
