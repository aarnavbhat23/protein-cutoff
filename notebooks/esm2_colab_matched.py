# Matched test: is a protein LM's familiarity about the family or the sequence?
# All proteins first seen in UniParc on/after 2023-01-01 (none can be in ESM-2's UniRef 2021_04 data),
# drawn from the same 2024-26 UniProt additions. new_family: whole UniRef50 cluster new since 2023;
# old_family_matched: cluster has pre-2023 members, matched 2:1 on length and reviewed status.
import json, math, subprocess, sys, time, urllib.parse, urllib.request
subprocess.run([sys.executable, "-m", "pip", "-q", "install", "fair-esm"], check=True)
import esm, numpy as np, torch
NEW = "A0AA51Z0J6 C0HME6 C0HMF8 P0DYA9 A0AAW1I778 A0AAW1LG41 P0DRM7 P0DXE4 P0DXE6 P0DXL5 P0DY42 A0AA38C2V5 A0AA38GMP9 A0AAD6NFM2 A0AAD9S2F9 A0AAU7VAB8 C0HMD0 C0HMD6 C0HME5 C0HMG8 P0DRB9 P0DRC0 P0DRC3 P0DRD4 P0DRD5 P0DRH7 P0DRI6 P0DRI7 P0DXH9 P0DXK7 P0DXK8 P0DXK9 P0DXL1 P0DXL3 P0DXL6 P0DXN9 P0DXW7 P0DXX1 P0DYE9 P9WEG8 P9WEK7 P9WEL0 P9WEL5 A0A1D0BN92 C0HMC8 C0HML8 E0D209 P0DRC6 P0DRI0 P0DRI1 P0DRI2 P0DRM2 P0DXA0 P0DXA1 P0DXO0 P0DXW3 P0DXW4 P0DXW5 P0DXW8 P0DXX6 P0DY00 P0DY21 P0DY38 P0DY40 P0DY41 P0DY43 P0DY44 P9WEI1 P9WEK4 P9WEL4 A0A410N6V5 A0AAN2H6H5 A0AAN2H888 A0AAN2L0I3 C0HMD7 P0DRI5 P0DYD0 P9WEI0 P9WEJ2 P9WEJ5 A0AA35RBR7 A0AA35X3E9 A0AA40LVM8 A0AA35X2A1 A0AA40LFM2 A0AA35CN54 A0AA35QXR6 A0AA35W5V5 A0AA36BJS0 A0AA36NFT9 A0AA40C8I9 A0AA40FEV9 A0AA40HUQ8 A0AA40HV03 A0AA88L9Q4 A0AA89BZV0 A0AA88Y5F9 A0AA89BMD2 A0AA86QHB7 A0AA86TG18 A0AA88IRP0 A0AA88YC86 A0AA89C6I6 A0AA89C9N1 A0AA88I4V7 A0AAA9S8H4 A0AAA9TLL4 A0AAA9SUZ2 A0AAA9YHZ6 A0AAF7A4B6 A0AAF6Z103 A0AAF6ZA81 A0AAF6YLR3 A0AAJ7THV5 A0AAQ5BHK0 A0AAQ5BGQ7 A0AAQ5BHJ4 A0AAQ5BI73 A0AAQ5BHX1 A0AAQ5BHJ5 A0AAN7VGE2 A0AAQ4EI19 A0AAQ4P8F5 A0AAN9TJF4 A0AAN9TLU8 A0AAN9Z068 A0AAN9ZIB5 A0AAN7XHV4 A0AAN7XQD9 A0AAN8JV67 A0AAW0N8T1 A0AAV2JLF8 A0AAV4NI31 A0AAV7JAP6 A0AAV7T028 A0AAW0P6Z4 A0AB32TAJ9 A0AB32TTZ9 A0AB39TZ30 A0AB34GI97 A0ABB0MVI8 A0ABB0MV20 A0ABD0LEB6 A0ABD0LYH4 A0ABD2QNE0 A0ABD2WFQ1 A0ABD3VJ77 A0ABK0LP96 A0ABK0LN78 A0AC58HZU9 A0AC58GR43 A0AC58IPS6 A0AC58JG01 A0AC58I5K5 A0ACI8S133 A0ACI8P6N6 A0ACI8QP07".split()
OLD = "A0AAD6NFV5 P0DYC0 A0A9E8S4N7 G2WYE8 A0A9W9TFI6 P0DRM5 P0DRL9 P0DRD6 A0AAV8T0T3 A0A9E7LUL3 P0DXL0 P0DY11 P0DY15 E2ITZ3 A0A9E7S4M3 P0DXJ6 C0HME7 C0HMF1 A0AAD6J2M3 A0A1L8ENT6 P0DO64 C0HMB2 P9WEI6 P9WEL6 A0AAD6J5F2 P0DXH7 P0DXE5 P9WEE7 P0DRN5 P0DXV9 G1SV32 A0A9Y1LQX3 P0DRM4 P0DO63 C0HMF4 P0DRH9 C0HMD8 P0DY18 D2XUU9 P0DRD0 P0DY26 C0HM97 C0HM71 C0HMC5 P0DRI9 P0DRC9 P0DRL5 P0DXT0 P0DX56 C0HMC6 P0DRH5 C0HMB6 P0DRC4 P0DRM3 C0HM61 P0DY39 P0DXY3 P9WEK6 A0AAD6NH75 P0DXJ3 P0DXW2 P0DY10 A0AA51Z3S7 Q47595 P0DRJ4 G0SF39 P9WEE6 P9WEH7 A0AA51VJB7 A0A9E8M4C5 G1TTQ5 C0HMD2 P0DRM1 P0DYB2 Q50KZ9 C0HMC9 A0A9E7S518 A0A9E9KUX0 P9WEF2 A0A9E7S5E3 A0AAW1M2U7 P0DXJ0 A0A9Y1YY28 A0AAD6NMH5 A0A9W9TGP7 A0A9E7LUR3 P0DY04 P0DRH3 P0DXJ4 C0HMB0 P0DRB0 P0DYB1 P0DRM6 A0AAD9S265 P0DX64 P0DYB0 P0DRL6 P0DRM0 F2ST24 P0DXO4 P0DYC1 P0DXV1 P0DXH3 P9WEK5 P9WEM1 G1SWM1 P9WEK8 A0A1V0QSH7 A0A9E9KU69 P0DXW0 P0DO76 A0A9E9KU71 A0AA51VIL5 A0A2Z5P0W7 A0AA51VIL7 A0A9W3N2X3 A0A9E8G339 P0DXK6 P9WEH5 P0DYC2 A0ABQ8WBK9 J9XYC6 A0AAD6ISS3 A0AAD6ITA8 A0AAD6NLP7 P0DY83 P0DY12 A0A9W3N256 P0DXH2 A0AAW1MQS1 P0DXK4 A0A2I1CSH6 A0AAD6IYT8 A0AAW1IQ05 A0AAD6J5T3 P9WEL2 C0HM44 P0DXH6 P9WEL7 P0DO65 P0DYB5 P9WEK9 A0AAD6NF83 P0DXH8 A0AAW1NEA3 A0AAV8SNK8 P9WEH8 A0AAW1J9R8 A0A5B8ND22 A0AAW1J1B8 A0AAD9S169 P9WEI5 A0A291FB00 A0AAW1JA93 P0DXF1 P0DO78 Q7S9V5 P9WEE8 P9WEK2 P0DY68 A0AC58HZV4 A0ACI8V743 A0AAA9YHL0 A0AAQ5BH68 A0AAA9U040 A0ABB0MV68 A0AB39TZ26 A0ABK0M239 A0AA88PJF5 A0ABK0LH83 A0AA35NZ37 A0AAA9SWA1 A0AA39M7X6 A0AA86V8N0 A0AA88Q631 A0ABD2DP67 A0AAU8J6L0 A0ACI8UUH4 A0AA97JDY6 A0AAF6ZCD0 A0AA38I268 A0AA39VLF3 A0AAQ5BII0 A0ACI8SPQ4 A0AAF6YM64 A0ABB0MVJ0 A0ABC8K0S2 A0ABK0LKE5 A0AA90H4N6 A0AAF6YNF5 A0AAA9SIR3 A0AB32TEC1 A0AAA9TYC3 A0AAF6YZK2 A0AA41MVY2 A0AAA9TFI4 A0AA88PZ64 A0AA88SFJ8 A0AA36DPU6 A0AA38UHY2 A0AA35KD37 A0AA40FSY4 A0AA88T406 A0AA88M8M1 A0ABK0LH29 A0AAF6YP72 A0AAA9SXA6 A0AAQ5BH33 A0AAF6ZGR8 A0AAW0XNH1 A0AAA9TD36 A0AB72X5E8 A0AB72XBE0 A0ABC9XAN0 A0AAA9TPX3 A0AAA9TKF7 A0AA96LIY1 A0AAA9S0M9 A0AAF6ZG78 A0AB32TEV5 A0AA47MJF2 A0AA40I0V5 A0ABK0M5W4 A0ABK0LLA4 A0AB32TLG9 A0ABK0L5Q0 A0AA97KZB7 A0AB74U7F0 A0AAZ3S3G4 A0ACI8UYA1 A0AAR2LYL8 A0AA35L3A6 A0AAQ5BH99 A0ABB0MV42 A0AA88MJ09 A0AAA9YHK8 A0ABK0L4V6 A0AC58HXB1 A0AAQ5BH58 A0AB12YBX2 A0AA35LJ55 A0AA36FBW7 A0ACI8V6A4 A0AA39HMR6 A0AAA9TNA2 A0AA47NWU0 A0AA34QVI4 A0AA39CAW8 A0AAN8CLN4 A0AAQ5XQS3 A0AAA9SL17 A0AAQ5BI04 A0AA34QVV2 A0AAA9SUD1 A0ACI8TDQ7 A0AA96LMB0 A0AA35NZI6 A0AA96NL62 A0AA37GCQ2 A0AAF6Z9D1 A0AAW2ATI2 A0AB32TYW9 A0AA88PFI6 A0AAJ7WY93 A0AAA9TT54 A0AAA9TBT5 A0AAA9T6V1 A0AAN8KPI9 A0AA35NZL8 A0AA47N9S7 A0AA88SY32 A0AAA9SG76 A0AAQ5BI51 A0AB72X7E3 A0AAT9PAB7 A0AB74U6W3 A0AAF6Z2F7 A0AAT9J3Z5 A0AB72WYH6 A0AB32T459 A0AA34QW04 A0ABB0MV96 A0ABB0MV00 A0AB72WZ63 A0AA39H385 A0AA97KCZ2 A0AA41MX00 A0AA40HPF2 A0AA88LWY0 A0ABI0P1A6 A0AA88SXR8 A0AA86RZM8 A0AA37T108 A0AAA9S698 F6QMY2 A0AAF6ZJ25 A0AA88MBT9 A0ABB0MVH8 A0AA38TLP7 A0AA88N844 A0AA97JL19 A0AAA9TJ11 A0AA40I5Q4 A0AA97K388 A0AA50HPG1 A0AA39FDC4 A0AAA9SPF0 A0AAF6YGY8 A0AAJ2JZI4 A0AAQ5BHS4 A0ACI8QKD0 A0AAA9U024 A0ABK0LY56 A0ABK0KZC3".split()
MODELS = ["esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D", "esm2_t30_150M_UR50D"]
SAMPLE = 64


def get(url):
    for i in range(7):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "protein-cutoff-audit"}), timeout=180) as r:
                return r.read()
        except Exception:
            if i == 6:
                raise
            time.sleep(2 ** i)


seqs = {}
accs = NEW + OLD
for i in range(0, len(accs), 50):
    q = " OR ".join(f"accession:{a}" for a in accs[i:i + 50])
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": q, "format": "tsv", "fields": "accession,sequence", "size": 500})
    for line in get(url).decode().strip().split("\n")[1:]:
        a, s = line.split("\t")
        seqs[a] = s
rows = [(a, seqs[a], "new_family") for a in NEW if a in seqs] + [(a, seqs[a], "old_family") for a in OLD if a in seqs]
print("proteins:", {g: sum(r[2] == g for r in rows) for g in ("new_family", "old_family")}, flush=True)
dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")


@torch.no_grad()
def ppl(model, alphabet, seq, seed=0):
    seq = seq[:1022]
    _, _, toks = alphabet.get_batch_converter()([("p", seq)])
    toks = toks.to(dev)
    pos_all = np.sort(np.random.default_rng(seed).choice(len(seq), size=min(len(seq), SAMPLE), replace=False))
    pos = torch.from_numpy(pos_all).to(dev)
    x = toks.repeat(len(pos), 1)
    x[torch.arange(len(pos), device=dev), pos + 1] = alphabet.mask_idx
    with torch.autocast("cuda", dtype=torch.float16, enabled=dev.type == "cuda"):
        logits = model(x)["logits"][torch.arange(len(pos), device=dev), pos + 1]
    nll = -torch.log_softmax(logits.float(), -1).gather(1, toks[0, pos + 1][:, None])[:, 0].sum().item()
    return math.exp(nll / len(pos_all)), len(seq)


ALL = {"design": "first seen >= 2023-01-01; new vs old UniRef50 family; matched on length and reviewed status", "models": {}}
for name in MODELS:
    model, alphabet = getattr(esm.pretrained, name)()
    model = model.eval().to(dev)
    res = [(g, *ppl(model, alphabet, s)) for a, s, g in rows]
    rng = np.random.default_rng(0)
    v = {g: np.array([p for gg, p, L in res if gg == g]) for g in ("new_family", "old_family")}
    y = np.log([p for g, p, L in res])
    X = np.column_stack([np.ones(len(res)), np.log([L for g, p, L in res]), [g == "old_family" for g, p, L in res]]).astype(float)
    beta = np.linalg.lstsq(X, y, rcond=None)[0][2]
    boots = []
    for _ in range(2000):
        i = rng.integers(0, len(y), len(y))
        boots.append(np.linalg.lstsq(X[i], y[i], rcond=None)[0][2])
    d = [np.median(rng.choice(v["old_family"], len(v["old_family"]))) - np.median(rng.choice(v["new_family"], len(v["new_family"]))) for _ in range(4000)]
    ALL["models"][name] = {g: {"n": int(len(x)), "median": round(float(np.median(x)), 3), "q25": round(float(np.percentile(x, 25)), 3),
                               "q75": round(float(np.percentile(x, 75)), 3)} for g, x in v.items()}
    ALL["models"][name]["median_old_minus_new"] = [round(float(np.median(v["old_family"]) - np.median(v["new_family"])), 3),
                                                    round(float(np.percentile(d, 2.5)), 3), round(float(np.percentile(d, 97.5)), 3)]
    ALL["models"][name]["length_adjusted_ratio_old_vs_new"] = [round(float(np.exp(beta)), 3), round(float(np.exp(np.percentile(boots, 2.5))), 3),
                                                                round(float(np.exp(np.percentile(boots, 97.5))), 3)]
    print(name, json.dumps(ALL["models"][name]), flush=True)
    del model
    torch.cuda.empty_cache()
print("===RESULT_JSON===")
print(json.dumps(ALL))
print("===END===")
