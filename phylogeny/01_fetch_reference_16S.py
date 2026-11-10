"""Build the species list and fetch one reference 16S sequence for each.

WHAT IT DOES
    Takes the most abundant species detected in this study, retrieves a
    full-length 16S rRNA record for each from NCBI, and records exactly which
    accession every tip came from.

INPUT
    data/counts.csv

OUTPUT
    output/phylo_species_list.txt       the species that entered the tree
    output/phylo_sequences_raw.fasta    one sequence per tip, named by species
    output/phylo_tip_provenance.csv     query_species, accession, fetched_as,
                                        retrieved_title, seq_length
    output/phylo_failed_species.txt     species with no usable record

    phylo_tip_provenance.csv is not optional bookkeeping.  The NCBI searches
    are ranked by relevance, so the top hit for a species can change between
    runs; without the accessions the tree cannot be rebuilt by anyone,
    including you.  Deposit it with the manuscript.

SPECIES SELECTION
    The TOPN most abundant species by total read count across all 24
    libraries, after removing names that cannot be placed on a species tree:
    "sp.", "bacterium", "endosymbiont", "Candidatus", "uncultured", "str.",
    "cf.", "aff.".  The filtered list is written out so the selection is
    auditable.

THE GENUS FALLBACK
    If no species-level record is found, a type-strain record for the GENUS is
    used instead and the tip keeps the species label.  That tip is therefore
    NOT that species - it is a congener standing in for it.  Every such tip is
    marked fetched_as = "genus_fallback" here and gets its own ring in the iTOL
    output, and species-level resolution must not be claimed for it.

BEFORE RUNNING
    Set EMAIL to your own address; NCBI requires it and will block anonymous
    bulk queries.  An API key raises the rate limit from 3 to 10 requests per
    second and is worth setting for a list this long:
    https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

    Runtime is roughly 20-30 minutes for 500 species without a key.
"""
import os, sys, time
import pandas as pd
from Bio import Entrez, SeqIO

EMAIL = os.environ.get("NCBI_EMAIL", "")
API_KEY = os.environ.get("NCBI_API_KEY", "")
TOPN = 500
MINLEN, MAXLEN = 1000, 1800
RETRIES = 3
SAVE_EVERY = 50

if not EMAIL:
    sys.exit("set NCBI_EMAIL in the environment before running "
             "(NCBI requires a contact address)")
Entrez.email = EMAIL
if API_KEY:
    Entrez.api_key = API_KEY
PAUSE = 0.12 if API_KEY else 0.40

DROP = [" sp.", "bacterium", "endosymbiont", "candidatus", "uncultured",
        " str.", " cf.", " aff.", "unknown", "unclassified"]

# ----------------------------------------------------------- species list
d = pd.read_csv("data/counts.csv").dropna(subset=["species"])
samp = sorted([c for c in d.columns if c[0] in "RS" and len(c) == 4])
assert len(samp) == 24, f"expected 24 sample columns, found {len(samp)}"

tot = d[samp].fillna(0).sum(axis=1)
tot.index = d.species.values
tot = tot.groupby(level=0).sum().sort_values(ascending=False)

keep = [s for s in tot.index if not any(x in s.lower() for x in DROP)]
species_list = keep[:TOPN]
open("output/phylo_species_list.txt", "w").write("\n".join(species_list))
print(f"{len(tot)} species detected -> {len(keep)} placeable -> "
      f"top {len(species_list)} kept")
print(f"  most abundant: {species_list[0]} ({int(tot[species_list[0]]):,} reads)")
print(f"  least abundant kept: {species_list[-1]} "
      f"({int(tot[species_list[-1]]):,} reads)")

# ----------------------------------------------------------- NCBI retrieval
def esearch_one(term):
    h = Entrez.esearch(db="nucleotide", term=term, retmax=1, sort="relevance")
    r = Entrez.read(h); h.close()
    return r["IdList"][0] if r["IdList"] else None

def efetch_one(uid):
    h = Entrez.efetch(db="nucleotide", id=uid, rettype="fasta", retmode="text")
    rec = SeqIO.read(h, "fasta"); h.close()
    return rec

def fetch_16S(species):
    """Return (record, mode) or (None, None). mode is 'species' or 'genus_fallback'."""
    slen = f"{MINLEN}:{MAXLEN}[SLEN]"
    queries = [
        (f'"{species}"[Organism] AND 16S[Title] AND rRNA[Title] AND {slen}',
         "species"),
        (f'"{species.split()[0]}"[Organism] AND 16S[Title] AND rRNA[Title] '
         f'AND {slen} AND "type strain"[Title]', "genus_fallback"),
    ]
    for attempt in range(RETRIES):
        try:
            for term, mode in queries:
                uid = esearch_one(term)
                time.sleep(PAUSE)
                if uid:
                    rec = efetch_one(uid)
                    time.sleep(PAUSE)
                    return rec, mode
            return None, None
        except Exception as e:
            print(f"    attempt {attempt + 1} failed: {e}")
            time.sleep(2 ** attempt)
    return None, None

seqs, prov, failed = [], [], []
print(f"\nfetching {len(species_list)} sequences "
      f"({'with' if API_KEY else 'without'} an API key)\n")

for i, species in enumerate(species_list, 1):
    rec, mode = fetch_16S(species)
    if rec is None:
        failed.append(species)
        print(f"[{i:3d}/{len(species_list)}] {species:55s} NOT FOUND")
        continue
    accession = rec.id.split()[0]
    title = rec.description
    rec.id = rec.name = species.replace(" ", "_").replace("/", "_")
    rec.description = ""
    seqs.append(rec)
    prov.append(dict(query_species=species, tip_label=rec.id,
                     accession=accession, fetched_as=mode,
                     retrieved_title=title, seq_length=len(rec.seq)))
    flag = "" if mode == "species" else "   <-- GENUS FALLBACK"
    print(f"[{i:3d}/{len(species_list)}] {species:55s} {accession}{flag}")
    if i % SAVE_EVERY == 0:
        SeqIO.write(seqs, "output/phylo_sequences_raw.fasta", "fasta")
        pd.DataFrame(prov).to_csv("output/phylo_tip_provenance.csv", index=False)
        print(f"  -- progress saved: {len(seqs)} sequences --")

SeqIO.write(seqs, "output/phylo_sequences_raw.fasta", "fasta")
P = pd.DataFrame(prov)
P.to_csv("output/phylo_tip_provenance.csv", index=False)
open("output/phylo_failed_species.txt", "w").write("\n".join(failed))

nfb = int((P.fetched_as == "genus_fallback").sum()) if len(P) else 0
print(f"\nretrieved {len(seqs)} of {len(species_list)}")
print(f"  species-level records : {len(seqs) - nfb}")
print(f"  genus fallbacks       : {nfb}   <- not that species; see the README")
print(f"  not found             : {len(failed)}")
if len(P):
    dup = P[P.duplicated("accession", keep=False)].sort_values("accession")
    if len(dup):
        print(f"\n  WARNING: {dup.accession.nunique()} accessions are reused across "
              f"{len(dup)} tips - identical sequences will collapse to zero-length "
              f"branches:")
        for a, g in dup.groupby("accession"):
            print(f"    {a}: {', '.join(g.query_species)}")
print("\nnext: bash phylogeny/02_align_trim_tree.sh")
