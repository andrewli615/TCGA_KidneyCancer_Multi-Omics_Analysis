"""Build an offline GO biological-process reference from official source files.

Usage:
  python3 scripts/build_go_reference.py \
    --hgnc hgnc_complete_set.txt \
    --gaf HUMAN-uniprot.gaf.gz \
    --obo go-basic.obo \
    --out reference

The committed reference files are sufficient to rerun the R analysis; this
script records how they were derived from the HGNC and GO source snapshots.
"""

import argparse
import csv
import gzip
import io
from collections import defaultdict
from functools import lru_cache
from pathlib import Path


def write_tsv_gz(path, rows, header):
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)


def load_hgnc(path):
    ensembl_to_symbols = defaultdict(set)
    approved_symbols = set()
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] != "Approved" or not row["symbol"]:
                continue
            approved_symbols.add(row["symbol"])
            ensembl = row["ensembl_gene_id"]
            if ensembl:
                ensembl_to_symbols[ensembl].add(row["symbol"])
    unique_map = {
        ensembl: next(iter(symbols))
        for ensembl, symbols in ensembl_to_symbols.items()
        if len(symbols) == 1
    }
    return unique_map, approved_symbols


def load_obo(path):
    terms = {}
    current = None

    def save():
        if current and current.get("id") and current.get("namespace") == "biological_process" and not current.get("obsolete"):
            terms[current["id"]] = current

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "[Term]":
                save()
                current = {"parents": []}
            elif line.startswith("["):
                save()
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("namespace: "):
                    current["namespace"] = line[11:]
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
                elif line.startswith("is_a: "):
                    current["parents"].append(line[6:].split()[0])
                elif line.startswith("relationship: part_of "):
                    current["parents"].append(line.split()[2])
    save()
    return terms


def build_annotations(gaf_path, approved_symbols, terms):
    direct = set()
    with gzip.open(gaf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[8] != "P":
                continue
            symbol, qualifier, term_id = fields[2], fields[3], fields[4]
            if symbol in approved_symbols and term_id in terms and "NOT" not in qualifier.split("|"):
                direct.add((symbol, term_id))

    @lru_cache(maxsize=None)
    def ancestors(term_id):
        result = {term_id}
        for parent in terms[term_id]["parents"]:
            if parent in terms:
                result.update(ancestors(parent))
        return frozenset(result)

    term_to_symbols = defaultdict(set)
    for symbol, term_id in direct:
        for ancestor in ancestors(term_id):
            term_to_symbols[ancestor].add(symbol)
    return term_to_symbols


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hgnc", type=Path, required=True)
    parser.add_argument("--gaf", type=Path, required=True)
    parser.add_argument("--obo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    mapping, approved = load_hgnc(args.hgnc)
    terms = load_obo(args.obo)
    annotations = build_annotations(args.gaf, approved, terms)

    write_tsv_gz(
        args.out / "ensembl_to_hgnc.tsv.gz",
        sorted(mapping.items()),
        ("ensembl_gene_id", "symbol"),
    )
    write_tsv_gz(
        args.out / "go_bp_annotations.tsv.gz",
        ((term_id, symbol) for term_id in sorted(annotations) for symbol in sorted(annotations[term_id])),
        ("go_id", "symbol"),
    )
    with (args.out / "go_bp_terms.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("go_id", "term"))
        writer.writerows((term_id, terms[term_id]["name"]) for term_id in sorted(annotations))

    print(f"Mapped {len(mapping)} Ensembl IDs; {len(annotations)} GO BP terms")


if __name__ == "__main__":
    main()
