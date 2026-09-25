# Gene annotation snapshot

These tables let `project.Rmd` run GO biological-process enrichment offline.
They were built on 2026-09-24 with [`scripts/build_go_reference.py`](../scripts/build_go_reference.py)
from the following official downloads:

| Source | URL | SHA-256 of downloaded file |
| --- | --- | --- |
| HGNC complete set | https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt | `69bb5722d5a42bb355580deb2c9f197ce3d7f7d13807173b9674a7db65e52191` |
| Human GO annotations | https://current.geneontology.org/annotations/gaf/HUMAN-uniprot.gaf.gz | `a0afba19dfb1f8fa996bc1bdcd61fd0c9bd4cf0d2bf2509d09ac86993c0e70a2` |
| GO basic ontology | https://current.geneontology.org/ontology/go-basic.obo | `b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc` |

The builder keeps approved HGNC symbols, filters GO annotations to biological
processes without `NOT`, and propagates each annotation through `is_a` and
`part_of` parent terms. The saved files are the exact reference used for this
repository's enrichment results; current upstream downloads may differ.
