# Test fixtures

| file | what | source and licence | written by |
|---|---|---|---|
| `bodymap_counts_subset.csv.gz` | raw counts, genes × samples, every gene, of 4 rat BodyMap samples from 21-week-old animals (liver, brain, testes, thymus) | Yu et al., *Nature Communications* 5:3230 (2014), GEO GSE53960, via Bioconductor `bodymapRat` 1.28.0; **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/); technical runs summed, subset to 4 samples | `scripts/40_product_validation.py` (gzip with mtime 0, so reruns are byte-identical) |
| `product_parity.json` | the pipeline's log2 CPM (`src/tfp/io.py` `log_cpm`, total library) of the 20 panel genes for those 4 samples, and their library sizes | derived from the file above (same terms) | `scripts/40_product_validation.py` |

They let `tests/test_check_core.js` check the browser's count → log2 CPM conversion against Python on real counts.
`tools/check_flow.js` uses the same counts for the drag-and-drop check. See `NOTICE.md` for all third-party data.
