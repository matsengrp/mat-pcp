"""Compression policy for the PCP outputs, defined once (#134).

The four PCP families -- ``{clade}_pcps``, ``all_pcps``,
``{clade}_pcps_dasm_ready`` and ``all_pcps_dasm_ready`` -- are written as
``.csv.gz``. On ``gisaid-and-public_2026-04-22_v8`` they came to **504 GB
apparent (470 GiB), 93% of the whole 542 GB output tree**, because every row
stores a full parent sequence *and* its child, differing by a handful of
mutations, and adjacent rows within a clade are near-identical to each other.

Quote one measure or the other, never a mix: ``du --apparent-size`` gives
504/542 GB (93%), while ``du``'s allocated blocks give 403/415 GiB (97%) on this
GPFS filesystem. An earlier revision of this docstring paired the apparent
``pcps/`` figure with the allocated tree total, which cannot both hold -- 470 GiB
is not 97% of 415 GiB.

Measured gzip ratios, sampled at 25/50/75% through each aggregate file so the
numbers are not an artifact of the head of the file:

===========  =====  ===========
gene         level  ratio
===========  =====  ===========
ORF1a        9      ~125x
ORF1a        6      ~17x
ORF1a        1      ~3x
S            9      ~110x
S            6      ~100x
===========  =====  ===========

**Level 9 beats level 6 on ratio *and* wall-clock** (ORF1a: 125x at ~97 MiB/s vs
17x at 29 MiB/s), consistently across every sample. That looks backwards but is
not: at level 9 gzip's full match search finds the very long repeats -- parent
against child within a row, and adjacent near-identical rows -- so it emits far
less output and finishes sooner. At level 6 the truncated hash-chain search
misses them and spends the time encoding literals instead. Level 9 is already
pandas' default for gzip, so ``compresslevel`` is deliberately not set here.

``mtime=0`` is the part that is easy to get wrong. Two facts about pandas' gzip
output, both measured:

* the default embeds **mtime=now**, so writing the same DataFrame twice yields
  different bytes;
* pandas sets the gzip **FNAME** flag (``FLG=0x08``), so the compressed bytes
  also embed the *output filename* -- identical content under two names differs.

The decompressed bytes are identical in every case. So ``mtime=0`` makes a
re-run byte-reproducible for a given filename, which matters because
``data-central`` manifests pin ``md5``/``xxh3`` per file: without it an unchanged
re-run would look like a data change. It does **not** make two differently-named
files comparable, which is why
``validate_test.py``'s cross-implementation sentinel md5s the *decompressed*
stream rather than the file.

Reading needs nothing: ``pd.read_csv`` infers gzip from the extension, as does
``dasm2.data.registry.load_dataframe`` on the consumer side (an
extension-matched check accepting ``.csv`` and ``.csv.gz``, so the name must end
in exactly ``.csv.gz`` -- ``.gz`` alone or ``.csv.zst`` raises there).

This lives in one module rather than being repeated at each ``to_csv`` because a
writer that quietly omitted ``mtime=0`` would break reproducibility for that one
file only, which no output-level check would attribute correctly.
"""

# Passed as ``to_csv(..., compression=PCP_CSV_COMPRESSION)`` by every rule that
# writes one of the four PCP families.
PCP_CSV_COMPRESSION = {"method": "gzip", "mtime": 0}
