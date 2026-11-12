# Paired mode: supplementary/secondary records follow their primary

Upstream issue: https://github.com/Daniel-Liu-c0deb0t/UMICollapse/issues/38

## Problem

In `--paired` mode, when the input BAM contains supplementary alignments, R2 is written twice, and Picard `SamToFastq` fails with "Illegal mate state".

Root cause, in `src/umicollapse/main/DeduplicateSAM.java`:

1. Only R2 records (`getSecondOfPairFlag`) are skipped when building dedup groups. Supplementary (0x800) and secondary (0x100) R1 records are treated as pair representatives, each grouped under its own coordinate and TLEN.
2. When a representative survives, `Writer.write` adds `ReversedRead(name, mateRef, mateStart)` to a pending set. A supplementary R1's RNEXT/PNEXT points at the same primary R2, so it adds the same key.
3. Groups are written by reference. On each reference change, `writeReversed(false)` writes matching primary R2s and removes them from the set. If the supplementary R1 is on a different reference from the primary R1 (the chimeric case), the primary's R2 has already been written and removed by the time the supplementary is processed. The supplementary re-adds it, and the final `writeReversed(true)` in `close()` writes R2 again.

A second consequence of the same cause: the primary and supplementary R1 are deduplicated independently, so a supplementary can survive while its primary is removed as a duplicate. This leaves an R2 whose only partner is a supplementary record, which is also an illegal mate state.

## Goal

In paired mode, supplementary and secondary records (from R1 or R2) are excluded from deduplication. They appear in the output exactly once if and only if their pair's primary R1 was written. Each primary R2 is written at most once.

## Scope

- In scope: the default paired path (`deduplicateAndMerge`) and `--two-pass` (`deduplicateAndMergeTwoPass`). Both share `Writer`.
- `--tag`: must not crash. Supplementary/secondary records pass through untagged (no `MI`/`RX`/`cs`/`su`, no duplicate flag). Giving them their primary's tags is a follow-up.
- Out of scope: single-end mode behaviour is unchanged.

## Design

### 1. Reading and grouping

In paired mode, skip any record with `isSecondaryOrSupplementary()` before the alignment key is built, in every loop that feeds grouping or writing:

- `deduplicateAndMerge` read loop
- the `--tag` second pass in `deduplicateAndMerge`
- the `--two-pass` first pass
- the `--two-pass` second pass

The two `--two-pass` loops must skip exactly the same records. The second pass detects that a group is complete by comparing a running record index against the `latest` index recorded in the first pass, so any mismatch in filtering would flush groups early or never.

The skip happens after `totalReadCount++`, so the reported input count is unchanged. A new counter is printed in paired mode as `Number of secondary/supplementary reads`.

### 2. `Writer` (paired mode only)

1. New field `HashSet<String> survivors`. `write(record)` adds `record.getReadName()` for every record it is given (only primary R1s, after section 1). The set is kept until `close()`.
2. `writeReversed`: the primary R2 match condition also requires `!record.isSecondaryOrSupplementary()`.
3. `ReversedRead.equals` also compares `coord`, matching `hashCode`. The match key becomes name + reference + position.
4. `writeReversed(true)` (the full pass in `close()`) handles each record one of two ways:
   - primary R2: existing logic (write if pending, then remove from the pending set);
   - secondary or supplementary record (R1 or R2), mapped, with its name in `survivors`: write it.

   The partial per-reference passes (`fullPass == false`) still write only primary R2s.

Result for each surviving pair: primary R1 written once during grouping, primary R2 once by the R2 logic, and each supplementary/secondary record once in the single final full pass. A pair whose primary R1 was removed as a duplicate contributes no records.

### Output ordering

Emission order does not determine file order. `Writer` opens the output with `makeSAMOrBAMWriter(header, false, out)`, which is not presorted, so htsjdk re-sorts to the header's sort order. Paired mode needs a coordinate-sorted, indexed input because `writeReversed` uses `reader.query`, so the output is coordinate-sorted on close. This is to be confirmed during testing.

### Memory

One name string per surviving pair. This is new for `--two-pass`, where it is roughly the size of the pending-R2 set at its peak but kept until the end. If it becomes a problem, storing a hash of the name would reduce it.

## Testing

Environment: a conda env (from `~/miniforge`) with `openjdk=11` and `samtools`. Download `htsjdk-2.19.0.jar` and `snappy-java-1.1.7.3.jar` into `lib/` as described in the README.

New test `src/test/TestPairedSupplementary.java`, a `main` program in the existing test style and added to `test.sh`. It uses htsjdk to build a small coordinate-sorted, indexed BAM (refs `chr1` and `chr2`, `_UMI` read-name suffix):

| Pair | Setup | Expected output |
|---|---|---|
| A | Primary R1+R2 on chr1; supplementary R1 on chr2 with mate fields pointing to R2 on chr1 | primary R1 ×1, primary R2 ×1, supplementary ×1 |
| B | Same position, TLEN and UMI as A, lower MAPQ (loses dedup); supplementary on chr2 | no records |
| C | Primary pair plus a secondary R1 | primary R1 ×1, primary R2 ×1, secondary ×1 |
| D | Primary pair plus a supplementary R2 | primary R1 ×1, primary R2 ×1, R2 supplementary ×1 |

The test runs default, `--two-pass` and `--tag`, and checks per-name counts of primary R1, primary R2 and non-primary records. For `--tag` it checks that the run completes and that all primary records are present.

Test-first: run the new test against the current code and confirm it fails on pair A (duplicate R2) and pair B (orphaned supplementary).

Regression checks:
- the existing `./test.sh` suite passes;
- on UMI-tools' single-end `example.bam`, the output of the current and new jars is identical by `test.CompareDedupUMI`.

Optional: run the new jar on the original BAM on Iris through `hpcgate`, and confirm Picard `SamToFastq` no longer reports "Illegal mate state".

## Deliverables

Source changes, the new test, an updated `test.sh`, and a rebuilt `umicollapse.jar`, committed together following the repo convention.
