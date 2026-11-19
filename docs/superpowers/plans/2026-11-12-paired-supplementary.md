# Paired-mode Supplementary/Secondary Handling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** In `--paired` mode, stop UMICollapse from writing R2 twice (and from orphaning supplementary records) when the input has supplementary or secondary alignments. Fixes upstream issue #38.

**Architecture:** Supplementary and secondary records are excluded from dedup grouping in every paired read loop. `Writer` records the read name of each primary R1 it writes. Its final full pass writes each primary R2 at most once, plus every mapped supplementary or secondary record whose name is in that set. Spec: `docs/superpowers/specs/2026-11-12-paired-supplementary-design.md`.

**Tech Stack:** Java 11, htsjdk 2.19.0, plain `javac`/`jar` via `build.sh`. Tests are `main` programs run by `run.sh`/`test.sh` (no JUnit).

**Conventions used below:**
- All commands run from the repo root `/Users/anoronha/git/umicollapse_fork`.
- `CR` means `~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork`, which runs a command inside the project's conda env. Write the full form in each command; shell state does not persist between commands.
- `test/` and `lib/` are gitignored. Scratch files (example BAM, baseline jar and outputs) go in `test/`.

---

## File map

| File | Change |
|---|---|
| `src/umicollapse/main/DeduplicateSAM.java` | Skip non-primary records in 4 paired read loops; new stat line; `Writer` survivor set and final-pass emission; `ReversedRead.equals` compares `coord` |
| `src/test/TestPairedSupplementary.java` | New: builds a tiny indexed BAM, runs `Main` in default / `--two-pass` / `--tag`, checks per-read-name record counts and output sort order |
| `test.sh` | Run the new test |
| `umicollapse.jar` | Rebuilt (repo convention: jar committed with source) |

---

### Task 1: Dev environment and baselines

**Files:** none tracked. Creates `lib/` and `test/` (both gitignored) and a conda env.

- [ ] **Step 1: Create the conda env**

```bash
~/miniforge/bin/conda create -y -n umicollapse -c conda-forge -c bioconda openjdk=11 samtools
```
Expected: the env is created. If the sandbox blocks network access, retry with the hosts the `<sandbox_violations>` block names (e.g. `conda.anaconda.org`).

- [ ] **Step 2: Download dependency jars**

```bash
mkdir -p lib && curl -sSL -o lib/htsjdk-2.19.0.jar https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar && curl -sSL -o lib/snappy-java-1.1.7.3.jar https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar && ls -l lib
```
Expected: two non-empty jars.

- [ ] **Step 3: Save the current (pre-fix) jar as a baseline**

```bash
mkdir -p test && cp umicollapse.jar test/umicollapse_orig.jar
```

- [ ] **Step 4: Build and run the existing tests**

```bash
~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./build.sh && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./test.sh 2>&1 | grep -cE "^Passed" ; ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./test.sh 2>&1 | grep -E "Error|Exception|Unrecognized" | head
```
Expected: `Done!` from the build, a non-zero `Passed` count, and no Error lines.

If `test.sh` prints `Unrecognized option: --release`, the launcher in `run.sh` is rejecting a `javac`-only flag. Change `run.sh` to:
```bash
java -cp umicollapse.jar $@
```
Then rerun this step and commit the change separately: `git add run.sh && git commit -m "Remove javac-only --release flag from run.sh"`.

- [ ] **Step 5: Restore the committed jar**

`build.sh` overwrote `umicollapse.jar`. Rebuilding from unchanged source should give a working jar, but keep the tree clean until Task 5:
```bash
git checkout -- umicollapse.jar && git status --short
```
Expected: only `?? CLAUDE.md`.

- [ ] **Step 6: Single-end baseline on UMI-tools example data**

```bash
curl -sSL -o test/example.bam https://github.com/CGATOxford/UMI-tools/releases/download/1.0.0/example.bam && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork samtools index test/example.bam && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork java -cp test/umicollapse_orig.jar:lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.7.3.jar umicollapse.main.Main bam -i test/example.bam -o test/dedup_orig.bam | tail -3
```
Expected: ends with `UMI collapsing finished in ... seconds!`, and `test/dedup_orig.bam` exists.

---

### Task 2: Failing test that reproduces issue #38

**Files:**
- Create: `src/test/TestPairedSupplementary.java`
- Modify: `test.sh`

- [ ] **Step 1: Write the test**

Create `src/test/TestPairedSupplementary.java`:

```java
package test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.util.Map;
import java.util.TreeMap;
import java.util.Arrays;

import java.io.File;
import java.nio.file.Files;

import umicollapse.main.Main;

// paired-end mode with supplementary/secondary alignments (upstream issue #38)
public class TestPairedSupplementary{
    private static final String SEQ = "\tACGTACGTAC\tIIIIIIIIII";

    // QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN
    private static final String[] RECORDS = {
        // pair A: supplementary R1 on another reference, mate fields point at the primary R2
        "pairA_AAAA\t99\tchr1\t100\t60\t10M\t=\t300\t210" + SEQ,
        "pairA_AAAA\t147\tchr1\t300\t60\t10M\t=\t100\t-210" + SEQ,
        "pairA_AAAA\t2147\tchr2\t500\t60\t10M\tchr1\t300\t0" + SEQ,
        // pair B: same position, TLEN and UMI as pair A with lower MAPQ, so it is collapsed into pair A
        "pairB_AAAA\t99\tchr1\t100\t30\t10M\t=\t300\t210" + SEQ,
        "pairB_AAAA\t147\tchr1\t300\t30\t10M\t=\t100\t-210" + SEQ,
        "pairB_AAAA\t2147\tchr2\t700\t30\t10M\tchr1\t300\t0" + SEQ,
        // pair C: secondary R1
        "pairC_CCCC\t99\tchr1\t1000\t60\t10M\t=\t1200\t210" + SEQ,
        "pairC_CCCC\t147\tchr1\t1200\t60\t10M\t=\t1000\t-210" + SEQ,
        "pairC_CCCC\t355\tchr2\t2000\t0\t10M\tchr1\t1200\t0" + SEQ,
        // pair D: supplementary R2
        "pairD_GGGG\t99\tchr1\t3000\t60\t10M\t=\t3200\t210" + SEQ,
        "pairD_GGGG\t147\tchr1\t3200\t60\t10M\t=\t3000\t-210" + SEQ,
        "pairD_GGGG\t2195\tchr2\t4000\t60\t10M\tchr1\t3000\t0" + SEQ
    };

    public static void main(String[] args) throws Exception{
        File dir = Files.createTempDirectory("umicollapse_paired_supp").toFile();
        File in = writeInput(new File(dir, "in.bam"));

        // counts per read name: primary R1, primary R2, non-primary R1, non-primary R2
        Map<String, String> deduped = new TreeMap<>();
        deduped.put("pairA_AAAA", "p1=1 p2=1 s1=1 s2=0");
        deduped.put("pairC_CCCC", "p1=1 p2=1 s1=1 s2=0");
        deduped.put("pairD_GGGG", "p1=1 p2=1 s1=0 s2=1");
        // pairB lost to pairA, so none of its records may be written

        // --tag writes every primary R1, so every pair and its non-primary records are kept
        Map<String, String> tagged = new TreeMap<>(deduped);
        tagged.put("pairB_AAAA", "p1=1 p2=1 s1=1 s2=0");

        check("default", in, dir, new String[0], deduped);
        check("two-pass", in, dir, new String[]{"--two-pass"}, deduped);
        check("tag", in, dir, new String[]{"--tag"}, tagged);
    }

    private static File writeInput(File f){
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 10000));
        header.addSequence(new SAMSequenceRecord("chr2", 10000));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        SAMLineParser parser = new SAMLineParser(header);
        SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, f);

        for(String line : RECORDS)
            writer.addAlignment(parser.parseLine(line));

        writer.close();
        return f;
    }

    private static void check(String name, File in, File dir, String[] extra, Map<String, String> expected) throws Exception{
        File out = new File(dir, name + ".bam");
        String[] base = {"bam", "-i", in.getPath(), "-o", out.getPath(), "--paired"};
        String[] args = Arrays.copyOf(base, base.length + extra.length);
        System.arraycopy(extra, 0, args, base.length, extra.length);

        try{
            Main.main(args);
        }catch(Exception e){
            e.printStackTrace();
            System.out.println("Error: paired supplementary\t" + name + "\tthrew " + e);
            return;
        }

        Map<String, String> actual = counts(out);

        if(actual.equals(expected)){
            System.out.println("Passed: paired supplementary\t" + name);
        }else{
            System.out.println("Error: paired supplementary\t" + name);
            System.out.println("Expected\t" + expected);
            System.out.println("Actual\t" + actual);
        }

        if(isCoordinateSorted(out))
            System.out.println("Passed: output coordinate-sorted\t" + name);
        else
            System.out.println("Error: output not coordinate-sorted\t" + name);
    }

    private static Map<String, String> counts(File f) throws Exception{
        Map<String, int[]> c = new TreeMap<>();
        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(f);

        for(SAMRecord record : reader){
            int[] n = c.computeIfAbsent(record.getReadName(), k -> new int[4]);
            n[(record.isSecondaryOrSupplementary() ? 2 : 0) + (record.getSecondOfPairFlag() ? 1 : 0)]++;
        }

        reader.close();

        Map<String, String> res = new TreeMap<>();

        for(Map.Entry<String, int[]> e : c.entrySet()){
            int[] n = e.getValue();
            res.put(e.getKey(), "p1=" + n[0] + " p2=" + n[1] + " s1=" + n[2] + " s2=" + n[3]);
        }

        return res;
    }

    private static boolean isCoordinateSorted(File f) throws Exception{
        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(f);
        int prevRef = -1;
        int prevStart = -1;
        boolean sorted = true;

        for(SAMRecord record : reader){
            int ref = record.getReferenceIndex();
            int start = record.getAlignmentStart();

            if(ref < prevRef || (ref == prevRef && start < prevStart))
                sorted = false;

            prevRef = ref;
            prevStart = start;
        }

        reader.close();
        return sorted;
    }
}
```

- [ ] **Step 2: Add it to `test.sh`**

Append this line to `test.sh`:
```bash
./run.sh test.TestPairedSupplementary
```

- [ ] **Step 3: Build and confirm the test fails on the current code**

```bash
~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./build.sh && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./run.sh test.TestPairedSupplementary 2>&1 | grep -E "Passed|Error|Expected|Actual|Exception"
```
Expected: `Error: paired supplementary` for `default`, `two-pass` and `tag`. In the default `Actual` line, `pairA_AAAA` shows `p2=2` (the duplicated R2 from the issue) and `pairB_AAAA` is present (the orphaned supplementary and its R2).

If the test fails to compile instead (for example an htsjdk 2.19 signature differs from the one used above), fix the test code only, then re-run.

The coordinate-sorted lines are expected to be `Passed`. If any reports `Error: output not coordinate-sorted`, the spec's assumption that htsjdk re-sorts the output is wrong. Stop and report this before continuing.

- [ ] **Step 4: Commit the test (without the rebuilt jar)**

```bash
git checkout -- umicollapse.jar && git add src/test/TestPairedSupplementary.java test.sh && git commit -m "Add failing test for paired-mode supplementary records (#38)"
```

---

### Task 3: Exclude secondary/supplementary records from dedup grouping

**Files:**
- Modify: `src/umicollapse/main/DeduplicateSAM.java` (read loops at ~lines 46–56, ~193, ~283–292, ~373; stat prints at ~258 and ~457)

All four paired read loops must skip the same records. The `--two-pass` second pass detects that a group is complete by matching a record index against the first pass's `latest`, so the two passes must filter identically.

- [ ] **Step 1: `deduplicateAndMerge` read loop**

Replace:
```java
        int chimeric = 0;
        int readCount = 0;

        for(SAMRecord record : reader){
            // always skip the reversed read
            if(paired && record.getReadPairedFlag() && record.getSecondOfPairFlag())
                continue;

            totalReadCount++;

```
with:
```java
        int chimeric = 0;
        int nonPrimary = 0;
        int readCount = 0;

        for(SAMRecord record : reader){
            // always skip the reversed read
            if(paired && record.getReadPairedFlag() && record.getSecondOfPairFlag())
                continue;

            totalReadCount++;

            // written by Writer only if the primary pair is written
            if(paired && record.isSecondaryOrSupplementary()){
                nonPrimary++;
                continue;
            }

```

- [ ] **Step 2: `deduplicateAndMergeTwoPass` first pass**

Replace:
```java
        int chimeric = 0;
        int readCount = 0;

        // first pass to figure out where each alignment position ends
        for(SAMRecord record : firstPass){
            // always skip the reversed read
            if(paired && record.getReadPairedFlag() && record.getSecondOfPairFlag())
                continue;

            totalReadCount++;

```
with:
```java
        int chimeric = 0;
        int nonPrimary = 0;
        int readCount = 0;

        // first pass to figure out where each alignment position ends
        for(SAMRecord record : firstPass){
            // always skip the reversed read
            if(paired && record.getReadPairedFlag() && record.getSecondOfPairFlag())
                continue;

            totalReadCount++;

            // written by Writer only if the primary pair is written
            if(paired && record.isSecondaryOrSupplementary()){
                nonPrimary++;
                continue;
            }

```

- [ ] **Step 3: The two second-pass filters (`--tag` pass and `--two-pass` second pass)**

There are two occurrences of this condition, with different indentation. In both, replace the first line:
```java
if(paired && ((removeUnpaired && !record.getReadPairedFlag()) // discard unpaired
```
with these two lines, keeping the existing indentation and aligning the `||` with the lines below it:
```java
if(paired && (record.isSecondaryOrSupplementary() // written by Writer if the primary pair is written
            || (removeUnpaired && !record.getReadPairedFlag()) // discard unpaired
```
After the edit, the `--tag` pass (in `deduplicateAndMerge`) reads:
```java
                if(paired && (record.isSecondaryOrSupplementary() // written by Writer if the primary pair is written
                            || (removeUnpaired && !record.getReadPairedFlag()) // discard unpaired
                            || (record.getReadPairedFlag() && record.getSecondOfPairFlag()) // ignore reversed reads
                            || (record.getReadPairedFlag() && record.getMateUnmappedFlag()) // discard unmapped reads
                            || (removeChimeric && record.getReadPairedFlag()
                                && !record.getReferenceName().equals(record.getMateReferenceName())))){ // discard chimeric reads
                    continue;
                }
```
and the `--two-pass` second pass reads:
```java
            if(paired && (record.isSecondaryOrSupplementary() // written by Writer if the primary pair is written
                        || (removeUnpaired && !record.getReadPairedFlag()) // discard unpaired
                        || (record.getReadPairedFlag() && record.getSecondOfPairFlag()) // ignore reversed reads
                        || (record.getReadPairedFlag() && record.getMateUnmappedFlag()) // discard unmapped reads
                        || (removeChimeric && record.getReadPairedFlag()
                            && !record.getReferenceName().equals(record.getMateReferenceName())))){ // discard chimeric reads
                continue;
            }
```

- [ ] **Step 4: Report the new count (both methods)**

There are two occurrences of:
```java
            System.out.println("Number of chimeric reads\t" + chimeric);
```
In both, add this line directly after it:
```java
            System.out.println("Number of secondary/supplementary reads\t" + nonPrimary);
```

- [ ] **Step 5: Build and run the test**

```bash
~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./build.sh && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./run.sh test.TestPairedSupplementary 2>&1 | grep -E "Passed|Error|Expected|Actual|Exception"
```
Expected: no exceptions, and still `Error` for all three modes. In `Actual`, every pair now has `p2=1` and `pairB_AAAA` is absent in default and two-pass. All `s1`/`s2` are 0, because nothing writes non-primary records yet (Task 4). The program prints `Number of secondary/supplementary reads\t4`.

- [ ] **Step 6: Commit (without the rebuilt jar)**

```bash
git checkout -- umicollapse.jar && git add src/umicollapse/main/DeduplicateSAM.java && git commit -m "Exclude secondary/supplementary records from paired-mode dedup grouping (#38)"
```

---

### Task 4: `Writer` emits non-primary records for written pairs

**Files:**
- Modify: `src/umicollapse/main/DeduplicateSAM.java` (`ReversedRead.equals` ~line 480; `Writer` ~lines 521–607)

- [ ] **Step 1: `ReversedRead.equals` also compares `coord` (matching `hashCode`)**

Replace:
```java
            if(ref != a.ref)
                return false;

            if(!name.equals(a.name))
                return false;
```
with:
```java
            if(ref != a.ref)
                return false;

            if(coord != a.coord)
                return false;

            if(!name.equals(a.name))
                return false;
```

- [ ] **Step 2: Survivor set field and construction**

Replace:
```java
        private String ref = null;
        private HashSet<ReversedRead> set;

        public Writer(File in, File out, SamReader r, boolean paired){
            if(paired){
                this.in = in;
                this.set = new HashSet<ReversedRead>();
            }
```
with:
```java
        private String ref = null;
        private HashSet<ReversedRead> set;
        private HashSet<String> survivors; // names of written pairs, for their secondary/supplementary records

        public Writer(File in, File out, SamReader r, boolean paired){
            if(paired){
                this.in = in;
                this.set = new HashSet<ReversedRead>();
                this.survivors = new HashSet<String>();
            }
```

- [ ] **Step 3: Record survivors in `write`**

Replace:
```java
                if(record.getReadPairedFlag()){
                    set.add(new ReversedRead(
                            record.getReadName(),
                            record.getMateReferenceName(),
                            record.getMateAlignmentStart()
                    ));
                }
            }

            writer.addAlignment(record);
```
with:
```java
                if(record.getReadPairedFlag()){
                    set.add(new ReversedRead(
                            record.getReadName(),
                            record.getMateReferenceName(),
                            record.getMateAlignmentStart()
                    ));
                }

                survivors.add(record.getReadName());
            }

            writer.addAlignment(record);
```

- [ ] **Step 4: Emit non-primary records in the final full pass; match only primary R2s**

In `writeReversed`, replace:
```java
            while(iter.hasNext()){
                SAMRecord record = iter.next();

                if(!record.getReadUnmappedFlag()
```
with:
```java
            while(iter.hasNext()){
                SAMRecord record = iter.next();

                if(record.isSecondaryOrSupplementary()){
                    // only the full pass, so each one is written once
                    if(fullPass && !record.getReadUnmappedFlag() && survivors.contains(record.getReadName()))
                        writer.addAlignment(record);

                    continue;
                }

                if(!record.getReadUnmappedFlag()
```
This `continue` also stops a secondary or supplementary R2 from matching a pending primary-R2 entry.

- [ ] **Step 5: Build and run the test**

```bash
~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./build.sh && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./run.sh test.TestPairedSupplementary 2>&1 | grep -E "Passed|Error|Expected|Actual|Exception"
```
Expected: exactly these six lines, with no `Error`:
```
Passed: paired supplementary	default
Passed: output coordinate-sorted	default
Passed: paired supplementary	two-pass
Passed: output coordinate-sorted	two-pass
Passed: paired supplementary	tag
Passed: output coordinate-sorted	tag
```

- [ ] **Step 6: Commit (without the rebuilt jar)**

```bash
git checkout -- umicollapse.jar && git add src/umicollapse/main/DeduplicateSAM.java && git commit -m "Write secondary/supplementary records once, only for written pairs (#38)"
```

---

### Task 5: Regression checks and rebuilt jar

**Files:**
- Modify: `umicollapse.jar` (rebuilt)

- [ ] **Step 1: Rebuild and run the full test suite**

```bash
~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./build.sh && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./test.sh 2>&1 | grep -E "Error|Exception" ; echo "exit-grep=$?"
```
Expected: no output from grep, then `exit-grep=1` (meaning no Error or Exception lines).

- [ ] **Step 2: Single-end output unchanged on example.bam**

```bash
~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./umicollapse bam -i test/example.bam -o test/dedup_new.bam | tail -1 && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork bash -c 'diff <(samtools view test/dedup_orig.bam) <(samtools view test/dedup_new.bam) && echo IDENTICAL' && ~/miniforge/bin/conda run -n umicollapse --cwd /Users/anoronha/git/umicollapse_fork ./run.sh test.CompareDedupUMI test/dedup_orig.bam test/dedup_new.bam _ | tail -1
```
Expected: `IDENTICAL`, then `Wrong	0`.

The `./umicollapse` launcher requests a 12–15 GB heap (`-Xms12G -Xmx15G`). If the JVM can't reserve that much on this machine, run the same thing with `java -Xmx2G -jar umicollapse.jar bam -i test/example.bam -o test/dedup_new.bam` instead.

- [ ] **Step 3: Commit the rebuilt jar**

```bash
git add umicollapse.jar && git status --short && git commit -m "Rebuild jar with paired-mode supplementary fix (#38)"
```
Expected: `git status --short` shows `M  umicollapse.jar` staged and `?? CLAUDE.md` untracked (not committed).

---

## Optional follow-up (user's call, not part of this plan's acceptance)

- Real-data check on Iris through `hpcgate`: run the new jar with `--paired` on the BAM from the issue, then run Picard `SamToFastq` on the output and confirm there is no "Illegal mate state" error. Never ssh to Iris directly.
- `--tag` follow-up: give secondary/supplementary records their primary's `MI`/`RX` tags and duplicate flag.
