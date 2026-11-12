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
