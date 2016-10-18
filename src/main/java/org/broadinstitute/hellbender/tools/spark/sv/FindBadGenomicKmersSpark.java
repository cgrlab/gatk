package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * SparkTool to identify 63-mers in the reference that occur more than 3 times.
 */
@CommandLineProgramProperties(summary="Find the set of high copy number kmers in a reference.",
        oneLineSummary="find ref kmers with high copy number",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class FindBadGenomicKmersSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    @VisibleForTesting static final Long MAX_KMER_FREQ = 3L;
    private static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    private static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;

    @Argument(doc = "file for ubiquitous kmer output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "kmer size.", fullName = "kSize", optional = true)
    private int kSize = SVConstants.KMER_SIZE;

    @Argument(doc = "maximum kmer DUST score", fullName = "kmerMaxDUSTScore")
    private int maxDUSTScore = SVConstants.MAX_DUST_SCORE;

    @Argument(doc = "additional high copy kmers (mitochondrion, e.g.) fasta file name",
            fullName = "highCopyFasta", optional = true)
    private String highCopyFastaFilename;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /** Get the list of high copy number kmers in the reference, and write them to a file. */
    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final SAMFileHeader hdr = getHeaderForReads();
        SAMSequenceDictionary dict = null;
        if ( hdr != null ) dict = hdr.getSequenceDictionary();
        final PipelineOptions options = getAuthenticatedGCSOptions();
        final ReferenceMultiSource referenceMultiSource = getReference();
        Collection<SVKmer> killList = findBadGenomicKmers(ctx, kSize, maxDUSTScore, referenceMultiSource, options, dict);
        if ( highCopyFastaFilename != null ) {
            killList = uniquify(killList, processFasta(kSize, maxDUSTScore, highCopyFastaFilename, options));
        }
        SVUtils.writeKmersFile(kSize, outputFile, options, killList);
    }

    /** Find high copy number kmers in the reference sequence */
    public static List<SVKmer> findBadGenomicKmers( final JavaSparkContext ctx,
                                                    final int kSize,
                                                    final int maxDUSTScore,
                                                    final ReferenceMultiSource ref,
                                                    final PipelineOptions options,
                                                    final SAMSequenceDictionary readsDict ) {
        // Generate reference sequence RDD.
        final JavaRDD<byte[]> refRDD = getRefRDD(ctx, kSize, ref, options, readsDict);

        // Find the high copy number kmers
        return processRefRDD(kSize, maxDUSTScore, refRDD);
    }

    /**
     * Turn a text file of overlapping records from a reference sequence into an RDD, and do a classic map/reduce:
     * Kmerize, mapping to a pair <kmer,1>, reduce by summing values by key, filter out <kmer,N> where
     * N <= MAX_KMER_FREQ, and collect the high frequency kmers back in the driver.
     */
    @VisibleForTesting static List<SVKmer> processRefRDD( final int kSize,
                                                          final int maxDUSTScore,
                                                          final JavaRDD<byte[]> refRDD ) {
        return refRDD
                .flatMapToPair(seq ->
                        SVKmerizerWithLowComplexityFilter.stream(seq, kSize, maxDUSTScore)
                                .map(kmer -> new Tuple2<>(kmer.canonical(kSize), 1))
                                .collect(Collectors.toCollection(() -> new ArrayList<>(seq.length))))
                .reduceByKey(Integer::sum)
                .filter(kv -> kv._2 > MAX_KMER_FREQ)
                .map(kv -> kv._1)
                .collect();
    }

    /**
     * Create an RDD from the reference sequences.
     * The reference sequences are transformed into a single, large collection of byte arrays. The collection is then
     * parallelized into an RDD.
     * Each contig that exceeds a size given by REF_RECORD_LEN is broken into a series of REF_RECORD_LEN chunks with a
     * K-1 base overlap between successive chunks. (I.e., for K=63, the last 62 bases in chunk n match the first 62
     * bases in chunk n+1) so that we don't miss any kmers due to the chunking -- we can just kmerize each record
     * independently.
     */
    private static JavaRDD<byte[]> getRefRDD( final JavaSparkContext ctx,
                                              final int kSize,
                                              final ReferenceMultiSource ref,
                                              final PipelineOptions options,
                                              final SAMSequenceDictionary readsDict ) {
        final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(readsDict);
        if ( dict == null ) throw new GATKException("No reference dictionary available");

        final int effectiveRecLen = REF_RECORD_LEN - kSize + 1;
        final List<byte[]> sequenceChunks = new ArrayList<>();
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final String seqName = rec.getSequenceName();
            final int seqLen = rec.getSequenceLength();
            final SimpleInterval interval = new SimpleInterval(seqName, 1, seqLen);
            try {
                final byte[] bases = ref.getReferenceBases(options, interval).getBases();
                for ( int start = 0; start < seqLen; start += effectiveRecLen ) {
                    sequenceChunks.add(Arrays.copyOfRange(bases, start, Math.min(start+REF_RECORD_LEN, seqLen)));
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't get reference sequence bases for " + interval, ioe);
            }
        }

        return ctx.parallelize(sequenceChunks, sequenceChunks.size()/REF_RECORDS_PER_PARTITION+1);
    }

    @VisibleForTesting static List<SVKmer> processFasta( final int kSize,
                                                         final int maxDUSTScore,
                                                         final String fastaFilename,
                                                         final PipelineOptions options ) {
        try ( BufferedReader rdr = new BufferedReader(new InputStreamReader(BucketUtils.openFile(fastaFilename, options))) ) {
            List<SVKmer> kmers = new ArrayList<>((int)BucketUtils.fileSize(fastaFilename, options));
            String line;
            StringBuilder sb = new StringBuilder();
            while ( (line = rdr.readLine()) != null ) {
                if ( line.charAt(0) != '>' ) sb.append(line);
                else if ( sb.length() > 0 ) {
                    SVKmerizerWithLowComplexityFilter.stream(sb,kSize,maxDUSTScore).forEach(kmers::add);
                    sb.setLength(0);
                }
            }
            if ( sb.length() > 0 ) {
                SVKmerizerWithLowComplexityFilter.stream(sb, kSize, maxDUSTScore).forEach(kmers::add);
            }
            return kmers;
        }
        catch ( IOException ioe ) {
            throw new GATKException("Can't read high copy kmers fasta file "+fastaFilename, ioe);
        }
    }

    private static Collection<SVKmer> uniquify(final Collection<SVKmer> coll1, final Collection<SVKmer> coll2 ) {
        final HopscotchSet<SVKmer> kmers = new HopscotchSet<>(coll1.size() + coll2.size());
        kmers.addAll(coll1);
        kmers.addAll(coll2);
        return kmers;
    }
}
