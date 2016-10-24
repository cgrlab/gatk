package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.*;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * SparkTool to identify kmers in a given reference
 */
@CommandLineProgramProperties(summary="Builds library of reference kmers to be used for in the PathSeq subtraction phase.",
        oneLineSummary="Builds library of reference kmers",
        programGroup = SparkProgramGroup.class)
public final class PathSeqKmerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    @VisibleForTesting
    private static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    private static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;

    @Argument(doc = "file for kmer output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "kmer size", fullName = "kSize", optional = true)
    private int kSize = 31;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /** Get the list of distinct kmers in the reference, and write them to a file. */
    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final SAMFileHeader hdr = getHeaderForReads();
        SAMSequenceDictionary dict = null;
        if ( hdr != null ) dict = hdr.getSequenceDictionary();
        final PipelineOptions options = getAuthenticatedGCSOptions();
        final ReferenceMultiSource referenceMultiSource = getReference();
        List<SVKmerSmall> kmerList = findKmers(ctx, kSize, referenceMultiSource, options, dict);

        HopscotchSet<SVKmerSmall> kmerSet = new HopscotchSet<>(kmerList);

        Output output = new Output(BucketUtils.createFile(outputFile, options));
        Kryo kryo=new Kryo();
        kryo.writeClassAndObject(output, kmerSet);
        output.close();
    }

    /** Find kmers in the reference sequence */
    public static List<SVKmerSmall> findKmers(final JavaSparkContext ctx,
                                              final int kSize,
                                              final ReferenceMultiSource ref,
                                              final PipelineOptions options,
                                              final SAMSequenceDictionary readsDict ) {
        // Generate reference sequence RDD.
        final JavaRDD<byte[]> refRDD = SVUtils.getRefRDD(ctx, kSize, ref, options, readsDict, REF_RECORD_LEN, REF_RECORDS_PER_PARTITION);

        // Non-redundant list of kmers
        return processRefRDD(kSize, refRDD);
    }

    /**
     * Turn a text file of overlapping records from a reference sequence into an RDD, and do a classic map/reduce:
     * Kmerize, mapping map to the collection of all canonicalized kmers, removing duplicates, and returning result to
     * the driver.
     */
    @VisibleForTesting static List<SVKmerSmall> processRefRDD( final int kSize,
                                                          final JavaRDD<byte[]> refRDD ) {
        return refRDD.flatMap(seq -> SVKmerizer.stream(seq, kSize, new SVKmerSmall(kSize))
                .map(kmer -> kmer.canonical(kSize))
                .collect(Collectors.toList())).distinct().collect();
    }

    private static Collection<SVKmer> uniquify(final Collection<SVKmer> coll1, final Collection<SVKmer> coll2 ) {
        final HopscotchSet<SVKmer> kmers = new HopscotchSet<>(coll1.size() + coll2.size());
        kmers.addAll(coll1);
        kmers.addAll(coll2);
        return kmers;
    }
}
