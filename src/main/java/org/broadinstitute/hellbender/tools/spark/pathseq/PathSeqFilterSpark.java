package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.github.lindenb.jbwa.jni.ShortRead;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.collect.Iterators;
import htsjdk.samtools.*;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerSmall;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.transformers.BaseQualityClipReadTransformer;
import org.broadinstitute.hellbender.transformers.DUSTReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.seqdoop.hadoop_bam.BAMInputFormat;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Read preprocessing and host filtering on reads from a BAM file",
        oneLineSummary = "PathSeqFilter on Spark",
        programGroup = SparkProgramGroup.class)
public final class PathSeqFilterSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;
    private SAMFileHeader header;

    private static final String maxLengthArgName = "pathseq_maxReadLength";
    @Argument(fullName = maxLengthArgName,
            shortName = maxLengthArgName,
            doc="Keep only reads with length at most equal to the specified value",
            optional=true)
    public int MAX_READ_LENGTH = 1000;

    private static final String minLengthArg = "pathseq_minReadLength";
    @Argument(fullName = minLengthArg,
            shortName = minLengthArg,
            doc="Keep only reads with length at least equal to the specified value",
            optional=true)
    public int MIN_READ_LENGTH = 31;

    @Argument(fullName="pathseq_ambigFilterFrac", shortName="pathseq_ambigFilterFrac", optional=true)
    public float FRAC_N_THRESHOLD = 0.03f;

    @Argument(fullName="pathseq_baseQualThresh", shortName="pathseq_baseQualThresh", optional=true)
    public int QUAL_PHRED_THRESH = 15;

    @Argument(fullName="pathseq_baseQualFrac", shortName="pathseq_mbaseQualFrac", optional=true)
    public float QUAL_MAX_FRAC = 0.05f;

    @Argument(fullName="pathseq_kmerLibPath", shortName="pathseq_kmerLibPath", optional=false)
    public String KMER_LIB_PATH;

    @Argument(fullName="pathseq_kmerSize", shortName="pathseq_kmerSize", optional=true)
    public int KMER_SIZE = 31;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String output;

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @Override
    @SuppressWarnings("unchecked")
    protected void runTool(final JavaSparkContext ctx) {

        header = getHeaderForReads();
        ctx.hadoopConfiguration().setBoolean(BAMInputFormat.KEEP_PAIRED_READS_TOGETHER_PROPERTY, true);

        final JavaRDD<GATKRead> reads = getReads();

        //Filter secondary/supplementary reads and reads that fail the vendor quality check
        JavaRDD<GATKRead> primaryReads = reads.filter(read -> !(read.isSecondaryAlignment() || read.failsVendorQualityCheck() || read.isSupplementaryAlignment()));

        //Filter duplicate reads
        final OpticalDuplicateFinder finder =
                new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null);
        JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.mark(primaryReads, getHeaderForReads(),MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES, finder, getRecommendedNumReducers());

        ReadFilter markedFilter = new MarkedDuplicateReadFilter();
        final JavaRDD<GATKRead> markedFilteredReads = markedReads.filter(read -> markedFilter.test(read));

        //Apply DUST masking
        ReadTransformer dustTransformer = new DUSTReadTransformer();
        final JavaRDD<GATKRead> readsDUSTMasked = markedFilteredReads.map(read -> dustTransformer.apply(read));

        //Apply base quality clipping
        ReadTransformer bqClippingTransformer = new BaseQualityClipReadTransformer();
        final JavaRDD<GATKRead> readsClipped = readsDUSTMasked.map(read -> bqClippingTransformer.apply(read));

        //Filter using base quality scores
        ReadFilter bqFilter = new BaseQualityReadFilter(QUAL_PHRED_THRESH,QUAL_MAX_FRAC);
        final JavaRDD<GATKRead> readsBQFiltered = readsClipped.filter(read -> bqFilter.test(read));

        //Filter reads with ambiguous bases
        ReadFilter ambigFilter = new AmbiguousBaseReadFilter(FRAC_N_THRESHOLD);
        final JavaRDD<GATKRead> readsAmbigFiltered = readsBQFiltered.filter(read -> ambigFilter.test(read));

        //Filter reads with less than MIN_READ_LENGTH bases
        ReadFilter readLengthFilter = new ReadLengthReadFilter(MIN_READ_LENGTH,MAX_READ_LENGTH);
        final JavaRDD<GATKRead> readsLengthFiltered = readsAmbigFiltered.filter(read -> readLengthFilter.test(read));

        //Load Kmer hopscotch set and filter reads containing > 0 matching kmer
        final PipelineOptions options = getAuthenticatedGCSOptions();
        Input input = new Input(BucketUtils.openFile(KMER_LIB_PATH, options));
        Kryo kryo=new Kryo();
        HopscotchSet<SVKmerSmall> kmerLibSet = (HopscotchSet<SVKmerSmall>)kryo.readClassAndObject(input);
        Broadcast<HopscotchSet<SVKmerSmall>> kmerLibBroadcast = ctx.broadcast(kmerLibSet);
        ContainsKmerReadFilter readKmerFilter = new ContainsKmerReadFilter(kmerLibBroadcast,KMER_SIZE);
        final JavaRDD<GATKRead> readsKmerFiltered = readsLengthFiltered.filter(read -> !readKmerFilter.test(read));

        //TODO: bwa filtering against user-specified reference

        //"Clean" reads by resetting all flags and setting pairedness flags according to the remaining read set
        final JavaRDD<GATKRead> readsCleaned = readsKmerFiltered
                .groupBy(read -> read.getName()) //group pairs using read name into set of Iterable<GATKRead>'s
                .values()   //Get Iterables
                .map(p -> {     //Determine if Iterable represents 1 or 2 reads and set flags accordingly
                    int size = 0;
                    Iterator<GATKRead> iter = p.iterator();
                    while (iter.hasNext()) {
                        size++;
                        iter.next();
                    }
                    boolean notPair = (size==1);
                    ArrayList<GATKRead> newPair = new ArrayList<>(size);
                    iter = p.iterator();
                    while (iter.hasNext()) {
                        newPair.add(CleanRead(iter.next(),notPair));
                    }
                    return newPair;
                })
                .flatMap(p -> p);   //Flatten back from Iterables to GAKTReads

        writeReads(ctx, output, readsCleaned);
    }

    private boolean BaseQualityFilterFn(GATKRead read) {
        return false;
    }

    private GATKRead CleanRead(GATKRead read, boolean notPaired) {

        GATKRead new_read = new SAMRecordToGATKReadAdapter(new SAMRecord(header));

        if (notPaired) {
            new_read.setIsPaired(false);
        } else {
            if (read.isFirstOfPair()) {
                new_read.setIsFirstOfPair();
            } else {
                new_read.setIsSecondOfPair();
            }
        }

        new_read.setName(read.getName());
        new_read.setReadGroup(read.getReadGroup());
        new_read.setBases(read.getBases());
        new_read.setBaseQualities(read.getBaseQualities());
        return new_read;
    }

}
