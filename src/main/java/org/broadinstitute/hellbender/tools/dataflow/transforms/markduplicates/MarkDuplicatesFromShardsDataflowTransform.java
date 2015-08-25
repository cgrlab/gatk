package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.coders.CoderException;
import com.google.cloud.dataflow.sdk.coders.CustomCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.coders.StringUtf8Coder;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.Flatten;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.Partition;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionList;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.dataflow.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import javax.ejb.Local;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * The main trasform for MarkDuplicates. Takes reads, returns them with some reads marked as duplicates.
 */
public final class MarkDuplicatesFromShardsDataflowTransform extends PTransform<PCollection<KV<String, Iterable<GATKRead>>>, PCollection<GATKRead>> {
    private static final long serialVersionUID = 1l;


    public class LocalCollection<T> {
        ArrayList<T> data;

        public LocalCollection() {
            data = new ArrayList<>();
        }

        // don't change ts after passing it in.
        public LocalCollection(ArrayList<T> ts) {
            data = ts;
        }

        public ArrayList<LocalCollection<T>> partitionBy(ToIntFunction<T> keyer) {
            ArrayList<LocalCollection<T>> ret = new ArrayList<LocalCollection<T>>(2);
            for (T el : data) {
                int key = keyer.applyAsInt(el);
                while (ret.size() < key) {
                    ret.add(new LocalCollection<T>());
                }
                // we're mutating the LocalCollection but that's OK because it's not closed yet.
                ret.get(key).data.add(el);
            }
            // no more mutation from this point on.
            return ret;
        }

        public LocalGroupedCollection<T> groupBy(Function<T,String> keyer) {
            return new LocalGroupedCollection<>(data,keyer);
        }

        // apply a function to every element of the collection
        public <U> LocalCollection<U> map(Function<T,U> func) {
            ArrayList<U> results = new ArrayList<U>(data.size());
            for (T el : data) {
                results.add(func.apply(el));
            }
            return new LocalCollection<>(results);
        }

        // apply a function to the whole collection
        public <U> LocalCollection<U> transform(Function<Iterable<T>, ArrayList<U>> func) {
            return new LocalCollection<>(func.apply(data));
        }

        public Iterable iterable() {
            return data;
        }
    }

    public class LocalGroupedCollection<T> {
        HashMap<String, LocalCollection<T>> groups;

        public LocalGroupedCollection(ArrayList<T> data, Function<T,String> keyer) {
            groups = new HashMap<>();
            for (T el : data) {
                String key = keyer.apply(el);
                LocalCollection<T> c;
                c = groups.getOrDefault(key, null);
                if (c==null) {
                    c = new LocalCollection<>();
                    groups.put(key, c);
                }
                c.data.add(el);
            }
        }

        private LocalGroupedCollection(HashMap<String, LocalCollection<T>> data) {
            groups = data;
        }

        // get a particular group
        public LocalCollection<T> get(String key) {
            return groups.get(key);
        }

        // apply a function to each group
        public <U> LocalGroupedCollection<U> map(Function<LocalCollection<T>, LocalCollection<U>> mapper) {
            HashMap<String, LocalCollection<U>> groups2 = new HashMap<String, LocalCollection<U>>();
            for (String k : groups.keySet()) {
                groups2.put(k, mapper.apply(groups.get(k)) );
            }
            return new LocalGroupedCollection(groups2);
        }

        // map and transform together
        public <U> LocalGroupedCollection<U> mapTransform(Function<Iterable<T>, ArrayList<U>> transformFn) {
            return map( (collection) -> collection.transform(transformFn) );
        }

        // return a single collection that contains the union of all the groups' content.
        public LocalCollection<T> flatten() {
            ArrayList<T> data = new ArrayList<>();
            for (LocalCollection<T> c : groups.values()) {
                data.addAll(c.data);
            }
            return new LocalCollection<>(data);
        }
    }


    private void newApply(final SAMFileHeader header, Iterable<GATKRead> readsIterable) {

        ArrayList<GATKRead> readsArray = Lists.newArrayList(readsIterable);
        LocalCollection<GATKRead> reads = new LocalCollection<>(readsArray);

        ArrayList<LocalCollection<GATKRead>> readsPartitioned = reads.partitionBy(MarkDuplicates::primaryNonPrimaryAlignment);
        LocalCollection<GATKRead> fragments = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        LocalCollection<GATKRead> fragmentsTransformed = transformFragments(header, fragments);

        LocalCollection<GATKRead> pairs = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        //LocalCollection<GATKRead> pairsTransformed = ...

//
//    final PCollection<GATKRead> pairs = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
//    final PCollection<GATKRead> pairsTransformed = MarkDuplicatesUtils.transformReads(header, opticalDuplicateFinder, pairs);
//
//    //no work on those
//    final PCollection<GATKRead> not_primary = readsPartitioned.get(ReadsPartition.NOT_PRIMARY.ordinal());
//
//    return PCollectionList.of(fragmentsTransformed).and(pairsTransformed).and(not_primary).apply(Flatten.<GATKRead>pCollections());
    }

    private LocalCollection<GATKRead> transformFragments(final SAMFileHeader header, LocalCollection<GATKRead> fragments) {
        return fragments
            .map(MarkDuplicatesUtils::markNonDuplicate)
            .groupBy((read) -> MarkDuplicatesUtils.keyForFragment(header, read))
            .mapTransform(MarkDuplicatesUtils::transformFragmentsFn)
            .flatten();
    }

    /**
     * (1) keyReadsByName: label each read with its read group and read name.
     * (2) GroupByKey: group together reads with the same group and name.
     * (3) keyPairedEndsWithAlignmentInfo:
     *   (a) Sort each group of reads (see GATKOrder below).
     *   (b) Pair consecutive reads into PairedEnds. In most cases there will only be two reads
     *       with the same name. TODO: explain why there might be more.
     *   (c) Label each read with alignment information: Library, reference index,
     *       stranded unclipped start and reverse strand.
     *   (d) Leftover reads are emitted, unmodified, as an unpaired end.
     * (4) GroupByKey: Group PairedEnds that share alignment information. These pairs
     *     are duplicates of each other.
     * (5) markDuplicatePairs:
     *   (a) For each group created by (4), sort the pairs by score and mark all but the
     *       highest scoring as duplicates.
     *   (b) Determine which duplicates are optical duplicates and increase the overall count.
     */
    static LocalCollection<GATKRead> transformReads(final SAMFileHeader header, final OpticalDuplicateFinder finder, final LocalCollection<GATKRead> pairs) {

        pairs
            // group reads by genomic location etc.
            .groupBy((read) -> MarkDuplicatesUtils.keyForRead(header, read))
            // use the groups to build paired reads
            .mapTransform( (reads) -> MarkDuplicatesUtils.pairTheReads(header,reads) )
            // merge the groups together
            .flatten()
            // now group by pair key
            .groupBy((pair) -> MarkDuplicatesUtils.keyForPair(header, pair));

        return null;

/*
        final PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> keyReadsByName = keyReadsByName(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> keyPairedEndsWithAlignmentInfo =
            keyPairedEndsWithAlignmentInfo(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markDuplicatePairedEnds = markDuplicatePairedEnds(finderPcolView);

        return pairs
            .apply(keyReadsByName)
            .apply(GroupByKey.<String, GATKRead>create())
            .apply(keyPairedEndsWithAlignmentInfo)
            .setCoder(KvCoder.of(StringUtf8Coder.of(), new PairedEndsCoder()))
            .apply(GroupByKey.<String, PairedEnds>create())
            .apply(markDuplicatePairedEnds);
        */
    }

    // primary vs non-primary alignments
    private enum ReadsPartition {
        PRIMARY, NOT_PRIMARY
    }

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    private static final int MIN_BASE_QUAL = 15;


    /**
     * Struct-like class to store information about the paired reads for mark duplicates.
     */
    private static final class PairedEnds {
        private GATKRead first, second;

        PairedEnds(final GATKRead first) {
            this.first = first;
        }

        public static PairedEnds of(final GATKRead first) {
            return new PairedEnds(first);
        }

        public PairedEnds and(final GATKRead second) {
            if (second != null &&
                ReadUtils.getStrandedUnclippedStart(first) > ReadUtils.getStrandedUnclippedStart(second)) {

                this.second = this.first;
                this.first = second;
            } else {
                this.second = second;
            }
            return this;
        }

        public String key(final SAMFileHeader header) {
            throw new RuntimeException("in progress");
            //return MarkDuplicatesReadsKey.keyForPairedEnds(header, first, second);
        }

        public GATKRead first() {
            return first;
        }

        public GATKRead second() {
            return second;
        }

        public int score() {
            return MarkDuplicatesFromShardsDataflowTransform.score(first) + MarkDuplicatesFromShardsDataflowTransform.score(second);
        }
    }

    /**
     * Special coder for the PairedEnds class.
     */
    private static final class PairedEndsCoder extends CustomCoder<PairedEnds> {
        private static final long serialVersionUID = 1L;

        private static final CustomCoder<GATKRead> readCoder = new GATKReadCoder();

        @Override
        public void encode( PairedEnds value, OutputStream outStream, Context context ) throws CoderException, IOException {
            if ( value == null || value.first() == null ) {
                throw new IOException("nothing to encode");
            }
            final boolean isCompletePair = value.second() != null ;
            SerializableCoder.of(Boolean.class).encode(isCompletePair, outStream, context);

            readCoder.encode(value.first, outStream, context);
            if ( isCompletePair ) {
                readCoder.encode(value.second, outStream, context);
            }
        }

        @Override
        public PairedEnds decode( InputStream inStream, Context context ) throws CoderException, IOException {
            final boolean isCompletePair = SerializableCoder.of(Boolean.class).decode(inStream, context);

            final PairedEnds pairedEnds = PairedEnds.of(readCoder.decode(inStream, context));
            if ( isCompletePair ) {
                pairedEnds.and(readCoder.decode(inStream, context));
            }
            return pairedEnds;
        }
    }

    // TODO: extract this into an independent class, and unify with other comparators in the codebase
    private static final class CoordinateOrder implements Comparator<GATKRead>, Serializable {
        private static final long serialVersionUID = 1l;
        private final SAMFileHeader header;

        public CoordinateOrder(final SAMFileHeader header) {
            this.header = header;
        }

        @Override
        public int compare(final GATKRead lhs, final GATKRead rhs) {
            if (rhs == lhs) return 0; //shortcut

            final int res1 = Integer.compare(ReadUtils.getReferenceIndex(lhs, header), ReadUtils.getReferenceIndex(rhs, header));
            if (res1 != 0) return res1;

            final int res2 = Long.compare(lhs.getStart(), rhs.getStart());
            if (res2 != 0) return res2;

            final int res3 = Boolean.compare(lhs.isDuplicate(), rhs.isDuplicate());
            if (res3 != 0) return res3;

            final int res4 = Boolean.compare(lhs.failsVendorQualityCheck(), rhs.failsVendorQualityCheck());
            if (res4 != 0) return res4;

            final int res5 = Boolean.compare(lhs.isPaired(), rhs.isPaired());
            if (res5 != 0) return res5;

            final int res6 = Boolean.compare(lhs.isProperlyPaired(), rhs.isProperlyPaired());
            if (res6 != 0) return res6;

            final int res7 = Boolean.compare(lhs.isFirstOfPair(), rhs.isFirstOfPair());
            if (res7 != 0) return res7;

            final int res8 = Boolean.compare(lhs.isSecondaryAlignment(), rhs.isSecondaryAlignment());
            if (res8 != 0) return res8;

            final int res9 = Boolean.compare(lhs.isSupplementaryAlignment(), rhs.isSupplementaryAlignment());
            if (res9 != 0) return res9;

            final int res10 = Integer.compare(lhs.getMappingQuality(), rhs.getMappingQuality());
            if (res10 != 0) return res10;

            final int res11 = Integer.compare(ReadUtils.getMateReferenceIndex(lhs, header), ReadUtils.getMateReferenceIndex(rhs, header));
            if (res11 != 0) return res11;

            final int res12 = Long.compare(lhs.getMateStart(), rhs.getMateStart());
            return res12 ;
        }
    }

    /**
     *  The header for all reads processed.
     */
    private final PCollectionView<SAMFileHeader> header;

    public MarkDuplicatesFromShardsDataflowTransform(final PCollectionView<SAMFileHeader> header) {
        this.header = header;
    }

    @Override
    public PCollection<GATKRead> apply( final PCollection<KV<String, Iterable<GATKRead>>> preads ) {
        return preads.apply(ParDo.named("MarkDuplicatesFromShardsDataflow")
            .withSideInputs(header)
            .of(new DoFnWLog<KV<String, Iterable<GATKRead>>, GATKRead>("MarkDuplicatesFromShardsDataflow") {
                private static final long serialVersionUID = 1l;

                @Override
                public void processElement(ProcessContext c) throws Exception {
                    HashMap<String, List<GATKRead>> pairsTogether = new HashMap<String, List<GATKRead>>();
                    HashMap<String, List<PairedEnds>> pairByKey = new HashMap<String, List<PairedEnds>>();
                    HashMap<String, List<GATKRead>> fragmentByKey = new HashMap<String, List<GATKRead>>();
                    final SAMFileHeader h = c.sideInput(header);
                    String shard = c.element().getKey();
                    for (GATKRead r : c.element().getValue()) {
                        if (r instanceof SAMRecordToGATKReadAdapter) {
                            SAMRecordToGATKReadAdapter sam = (SAMRecordToGATKReadAdapter)r;
                            sam.setHeader(h);
                        }
                        if (!isPrimary(r)) {
                            // no work on those
                            c.output(r);
                            continue;
                        }
                        // mutating the input! That's only OK because we know there's no sibling transform to this one
                        r.setIsDuplicate(false);
                        final String key = MarkDuplicatesUtils.keyForRead(h,r);
                        if (null!=key) addToList(fragmentByKey, key, r);

                        if (ReadUtils.readHasMappedMate(r)) {
                            throw new RuntimeException("in progress");
                            //final String key2 = MarkDuplicatesReadsKey.keyForPair(h,r);
                            //addToList(pairsTogether, key2, r);
                        }
                    }
                    bunny.stepEnd("first grouping");
                    // now we have all the data, process it

                    // transform fragments
                    for (List<GATKRead> colocatedFragments : fragmentByKey.values()) {

                        // "transform fragment r"

                        //split reads by paired vs unpaired
                        final Map<Boolean, List<GATKRead>> byPairing = StreamSupport.stream(colocatedFragments.spliterator(), false).collect(Collectors.partitioningBy(
                            read -> ReadUtils.readHasMappedMate(read)
                        ));

                        // Note the we emit only fragments from this mapper.
                        if (byPairing.get(true).isEmpty()) {
                            // There are no paired reads, mark all but the highest scoring fragment as duplicate.
                            final List<GATKRead> frags = Ordering.natural().reverse().onResultOf((GATKRead read) -> score(read)).immutableSortedCopy(byPairing.get(false));
                            if (!frags.isEmpty()) {
                                c.output(frags.get(0));                         //highest score - just emit
                                for (final GATKRead record : Iterables.skip(frags, 1)) {  //lower   scores - mark as dups and emit
                                    record.setIsDuplicate(true);
                                    c.output(record);
                                }
                            }
                        } else {
                            // There are paired ends so we mark all fragments as duplicates.
                            for (final GATKRead record : byPairing.get(false)) {
                                record.setIsDuplicate(true);
                                c.output(record);
                            }
                        }
                    }
                    bunny.stepEnd("transform fragments");

                    // transform reads
                    // markGroupedDuplicatePairs
                    for (List<GATKRead> pairs : pairsTogether.values()) {
                        pairs.sort(new CoordinateOrder(h));
                        PairedEnds pair = null;
                        //Records are sorted, we iterate over them and pair them up.
                        for (final GATKRead record : pairs) {
                            if (pair == null) {                                //first in pair
                                pair = PairedEnds.of(record);
                            } else {                                           //second in pair
                                pair.and(record);
                                addToList(pairByKey, pair.key(h), pair);
                                pair = null;                                   //back to first
                            }
                        }
                        if (pair != null) {                                    //left over read
                            addToList(pairByKey, pair.key(h), pair);
                        }
                    }
                    bunny.stepEnd("transform reads (1/2)");
                    // markPairedEnds
                    for (List<PairedEnds> pairs : pairByKey.values()) {
                        final ImmutableListMultimap<Boolean, PairedEnds> paired = Multimaps.index(pairs, pair -> pair.second() != null);

                        // As in Picard, unpaired ends left alone.
                        for (final PairedEnds pair : paired.get(false)) {
                            c.output(pair.first());
                        }

                        //order by score
                        final List<PairedEnds> scored = Ordering.natural().reverse().onResultOf((PairedEnds pair) -> pair.score()).immutableSortedCopy(paired.get(true));
                        final PairedEnds best = Iterables.getFirst(scored, null);
                        if (best != null) {
                            c.output(best.first());
                            c.output(best.second());
                        }
                        //Mark everyone who's not best as a duplicate
                        for (final PairedEnds pair : Iterables.skip(scored, 1)) {
                            GATKRead record = pair.first();
                            record.setIsDuplicate(true);
                            c.output(record);

                            record = pair.second();
                            record.setIsDuplicate(true);
                            c.output(record);
                        }
                    }
                    bunny.stepEnd("transform reads (2/2)");

                }
            }));
    }

    public static <T> void addToList(HashMap<String,List<T>> map, String key, T value) {
        if (map.containsKey(key)) {
            map.get(key).add(value);
        }
        else {
            ArrayList<T> list = new ArrayList<T>();
            list.add(value);
            map.put(key, list);
        }
    }

    private static boolean isPrimary(GATKRead read) {
        return !(read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped());
    }

    // lots of below is unneeded, will remove once I'm sure we're doing exactly the same thing


    /**
     * (1) Reads are grouped by read group and read name.
     * (2) Reads are then paired together as follows:
     *   (a) The remaining reads (one per fragment key) are coordinate-sorted and paired
     *       consecutively together and emitted.
     *   (b) If a read is leftover, it is emitted, unmodified, as an unpaired end.
     * (3) Paired ends are grouped by a similar key as above but using both reads.
     * (4) Any unpaired end is emitted, unmodified.
     * (5) The remained paired ends are scored and all but the highest scoring are marked as
     *     duplicates. Both reads in the pair are emitted.
     */
    private PCollection<GATKRead> transformReads(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<GATKRead> pairs) {
        final PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForPairs = makeKeysForPairs(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> markGroupedDuplicatePairs = markGroupedDuplicatePairs(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markPairedEnds = markPairedEnds();

        return pairs
            .apply(makeKeysForPairs)
            .apply(GroupByKey.<String, GATKRead>create())
            .apply(markGroupedDuplicatePairs)
            .setCoder(KvCoder.of(StringUtf8Coder.of(), new PairedEndsCoder()))
            .apply(GroupByKey.<String, PairedEnds>create())
            .apply(markPairedEnds);
    }

    /**
     * Makes keys for read pairs. To be grouped by in the next step.
     */
    private PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForPairs(final PCollectionView<SAMFileHeader> headerPcolView) {
        return ParDo
            .named("make keys for pairs")
            .withSideInputs(headerPcolView)
            .of(new DoFn<GATKRead, KV<String, GATKRead>>() {
                private static final long serialVersionUID = 1l;

                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final GATKRead record = context.element();
                    if (ReadUtils.readHasMappedMate(record)) {
                        final SAMFileHeader h = context.sideInput(headerPcolView);
                        throw new RuntimeException("in progress");
                        /*
                        final String key = ""; //MarkDuplicatesReadsKey.keyForPair(h, record);
                        final KV<String, GATKRead> kv = KV.of(key, record);
                        context.output(kv);
                        */
                    }
                }
            });
    }

    private PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> markGroupedDuplicatePairs(final PCollectionView<SAMFileHeader> headerPcolView) {
        return ParDo
            .named("pair ends")
            .withSideInputs(headerPcolView)
            .of(new DoFn<KV<String, Iterable<GATKRead>>, KV<String, PairedEnds>>() {
                private static final long serialVersionUID = 1L;
                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final SAMFileHeader header = context.sideInput(headerPcolView);
                    final List<GATKRead> sorted = Lists.newArrayList(context.element().getValue());
                    sorted.sort(new CoordinateOrder(header));
                    PairedEnds pair = null;
                    //Records are sorted, we iterate over them and pair them up.
                    for (final GATKRead record : sorted) {
                        if (pair == null) {                                //first in pair
                            pair = PairedEnds.of(record);
                        } else {                                           //second in pair
                            pair.and(record);
                            context.output(KV.of(pair.key(header), pair));
                            pair = null;                                   //back to first
                        }
                    }
                    if (pair != null) {                                    //left over read
                        context.output(KV.of(pair.key(header), pair));
                    }
                }
            });
    }


    private PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markPairedEnds() {
        return ParDo
            .named("mark paired ends")
            .of(new DoFn<KV<String, Iterable<PairedEnds>>, GATKRead>() {
                private static final long serialVersionUID = 1l;

                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final ImmutableListMultimap<Boolean, PairedEnds> paired = Multimaps.index(context.element().getValue(), pair -> pair.second() != null);

                    // As in Picard, unpaired ends left alone.
                    for (final PairedEnds pair : paired.get(false)) {
                        context.output(pair.first());
                    }

                    //order by score
                    final List<PairedEnds> scored = Ordering.natural().reverse().onResultOf((PairedEnds pair) -> pair.score()).immutableSortedCopy(paired.get(true));
                    final PairedEnds best = Iterables.getFirst(scored, null);
                    if (best != null) {
                        context.output(best.first());
                        context.output(best.second());
                    }
                    //Mark everyone who's not best as a duplicate
                    for (final PairedEnds pair : Iterables.skip(scored, 1)) {
                        GATKRead record = pair.first();
                        record.setIsDuplicate(true);
                        context.output(record);

                        record = pair.second();
                        record.setIsDuplicate(true);
                        context.output(record);
                    }
                }
            });
    }


    /**
     * Takes the reads,
     * group them by library, contig, position and orientation,
     * within each group
     *   (a) if there are only fragments, mark all but the highest scoring as duplicates, or,
     *   (b) if at least one is marked as paired, mark all fragments as duplicates.
     *  Note: Emit only the fragments, as the paired reads are handled separately.
     */
    PCollection<GATKRead> transformFragments(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<GATKRead> fragments) {
        final PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForFragments =  makeKeysForFragments(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<GATKRead>> markGroupedDuplicateFragments = markGroupedDuplicateFragments();
        return fragments
            .apply(makeKeysForFragments)
            .apply(GroupByKey.<String, GATKRead>create())
            .apply(markGroupedDuplicateFragments);//no need to set up coder for Read (uses GenericJsonCoder)
    }

    private PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<GATKRead>> markGroupedDuplicateFragments() {
        return ParDo.named("mark dups")
            .of(new DoFn<KV<String, Iterable<GATKRead>>, GATKRead>() {
                private static final long serialVersionUID = 1l;

                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    //split reads by paired vs unpaired
                    final Map<Boolean, List<GATKRead>> byPairing = StreamSupport.stream(context.element().getValue().spliterator(), false).collect(Collectors.partitioningBy(
                        read -> ReadUtils.readHasMappedMate(read)
                    ));

                    // Note the we emit only fragments from this mapper.
                    if (byPairing.get(true).isEmpty()) {
                        // There are no paired reads, mark all but the highest scoring fragment as duplicate.
                        final List<GATKRead> frags = Ordering.natural().reverse().onResultOf((GATKRead read) -> score(read)).immutableSortedCopy(byPairing.get(false));
                        if (!frags.isEmpty()) {
                            context.output(frags.get(0));                         //highest score - just emit
                            for (final GATKRead record : Iterables.skip(frags, 1)) {  //lower   scores - mark as dups and emit
                                record.setIsDuplicate(true);
                                context.output(record);
                            }
                        }
                    } else {
                        // There are paired ends so we mark all fragments as duplicates.
                        for (final GATKRead record : byPairing.get(false)) {
                            record.setIsDuplicate(true);
                            context.output(record);
                        }
                    }
                }
            });
    }

    /**
     * How to assign a score to the read in MarkDuplicates (so that we pick the best one to be the non-duplicate).
     */
    //Note: copied from htsjdk.samtools.DuplicateScoringStrategy
    private static int score(final GATKRead record) {
        if (record == null) {
            return 0;
        } else {
            int sum = 0;
            for ( byte b : record.getBaseQualities() ) {
                int i = (int)b;
                if ( i >= MIN_BASE_QUAL ) {
                    sum += i;
                }
            }
            return sum;
        }
    }


    /**
     * Groups reads by keys - keys are tuples of (library, contig, position, orientation).
     */
    PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForFragments(final PCollectionView<SAMFileHeader> headerPcolView) {
        return ParDo
            .named("make keys for reads")
            .withSideInputs(headerPcolView)
            .of(new DoFn<GATKRead, KV<String, GATKRead>>() {
                private static final long serialVersionUID = 1L;
                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final GATKRead record = context.element();
                    record.setIsDuplicate(false);
                    final SAMFileHeader h = context.sideInput(headerPcolView);
                    throw new RuntimeException("in progress");
                    /*
                    final String key = MarkDuplicatesReadsKey.keyForFragment(h, record);
                    final KV<String, GATKRead> kv = KV.of(key, record);
                    context.output(kv);
                    */
                }
            });
    }
}