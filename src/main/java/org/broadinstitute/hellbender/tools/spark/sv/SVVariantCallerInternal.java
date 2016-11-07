package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections4.IterableUtils;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_3_TO_5;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_5_TO_3;

/**
 * Internal caller for calling structural variants.
 */
class SVVariantCallerInternal implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * First step in calling variants: parse all alignment records for a single locally-assembled contig and generate
     * chimeric alignments if available.
     *
     * @param contigSequence                the assembled sequence
     * @param alignmentRegionsIterable      alignment records given by aligner for the sequence
     * @param minAlignLength                minimum alignment length that will be used in filtering alignments
     * @return                              the chimeric alignments of this sequence (empty if the sequence does not have any alignments)
     */
    static Iterable<ChimericAlignment> assembleBreakpointsFromAlignmentRegions(final byte[] contigSequence, final Iterable<AlignmentRegion> alignmentRegionsIterable, final Integer minAlignLength) {
        final List<AlignmentRegion> alignmentRegions = IterableUtils.toList(alignmentRegionsIterable);
        if (alignmentRegions.size() > 1) { // todo: should remove this if checking and always sort when the input is guaranteed to be of size > 1
            alignmentRegions.sort(Comparator.comparing(a -> a.startInAssembledContig));
        }
        return getBreakpointAlignmentsFromAlignmentRegions(contigSequence, alignmentRegions, minAlignLength);
    }

    // TODO: this is looking at one chimeric alignment, i.e. a pair of alignment records at a time, potentially making calling complex variant difficult
    /**
     * Second step in calling variants: produce a single {@link BreakpointAllele} from a single {@link ChimericAlignment} of a contig.
     *
     * @param breakpointIdAndAssembledBreakpoint    the inner Tuple2 is a pair of (asmID, contigID) for the contig that generated this chimeric alignment
     */
    static Tuple2<BreakpointAllele, Tuple2<Tuple2<String,String>, ChimericAlignment>> keyByBreakpointAllele(final Tuple2<Tuple2<String, String>, ChimericAlignment> breakpointIdAndAssembledBreakpoint) {
        final ChimericAlignment chimericAlignment = breakpointIdAndAssembledBreakpoint._2();
        // TODO: 11/4/16 INSDEL calling should add another kind of BreakpointAllele here
        final BreakpointAllele breakpointAllele = SVVariantCallerUtils.involvesStrandSwitch(chimericAlignment) ? new BreakpointAlleleInversion(chimericAlignment) : new BreakpointAllele(chimericAlignment);
        return new Tuple2<>(breakpointAllele, new Tuple2<>(breakpointIdAndAssembledBreakpoint._1, breakpointIdAndAssembledBreakpoint._2));
    }

    /**
     * Third step in calling variants: produce a VC from a {@link BreakpointAllele} (consensus among different assemblies if they all point to the same breakpoint).
     *
     * @param assembledBreakpointsPerAllele     consensus among different assemblies if they all point to the same breakpoint
     * @param broadcastReference                broadcasted reference
     * @throws IOException                      due to read operations on the reference
     */
    static VariantContext getVariantFromBreakpointAlleleAlignments(final Tuple2<BreakpointAllele, Iterable<Tuple2<Tuple2<String, String>, ChimericAlignment>>> assembledBreakpointsPerAllele,
                                                                   final Broadcast<ReferenceMultiSource> broadcastReference) throws IOException {

        final BreakpointAllele breakpointAllele = assembledBreakpointsPerAllele._1;
        final String contig = breakpointAllele.leftAlignedLeftBreakpoint.getContig();
        final int start = breakpointAllele.leftAlignedLeftBreakpoint.getStart();
        final int end = breakpointAllele.leftAlignedRightBreakpoint.getStart();
        VariantContextBuilder vcBuilder = new VariantContextBuilder().chr(contig).start(start).stop(end);

        vcBuilder = vcBuilder.alleles(produceAlleles(broadcastReference.getValue(), contig, start, end))
                             .id(produceVariantId(breakpointAllele));
        vcBuilder = updateAttributes(vcBuilder, breakpointAllele, assembledBreakpointsPerAllele._2, start, end);

        return vcBuilder.make();
    }

    // -----------------------------------------------------------------------------------------------
    // Group 1: from contig alignments to breakpoints
    // -----------------------------------------------------------------------------------------------

    @VisibleForTesting
    static List<ChimericAlignment> getBreakpointAlignmentsFromAlignmentRegions(final byte[] sequence, final List<AlignmentRegion> alignmentRegionList, final Integer minAlignLength) {
        if (alignmentRegionList.isEmpty()) {
            return new ArrayList<>();
        }
        final List<ChimericAlignment> results = new ArrayList<>(alignmentRegionList.size() - 1);
        final Iterator<AlignmentRegion> iterator = alignmentRegionList.iterator();
        final List<String> insertionAlignmentRegions = new ArrayList<>();
        if ( iterator.hasNext() ) {
            AlignmentRegion current = iterator.next();
            while (treatAlignmentRegionAsInsertion(current) && iterator.hasNext()) {
                current = iterator.next();
            }
            while ( iterator.hasNext() ) {
                final AlignmentRegion next = iterator.next();
                if (currentAlignmentRegionIsTooSmall(current, next, minAlignLength)) {
                    continue;
                }

                if (treatNextAlignmentRegionInPairAsInsertion(current, next, minAlignLength)) {
                    if (iterator.hasNext()) {
                        insertionAlignmentRegions.add(next.toPackedString());
                        continue;
                    } else {
                        break;
                    }
                }

                final AlignmentRegion previous = current;
                current = next;

                final byte[] sequenceCopy = Arrays.copyOf(sequence, sequence.length);

                String homology = getHomology(current, previous, sequenceCopy);
                String insertedSequence = getInsertedSequence(current, previous, sequenceCopy);

                final ChimericAlignment chimericAlignment = new ChimericAlignment(previous, current, insertedSequence, homology, insertionAlignmentRegions);

                results.add(chimericAlignment);
            }
        }
        return results;
    }
    // TODO: 11/2/16 test
    static boolean currentAlignmentRegionIsTooSmall(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return current.referenceInterval.size() - current.overlapOnContig(next) < minAlignLength;
    }
    // TODO: 11/2/16 test
    static boolean treatAlignmentRegionAsInsertion(final AlignmentRegion next) {
        return next.mapQual < 60;
    }
    @VisibleForTesting
    static boolean treatNextAlignmentRegionInPairAsInsertion(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return treatAlignmentRegionAsInsertion(next) ||
                (next.referenceInterval.size() - current.overlapOnContig(next) < minAlignLength) ||
                current.referenceInterval.contains(next.referenceInterval) ||
                next.referenceInterval.contains(current.referenceInterval);
    }
    // TODO: 11/2/16 test
    static String getHomology(final AlignmentRegion current, final AlignmentRegion previous, final byte[] sequenceCopy) {
        String homology = "";
        if (previous.endInAssembledContig >= current.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(sequenceCopy, current.startInAssembledContig - 1, previous.endInAssembledContig);
            if (previous.referenceInterval.getStart() > current.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
            }
            homology = new String(homologyBytes);
        }
        return homology;
    }
    // TODO: 11/2/16 test
    static String getInsertedSequence(final AlignmentRegion current, final AlignmentRegion previous, final byte[] sequenceCopy) {
        String insertedSequence = "";
        if (previous.endInAssembledContig < current.startInAssembledContig - 1) {

            final int insertionStart;
            final int insertionEnd;

            insertionStart = previous.endInAssembledContig + 1;
            insertionEnd = current.startInAssembledContig - 1;

            final byte[] insertedSequenceBytes = Arrays.copyOfRange(sequenceCopy, insertionStart - 1, insertionEnd);
            if (previous.referenceInterval.getStart() > current.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
            }
            insertedSequence = new String(insertedSequenceBytes);
        }
        return insertedSequence;
    }

    // -----------------------------------------------------------------------------------------------
    // Group 2: from consensus (novel adjacency/disjoint) alleles to standard VC
    // -----------------------------------------------------------------------------------------------

    /**
     * Auxiliary struct for extracting information from {@link BreakpointAllele.BreakpointAlleleInversion} to make the call.
     */
    static final class AuxStruct implements Serializable {
        private static final long serialVersionUID = 1L;

        int numAssembledBreakpoints;
        int highMqMappings;
        int maxAlignLength;

        final List<Integer> mqs;
        final List<Integer> alignLengths;
        final List<String> breakpointIds;
        final List<String> assembledContigIds;
        final List<String> insertionMappings;

        AuxStruct(final List<Tuple2<Tuple2<String, String>, ChimericAlignment>> assembledBreakpoints){
            numAssembledBreakpoints = 0;
            highMqMappings = 0;
            maxAlignLength = 0;

            final int numBreakpoints = assembledBreakpoints.size();

            mqs = new ArrayList<>(numBreakpoints);
            alignLengths = new ArrayList<>(numBreakpoints);
            breakpointIds = new ArrayList<>(numBreakpoints);
            assembledContigIds = new ArrayList<>(numBreakpoints);
            insertionMappings = new ArrayList<>(numBreakpoints);

            collectInfo(assembledBreakpoints);
        }

        private void collectInfo(final List<Tuple2<Tuple2<String, String>, ChimericAlignment>> assembledBreakpoints){
            for (final Tuple2<Tuple2<String,String>, ChimericAlignment> assembledBreakpointPair : assembledBreakpoints) {
                final ChimericAlignment chimericAlignment = assembledBreakpointPair._2;
                numAssembledBreakpoints++;
                final int assembledBreakpointMapq = Math.min(chimericAlignment.region1.mapQual, chimericAlignment.region2.mapQual);
                if (assembledBreakpointMapq == 60) {
                    highMqMappings++;
                }
                mqs.add(assembledBreakpointMapq);
                final int assembledBreakpointAlignmentLength =
                        Math.min(chimericAlignment.region1.referenceInterval.size(),
                                chimericAlignment.region2.referenceInterval.size()) - chimericAlignment.region1.overlapOnContig(chimericAlignment.region2);
                alignLengths.add(assembledBreakpointAlignmentLength);
                maxAlignLength = Math.max(maxAlignLength, assembledBreakpointAlignmentLength);
                breakpointIds.add(assembledBreakpointPair._1._1);
                assembledContigIds.add(chimericAlignment.contigId);
                insertionMappings.addAll(chimericAlignment.insertionMappings);
            }
        }
    }

    static List<Allele> produceAlleles(final ReferenceMultiSource reference, final String contig, final int start, final int end) throws IOException {
        return new ArrayList<>(Arrays.asList(Allele.create(new String(reference.getReferenceBases(null, new SimpleInterval(contig, start, start)).getBases()), true), Allele.create("<INV>")));
    }

    static String produceVariantId(final BreakpointAllele breakpointAllele) {
        final BreakpointAllele.BreakpointAlleleInversion inversionAllele = (BreakpointAllele.BreakpointAlleleInversion) breakpointAllele;
        return inversionAllele.getInversionType().name() + "_" + inversionAllele.leftAlignedLeftBreakpoint.getContig() + "_" + inversionAllele.leftAlignedLeftBreakpoint.getStart() + "_" + inversionAllele.leftAlignedRightBreakpoint.getStart();
    }

    static VariantContextBuilder updateAttributes(VariantContextBuilder vcBuilder, final BreakpointAllele breakpointAllele,
                                           final Iterable<Tuple2<Tuple2<String, String>, ChimericAlignment>> alignments,
                                           final int start, final int end){

        final BreakpointAllele.BreakpointAlleleInversion inversionAllele = (BreakpointAllele.BreakpointAlleleInversion) breakpointAllele;
        final AuxStruct auxStruct = new AuxStruct(IterableUtils.toList(alignments));

        vcBuilder = vcBuilder.attribute(VCFConstants.END_KEY, end)
                .attribute(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.INV.toString())
                .attribute(GATKSVVCFHeaderLines.SVLEN, end - start)
                .attribute(GATKSVVCFHeaderLines.TOTAL_MAPPINGS, auxStruct.numAssembledBreakpoints)
                .attribute(GATKSVVCFHeaderLines.HQ_MAPPINGS, auxStruct.highMqMappings)
                .attribute(GATKSVVCFHeaderLines.MAPPING_QUALITIES, auxStruct.mqs.stream().map(String::valueOf).sorted().collect(Collectors.joining(",")))
                .attribute(GATKSVVCFHeaderLines.ALIGN_LENGTHS, auxStruct.alignLengths.stream().map(String::valueOf).sorted().collect(Collectors.joining(",")))
                .attribute(GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH, auxStruct.maxAlignLength)
                .attribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS, auxStruct.breakpointIds.stream().sorted().collect(Collectors.joining(",")))
                .attribute(GATKSVVCFHeaderLines.CONTIG_IDS, auxStruct.assembledContigIds.stream().map(s -> s.replace(" ", "_")).sorted().collect(Collectors.joining(",")));

        if (inversionAllele.insertedSequence.length() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, inversionAllele.insertedSequence);
        }

        if (auxStruct.insertionMappings.size() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE_MAPPINGS, auxStruct.insertionMappings.stream().sorted().collect(Collectors.joining("|")));
        }

        if (inversionAllele.homology.length() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.HOMOLOGY, inversionAllele.homology);
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.HOMOLOGY_LENGTH, inversionAllele.homology.length());
        }

        if (inversionAllele.getInversionType() == INV_5_TO_3) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_5_TO_3, "");
        }

        if (inversionAllele.getInversionType() == INV_3_TO_5) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_3_TO_5, "");
        }

        return vcBuilder;
    }
}
