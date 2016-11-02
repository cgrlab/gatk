package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.List;


public class CallVariantsFromAlignedContigsSparkTest extends BaseTest {

    static final String LONG_CONTIG1 = "TTTTTTTTTTTTTTTCTGAGACCGAATCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACGATCTTGGCTTACTGCAAGCTCTGCCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCCCACCCCCCCACCCCCCCAGGTAGCTG" +
            "GGACTACAGGTGTCTGCCACCACACCTGGCTAAATTTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAAGATGGTTTCGCTCTCCTGACCTCGCGATCCGCCCACCTCGGCCTCTCAAAGTGCTGGGATTACAGGCCTGAGCCACTGCGCCC" +
            "TGCCTGACCCCTTCTTTTAAAACAAATCTTTTGTCTTTGTCTTCATTTCTGCATTCGTCCCCTTCGTTCAATCCTGTAGGAATTGACAGTGATATTGGGGAGTTCCCATCCCTGGATTTGGGATTTCCTCGAGTTTCCAGCCCTGTCCTTGTGGCCAAAAAGT" +
            "GGCACAGTTTGGTTTTTGAATCCTGTTCGTGACGCCAAAAATATTCTGTGTGGGAAACATTCAAGGGGAGAAGAAAAGACACACACACAATACCTTTAAGGGTAAATAAGCTTTGTGCCACGTAAATGGAAATGCAGATATAATTAGCAAATTATATAATAAG" +
            "CAAATCAATGTAATAACCAGATTGATATAATAAGCATGTTGATATAATAAGCAAATTGCAATGGGAAGAGGAGAAAGGAAAAGAGATATATATATTTACACTCACCAGACTATGGAGGATTCACCACCAGACTGGGAAGCAACAGCTTGGGCTCCAGAGTCAG" +
            "CTACTCATTCCTGCACAGATGAGGAGGGTCTAATGAAGCTTCAGCACAATCTGATACCCTAGCTCTTTTGTAACGAGTTGTTTGGCATAAGGCCCAGTCATGAGGGCCATTTGCAACTGGGGTCAAGGAACACAAAACTGTCAACTTGTTTTTGCGATTGTCT" +
            "ATTGTTTTTCAACAACTAATATATAGAAATAGATTGAAATAGAGATTTCTCTGAAACAGCGCTGAATGAATGCCTCAAGGGCCTCACACAACCTGTTTGAGAACTTGGTGACTACTGAGTTTGTCCACGTTCAATTAAGTTCAAATTTAGTATTTAACTTTTC" +
            "CTCCACAAATTACCCAGTCTTGATGAATCAGCTCTGTCTAGACATAGGGCAAGGTGAACCCCTTGGGCAGTTACACAACCTCCGCCTTCTGGGTTTAAGCAATTCTCCTGCCTCAGCCTCCGGACTAGCTGGGTCTACAGGTGTGCAGCACCACACCCAGCTA" +
            "GTTATTTGTACTTTTAGTAGAAATGGGGTTTCACCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAAGTGATCCACCCACCTTGGCCTCCCAAAGTGTTGCGATTACAGGCACTTGCCAGTGAACCTGGCCCTAAATGACTTCTTTCTATCTCCTATAA" +
            "CAATTTGAAATTACTTAAAGGTGGTTTCAAATTGAAAAAATAAAAAGAATTTGGATAAAAATAAAATGTAAACAGTTTTTAAAAATTACAAGAGATTACAAAATATACATGTAAAACCGGAGTGGTCAAAAATGACAAATTTGATTTATTTATAAGGTTTATT" +
            "AAAATTAGCTTTAGTATTGATAATACACTATTACCAAAGTAAAAGCTGATTTTCTCTTGAAAAAAATTTTATGTATTATTAATATGACAGCAAAATACTTCTGTTCACCTTTTGAATATATTCAAAAAGAGAGAGAGTAAAAAAAAAGATAGAATTTTCCCAT" +
            "GATCTGGGTTGGGCCTGGCTTAGCTCAGGGAGGAAGCCCTGCCTGAAAAACGCTGCAGCTTAGGCTGTGACTTTTTCTTCACTCAGTCCTGATCACATTTTCTGTTACTCAGGGCCTGAGGGGGCGGGGGCCTTAAGCATTATCCAATCAGAAACGCTGGGCT" +
            "GACGCCCCGTCCGGGAGGGAGGTGGGGGGGGTCAGCCCCCCGCCCAGCCAGCCGCCCCGTCCGGGAGGTGGGGGGTGCCTCTGCCCGGCCGCCCCTACTGGGAAGTGAGGAGCCCCTCTGCCCGGCCACCACCCCGTCTGGGAGGTGTACCCAACAGCTCATT" +
            "GAGAACGGGCCATGATGACAATGGCAGTTTTGTGGAATAGAAACGGGGGAAAGGTGGGGAAAAGATTGAGAAATCGGATGGTTGCTGTGTCTGTGTAGAAAGAGGTAGACATGGGAGACTTCATTTTGTTCTGTACTAAGAAAAATTCTTCTGCCTTGGGATC" +
            "CTGTTGACCTATGACCTTACCCCCAACCCTGTGCTCTCTGAAACATGTGCTGTGTCCACTCAGGGTTAAATGGATTAAGGGCGGTGCAAGATGTGCTTTGTTAAACAGATGCTTGAAGGCAGCATGCTCGTTAAGAGTCATCACCACTCCCTAATCTCAAGTA" +
            "CCCGGGGACACAAACACTGCGGAAGGCCGCAGGGTCCTCTGCCTAGGAAAACCAGAGACCTTTGTTCACTTGTTTATCTGCTGACCTTCCCTCCACTGTTGTCCTATGACCCTGCCAAATCCCCCTCTGCGAGAAACACCCAAGAATGATCAACTAAAAAAAA" +
            "AAAAGAAAAAAGAAAAAAGAAAAAATACACCCTGGCAGAAAACTTCTCGAATGTAAGTAAAGTGGGAGTTTTCAACTTGCCCACCTCCTTACAAAAACATGCAATAAGAATGTGCTAGATATAAAGCTTACAAATGTGAGAAATACGGGAAGTTTTTCAATTT" +
            "TAACAAATACTTTCAAAGTCATGTGAAAATTCTCACTAGAAAGAAATCCAATAAAGGCAAGTTGAGGCATAAAAAAAAAAAAAAGAAAAAGCTCAGAATACATATAAGATTGACTAAGTCAGCCAGAAAATAATCCCCTAAAAGAAATTTCTCTCTAAACACC" +
            "CAATGTGCACAGCTACTCTCAGCATGAGAAACATGAGCTTTATGAAGAAAGGTGGCAGATTTTCAGAAGAATTTTATAAAAGTTTCTTTTCCATCTCTGCTGCTCTCTCATCTCCTAGCCATTGAATGGGGGTTCTATATTGAAATACATCTGACAACTTCCA" +
            "ACAACACTTTTTGATCAAGAAATAGAATTTGACTATGTTCGTATAGTGGAATATATTAGAACTTGTAACACAGCTAACTGAATAGCTATTATGGTGTTTGGGTGGCCACATCACCTGTCTTTATTTGTCCGGTAATAGCAGCATTCCAATTTAAAGAAATAAA" +
            "AGATACCAAAATTGTGTTTACTTTTAATTATTCCTATTGAATAAAGTAATAAGCATGTCAGACTGATATCTATCATAACAATAAATTTTGTTTGGATATTATATTAGATATAAATATTTAAGTATGAATAATTTTAATGAACTAGTCATAATGTATGTAGCAT" +
            "TTTTAAAAATTGTAACTATACTTCAGTTAAAACACTTTATATTTCAAAAGCATAAATAACAATATTAAAATAACAATTTAGGTGATTCATTCAAAGTAAGTATTGGGGCTTTATGTTCATACTATTGTAGAAAATACTGCTTATGGCTCATGCCTGTAATCCC" +
            "AGCACATTGGGAGGCTGAGGTGGGTAGATCACCTGAGGTCAGGAGTTCCTGATCCCATCTTTACTAAAAATACAAAACTTACCCAGGGTGGTTGTGCACACTTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATTGCTTGAACAAGGGAGGAAATGG" +
            "TTGCAGTGAGCCATGATCATGCCACTGAACCCCAGCCTGGGCAAGAGAGTGAGACTGTCTCAAAAAAAAAAAAAACTGTTTAATTTTTATGAATGCAGGTTTTCTGCAAACACTACACATAACTATGCTAATTGTTCTGAAGTAATAAATAGAAAGCAAGGCA" +
            "CAACTACAGACTCCACTGTTCAGTTTATGCACTGAACTGTTCTTGCTTTTGCAGTGTAAGTATTTCTGCCTGCAAATACTGGATAATTACCTTGGATCATCAGATTTCTATCAAAGGAATTTAGTATCTTTTAGTCTTTATCATTTTGTATTGCTAAATTTAT" +
            "CTGTGTGTTAAGCTTCTGTGTGCTCTTAAAATGAGGTTTTATCTAAACAAACCTGTGTCTACTTTAAAAGACTAAACATGAAAAAACTAAACTTTTCAGAACCAAAAACAAAGCAATAAATCTGAAGTACTAGATAGTCTGGAGTGAGATTTATTTAGCTTTT" +
            "TTTTTTTTTTGAGATGGAGTCTCGCTCTGTCACCGAGGCTGGAGTGCAGTGGCACGAACTCGGCTCACTGCAAAAGCTCTGCCTCCCAGCTTCATGCCATTCTCCTACCTCAGCCTCCCAAGTAGCTGGGATTACAGGCAACTGCCACCACGCCCAGCTAATT" +
            "TTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTTGATCTCCTGACCTCATGATCTGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCATGCAAAACCGCGCCCTGCCCTTATTTAGCTCTTAATATGATTTACATATAT" +
            "TTCAAAAAAGCAGAGAAAAATATCTACATATAATCTAAATCCCTTAAGAGAGAATAGCAAAAAAATTTTTGGAACTTTTTAAAGAGTTTTGAACTCTTGGACATCTGAATTTTGCACACTATATGCACTTGAAAGAATGTTAATGGGGAAAAAGCAGAAGAGA" +
            "AAAAGACGTTATAAAAAAATCCATGAGTGCACAACACCTATGAAACAGAATAGAGAGCCCAGAAATAATGCCTTTCACCTGCAACCATCAGATTTCTGACAAAGCTGACAAGAGGAATGTGGGAAGAATTCTCTCTTTCATAAATGGTGCTGGAATAACTATC" +
            "TACCACTATGTAGAAGACTGAAGTGGACCCCTTCATTACACCATATAAAAAAATCAACTGAAGATAAATTAAGGACTTAAATGTAAAACTTAAAATTATAAGAAACCCTGCAAGATAACCTAGGAAATAGCATTCTAGACACAGAAACAGGTAAAGACTTCAT" +
            "GATGAAGCTACCAAAAGCAACTGCAACAGAAGTAAATTGACAAATGGGATGTATTTAAACTTAAGAGCTTCTTCACAGCAAAGGAAACTATCAACAGAGTAAACAGACAAACTAGAGAATAAAAGAATATATTTGTAAATTTTGCCTCTGAAAAAGGTCTAAT" +
            "ATACAGAATTTATTAGGAACTTAAACAAGTTTACAAGAAAAAAACACACTCATTAAAAAGTATGCAAAAAACATGAACAGATGCCTTTCAATAGAAGATACACATAGGGCTAACAAGCATATGAAAAAAAAAATGCTCATTGCTAATCATTAGAGAAATTAGA" +
            "AGAGAAATCACAGGCTGGGTGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGCAGATTACAAGGTCAGGAGATCAAGATCATCCTGGCTAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAAATCAGCCAGGTGTGGCAGTGT" +
            "GCACCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATTGCTTGAATCTGGTAGGCAGAGGTTGCAGTGAGCTGAGATCACACCACTGCACTCCTGCCTGGGCAACAGAGCAAGACTCCGTCTCAAACACACACACACAGACACACACACACACACAC" +
            "ACACACACACACACACACACGCAGAGAAACCACAATGAGATACCACCTCATACCAGTCAGAATGGCTATTTTTAAAAAGTCAAAAGATAACAGATGCTGACAAAGTTGCAGAGAAAAGGGAATGCTTATACTCCTCTGGTGGGAGTGTAAATTATTTCAACAA" +
            "CTGTAAAAATCAGTGTGACAATTCCTCACAGAATGAAAAACAGAATTATCATTCGACTGGGAAACCTCATAATTGAGTATATACCCAAAGAAATATGAAATATTATAAAGACACATCCACATGCATGTTCACTGCAGCACTATTCACAATAGCAAAGACACGG" +
            "ACAGACTAAATGCCTATCAATGGCAGACTGGATCAAGAAAATATGGTATGGTCAGATGCGGTGGCTCATGCCTGTAATTCCAGCCCTTTGGGAGGCTGAGGCAGGTGGATTGCCTGAGCTTAGAAGTTTGAGACCACTCTGGGCAACATGGCAAAATTTTGTC" +
            "TCCACAGAAGATACAAAAAAAAAAAAAAAAAA";

    @Test
    public void testInversionFilter() throws Exception {

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("100M"), true, new SimpleInterval("1", 10000, 10100), 60, 1, 100, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("100M"), false, new SimpleInterval("1", 20100, 20200), 60, 101, 200, 0);
        final ChimericAlignment breakpoint1 = new ChimericAlignment(region1, region2, "", "", new ArrayList<>());

        final Tuple2<BreakpointAllele, Tuple2<Tuple2<String, String>, ChimericAlignment>> breakpointTuple1 = new Tuple2<>(new BreakpointAllele.BreakpointAlleleInversion(breakpoint1), new Tuple2<>(new Tuple2<>("1", "contig-1"), breakpoint1));

        Assert.assertTrue(SVVariantCallerUtils.isInversion(breakpointTuple1._1()));

        final AlignmentRegion region3 = new AlignmentRegion("4", "contig-7", TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region4 = new AlignmentRegion("4", "contig-7", TextCigarCodec.decode("137S141M"), false, new SimpleInterval("10", 38342908, 38343049), 60, 138, 278, 0);
        final ChimericAlignment breakpoint2 = new ChimericAlignment(region3, region4, "", "", new ArrayList<>());

        final Tuple2<BreakpointAllele, Tuple2<Tuple2<String, String>, ChimericAlignment>> breakpointTuple2 = new Tuple2<>(new BreakpointAllele.BreakpointAlleleInversion(breakpoint2), new Tuple2<>(new Tuple2<>("14399","contig-7"), breakpoint2));

        Assert.assertFalse(SVVariantCallerUtils.isInversion(breakpointTuple2._1()));

        final AlignmentRegion region5 = new AlignmentRegion("3", "contig-7", TextCigarCodec.decode("137M141S"), true, new SimpleInterval("19", 38343346, 38343483), 60, 1, 137, 0);
        final AlignmentRegion region6 = new AlignmentRegion("3", "contig-7", TextCigarCodec.decode("137S141M"), false, new SimpleInterval("19", 38342908, 38343049), 60, 138, 278, 0);
        final ChimericAlignment breakpoint3 = new ChimericAlignment(region5, region6, "", "", new ArrayList<>());

        final Tuple2<BreakpointAllele, Tuple2<Tuple2<String, String>, ChimericAlignment>> breakpointTuple3 = new Tuple2<>(new BreakpointAllele.BreakpointAlleleInversion(breakpoint3), new Tuple2<>(new Tuple2<>("14399","contig-7"), breakpoint3));

        Assert.assertTrue(SVVariantCallerUtils.isInversion(breakpointTuple3._1()));
    }

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegions() throws Exception {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("532M87S"), true, new SimpleInterval("8", 118873207, 118873739), 60, 1, 532, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("518S29M72S"), false, new SimpleInterval("1", 175705642, 175705671), 3, 519, 547, 0);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("543S76M"), false, new SimpleInterval("1", 118875262, 118875338), 60, 544, 619, 0);
        final ArrayList<AlignmentRegion> alignmentRegionList = new ArrayList<>();
        alignmentRegionList.add(region1);
        alignmentRegionList.add(region2);
        alignmentRegionList.add(region3);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = SVVariantCallerInternal.getBreakpointAlignmentsFromAlignmentRegions(contigSequence, alignmentRegionList, 50);
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "contig-1");
        Assert.assertEquals(chimericAlignment.region1, region1);
        Assert.assertEquals(chimericAlignment.region2, region3);
        Assert.assertEquals(chimericAlignment.homology, "");
        Assert.assertEquals(chimericAlignment.insertedSequence, "GAGATAGAGTC");
    }

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegionsWithOverlappingAlignmentRegion() throws Exception {
        final byte[] contigSequence = "ACTAGAGCATCTACGTGTTCCTGTGGTTTTGGAGCAAGAGTGATTTGAGTTTCAGAGATTTTTACTAATTCTTCTTCCCCTACCAGAAAAAAAGATCTTACCATTTGAGAGTGAGATGTAAACCCAGCCCTGTCTGACCTGAGTCTGTGCCCTAAGCCTATGCTAAGCCAAGCAGTGCCTGGAGCCACCACAGGTCCACACAATTCGTTAACATGATGAAGCAAGGATGGAAATTGGACAAAATAGTGTGCCTACTGAATCTAAGAATGAAAAATGATTGCACTCCTACTCTGAGTGCTTTGGAGCACTGCCCAGTTGGGCAAAGGGTCAGCGCCTGGGCAGAGGTCCCCACAACCTGGCAGGAGTGTGGTCGGCCACCCTATGGGCCTCCATCATGTGCAGTGACAGCGGGGCTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGCCTATAATACAGTTCCGCATGGCCACGTCATACAACCCTGTCATATTGGTGAGCAATTGCTGTGTAGCCAAAGACCCCAAAACTCAAACAGCATTTATTATTATTGCCCCCATGTCTGAGAGTCAGATGTGCATTTGCTGATCTCAGCTTGTTTGAGCTGCTGCAGGGTTGGGGCTCTGCTCCAGGCAGGCTTAGCTGTCACCACATGCACACATACATTCTGGGCCTCTGCTGCGCGCGTCACGTTCACTGAAGATCTTGGGATTGGGAGTTAGGGCGGTGGGAGGGCCCAGCAAAGTCACCTGGCGATGGCAGGGACACAGGGAGGAATGTAGAATGGGGCCGATGATGGGACCCACACGTCTGCAAAGCTGCGGTCTCCTTGAGGGGTGGAGACAGCAACAACTCACCGCACGCGGTGCTTCAGTTCACCATCTCCCTGGGACATTAGGGGGCCCCGTGTTATCTCATTTTGCTCTGGTTTGCATTAGTTTTTTATCACTTCGTAGATGAAGCCACTGACACCCAGAGAGGGAAAGTGGCCTGACCAAGGGCCACAGCAGGGGAGCGAAGGAGCCCCACAGTTCGGCAGGAACACAGCCTCTCCCTGGCTTTCAGGTTCACTGACATCTTCTCATGGCCTCTGTAACTCACCAGGCATCAGGGTGTAGTCCTTAGACCAGTGTCCCACAGCTGCCACAGAGTGGGAGCTCACCATCAGTTATAAGTCACTAGAAAGGCTTTTGGACATTATAAGCTACAATGGAAAATAAGTCATCTGTGGATTTTTGTGACAGATTCCAAAAATTTGAATATTTTGTCTACTTAGGTTTTTGGTTAATTTTATCCTCAAAACTGTTCTGCAGTGATTAAGCTGTACAAACTGCATCATGGGCGAATTGGCATATTCAGAAATGACTGATATTCTTGATTTCAGTTTTTTACTTTGTATGTAGCTCCTCAAGGAAAC".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("487M1006S"), true, new SimpleInterval("20", 23102817, 23103304), 60, 1, 487, 1);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("483S42M968S"), false, new SimpleInterval("20", 23103196, 23103238), 60, 484, 525, 2);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("523S970M"), true, new SimpleInterval("20", 23103633, 23104603), 60, 524, 1493, 3);
        final ArrayList<AlignmentRegion> alignmentRegionList = new ArrayList<>();
        alignmentRegionList.add(region1);
        alignmentRegionList.add(region2);
        alignmentRegionList.add(region3);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = SVVariantCallerInternal.getBreakpointAlignmentsFromAlignmentRegions(contigSequence, alignmentRegionList, 50);
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "contig-1");
        Assert.assertEquals(chimericAlignment.region1, region1);
        Assert.assertEquals(chimericAlignment.region2, region3);
        Assert.assertEquals(chimericAlignment.homology, "");
        Assert.assertEquals(chimericAlignment.insertedSequence, "TGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCA");
        Assert.assertEquals(chimericAlignment.insertionMappings.size(), 1);
        Assert.assertEquals(chimericAlignment.insertionMappings.get(0), "1-contig-1:484-525:20,23103196,-,483S42M968S,60,2");
    }

    @Test
    public void testGetAssembledBreakpointFromAlignmentRegionsStrangeLeftBreakpoint() throws Exception {
        final byte[] sequence = LONG_CONTIG1.getBytes();
        AlignmentRegion region1 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, new SimpleInterval("chr19", 20138007, 20142231), 60, 1, sequence.length - 1986, 36);
        AlignmentRegion region2 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, new SimpleInterval("chr19", 20152030, 20154634), 60, 3604, sequence.length, 36);
        final ArrayList<AlignmentRegion> alignmentRegionList = new ArrayList<>();
        alignmentRegionList.add(region1);
        alignmentRegionList.add(region2);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = SVVariantCallerInternal.getBreakpointAlignmentsFromAlignmentRegions(sequence, alignmentRegionList, 50);
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "702700");
        Assert.assertEquals(chimericAlignment.region1, region1);
        Assert.assertEquals(chimericAlignment.region2, region2);

        final Tuple2<SimpleInterval, SimpleInterval> leftAndRightBreakpointsOnReferenceLeftAlignedForHomology = chimericAlignment.getLeftAndRightBreakpointsOnReferenceLeftAlignedForHomology();

        Assert.assertEquals(leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._1(), new SimpleInterval("chr19", 20138007, 20138007));
        Assert.assertEquals(leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._2(), new SimpleInterval("chr19", 20152651, 20152651));
    }

    @Test
    public void testTreatAlignmentRegionAsInsertion() throws Exception {
        AlignmentRegion overlappingRegion1 = new AlignmentRegion("overlap", "22", TextCigarCodec.decode("47S154M"), false, new SimpleInterval("19", 48699881, 48700035), 60, 1, 154, 0);
        AlignmentRegion overlappingRegion2 = new AlignmentRegion("overlap", "22", TextCigarCodec.decode("116H85M"), true, new SimpleInterval("19", 48700584, 48700669), 60, 117, 201, 0);

        Assert.assertTrue(SVVariantCallerInternal.treatNextAlignmentRegionInPairAsInsertion(overlappingRegion1, overlappingRegion2, 50));
    }

}