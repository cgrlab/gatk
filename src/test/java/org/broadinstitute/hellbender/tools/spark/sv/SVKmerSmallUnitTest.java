package org.broadinstitute.hellbender.tools.spark.sv;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer.toKmer;

/**
 * Unit tests for SVKmerSmall and SVKmerizer<SVKmerSmall>.
 */
public class SVKmerSmallUnitTest {
    @Test
    public void testDefaultConstruction() {
        Assert.assertEquals(new SVKmerSmall(10).toString(10), "AAAAAAAAAA");
        Assert.assertEquals(new SVKmerSmall(11).toString(11), "AAAAAAAAAAA");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDefaultConstructionWithTooLargeK() {
        final SVKmerSmall tooBigK = new SVKmerSmall(32);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDefaultConstructionWithTooSmallK() {
        final SVKmerSmall tooSmallK = new SVKmerSmall(0);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"ACGTACGTACGT"}, {"ACGTACGTACGTC"}, {"ACGTACGTACGTCC"}, {"ACGTACGTACGTCCC"}, {"ACGTACGTACGTCCCC"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testConstructionAndToString( final String str ) {
        Assert.assertEquals(str, SVKmerizer.toKmer(str,new SVKmerSmall(str.length())).toString(str.length()));
        Assert.assertEquals(str, SVKmerizer.toKmer(str.getBytes(),new SVKmerSmall(str.length())).toString(str.length()));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testSuccessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.delete(0, 1);
        sb.append('A');
        final SVKmerSmall kkk = SVKmerizer.toKmer(str,new SVKmerSmall(str.length()));
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.A, K).toString(K));
        sb.setCharAt(K-1, 'C');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.C, K).toString(K));
        sb.setCharAt(K-1, 'G');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.G, K).toString(K));
        sb.setCharAt(K-1, 'T');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.T, K).toString(K));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testPredecessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.insert(0, 'A');
        sb.setLength(K);
        final SVKmerSmall kkk = SVKmerizer.toKmer(str,new SVKmerSmall(str.length()));
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.A, K).toString(K));
        sb.setCharAt(0, 'C');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.C, K).toString(K));
        sb.setCharAt(0, 'G');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.G, K).toString(K));
        sb.setCharAt(0, 'T');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.T, K).toString(K));
    }

    @Test
    public void testReverseComplementation() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGA",new SVKmerSmall(8)).reverseComplement(8), SVKmerizer.toKmer("TCGTACGT",new SVKmerSmall(8)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC",new SVKmerSmall(9)).reverseComplement(9), SVKmerizer.toKmer("GACGTACGT",new SVKmerSmall(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTCCGTC",new SVKmerSmall(9)).reverseComplement(9), SVKmerizer.toKmer("GACGGACGT",new SVKmerSmall(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTGCGTC",new SVKmerSmall(9)).reverseComplement(9), SVKmerizer.toKmer("GACGCACGT",new SVKmerSmall(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTTCGTC",new SVKmerSmall(9)).reverseComplement(9), SVKmerizer.toKmer("GACGAACGT",new SVKmerSmall(9)));
    }

    @Test
    public void testCanonicalization() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC",new SVKmerSmall(9)).canonical(9), SVKmerizer.toKmer("ACGTACGTC",new SVKmerSmall(9)));
        Assert.assertEquals(SVKmerizer.toKmer("GACGTACGT",new SVKmerSmall(9)).canonical(9), SVKmerizer.toKmer("ACGTACGTC",new SVKmerSmall(9)));
    }

    @Test
    public void testComparison() {
        final SVKmerSmall kkk1 = SVKmerizer.toKmer("ACGTA",new SVKmerSmall(5));
        final SVKmerSmall kkk2 = SVKmerizer.toKmer("ACGTC",new SVKmerSmall(5));
        Assert.assertTrue(kkk1.compareTo(kkk1) == 0);
        Assert.assertTrue(kkk2.compareTo(kkk2) == 0);
        Assert.assertTrue(kkk1.compareTo(kkk2) < 0);
        Assert.assertTrue(kkk2.compareTo(kkk1) > 0);
    }

    @Test
    public void testHashCode() {
        Assert.assertNotEquals(SVKmerizer.toKmer("TAGCGTA",new SVKmerSmall(7)).hashCode(), SVKmerizer.toKmer("TAGCGTC",new SVKmerSmall(7)).hashCode());
        Assert.assertEquals(SVKmerizer.toKmer("TAGGGTC",new SVKmerSmall(7)).hashCode(), SVKmerizer.toKmer("TAGGGTC",new SVKmerSmall(7)).hashCode());
    }

    @Test
    public void testKmerization() {
        final SVKmerizer<SVKmerSmall> kmerizer = new SVKmerizer<SVKmerSmall>("AAAAATT", 5, new SVKmerSmall(7));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA",new SVKmerSmall(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAT",new SVKmerSmall(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAATT",new SVKmerSmall(5)));
        Assert.assertTrue(!kmerizer.hasNext());
    }

    @Test
    public void testKmerizationAcrossN() {
        final SVKmerizer<SVKmerSmall> kmerizer = new SVKmerizer<SVKmerSmall>("AAAAANTTTTT", 5,new SVKmerSmall(11));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA",new SVKmerSmall(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("TTTTT",new SVKmerSmall(5)));
        Assert.assertTrue(!kmerizer.hasNext());
    }
}