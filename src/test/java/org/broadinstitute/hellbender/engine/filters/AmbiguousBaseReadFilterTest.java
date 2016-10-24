package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.transformers.DUSTReadTransformer;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

import static org.testng.Assert.*;

public class AmbiguousBaseReadFilterTest {

    private SAMFileHeader header;

    @BeforeMethod
    public void before() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
            {"CGTTTTTCAGTGATTTCTTCATTTTTCAATTCGTCAAGTGGATGTTTCTCATTTTCCATGATTTTCAGTTTTCTTGCCATATTCCACGTCCTACAGTGGA", Boolean.TRUE},
            {"NNNNNTTCAGTGATTTCTTCATTTTTCAATTCGTCAAGTGGATGTTTCTCATTTTCCATGATTTTCAGTTTTCTTGCCATATTCCACGTCCTACAGTGGA", Boolean.FALSE},
            {"CGTTTTTCAGTGATTTCTTCATTTTTCAATTCGTCAAGTGGATGTTTCTCATTTTCCATGATTTTNNNNNTTCTTGCCATATTCCACGTCCTACAGTGGA", Boolean.FALSE},
            {"CGTTTTTCAGTGATTTCTTCATTNNNNAATTCGTCAAGTGGATGTTTCTCATTTTCCATGATTTTCAGTTTTCTTGCCATATTCCACGTCCTACAGTGGA", Boolean.TRUE}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testTest(String seq_in, Boolean test_out) throws Exception {
        AmbiguousBaseReadFilter filter = new AmbiguousBaseReadFilter(0.05f);
        byte[] quals = seq_in.getBytes().clone();
        Arrays.fill(quals,(byte)'I');
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 10, seq_in.getBytes(),quals);
        boolean test_i = filter.test(read_in);
        Assert.assertEquals(test_out.booleanValue(),test_i);
    }

}