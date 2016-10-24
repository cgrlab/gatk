package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import java.util.Arrays;

public class BaseQualityReadFilterTest {

    private SAMFileHeader header;

    @BeforeMethod
    public void before() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"/IIIIIIIIIIIIIIIIIII", Boolean.FALSE},
                {"IIIIIIIIIIIIIIIIIIII", Boolean.TRUE},
                {"00000000000000000000", Boolean.TRUE},
                {"0IIIIIIIIIIIIIIIIII", Boolean.TRUE}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testTest(String quals_in, Boolean test_out) throws Exception {
        BaseQualityReadFilter filter = new BaseQualityReadFilter(15,0.05f);
        byte[] bases = quals_in.getBytes().clone();
        Arrays.fill(bases,(byte)'A');
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 10, bases, SAMUtils.fastqToPhred(quals_in));
        boolean test_i = filter.test(read_in);
        Assert.assertEquals(test_out.booleanValue(),test_i);
    }

}