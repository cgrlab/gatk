package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

import static org.testng.Assert.*;

/**
 * Created by markw on 10/24/16.
 */
public class BaseQualityClipReadTransformerTest {

    private SAMFileHeader header;

    @BeforeMethod
    public void before() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
            {"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "CTCAAGTACAAGCTGATCCAGACCTACAGGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTGCAATCC",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"},

            {"####################################################################################################",
            "",
            ""},

            {"############################IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "GGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTGCAATCC",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"},

            {"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII######",
            "CTCAAGTACAAGCTGATCCAGACCTACAGGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTG",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testApply(String quals_in, String bases_out, String quals_out) throws Exception {
        String bases_in = "CTCAAGTACAAGCTGATCCAGACCTACAGGGTGATGTCATTAGAGGCACTGATAACACACACACTATGGGGTGGGGGTGGACAGTTCCCCACTGCAATCC";

        BaseQualityClipReadTransformer trans = new BaseQualityClipReadTransformer();
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 10, bases_in.getBytes(),SAMUtils.fastqToPhred(quals_in));
        GATKRead read_out = trans.apply(read_in);
        Assert.assertEquals(read_out.getBases(),bases_out.getBytes());
        Assert.assertEquals(read_out.getBaseQualities(),SAMUtils.fastqToPhred(quals_out));
    }

}