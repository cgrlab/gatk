package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.transformers.DUSTReadTransformer;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;

import static org.testng.Assert.*;

/**
 * Created by markw on 10/24/16.
 */
public class MarkedDuplicateReadFilterTest {

    private SAMFileHeader header;

    @BeforeMethod
    public void before() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
    }

    @Test
    public void testTest() throws Exception {

        MarkedDuplicateReadFilter filter = new MarkedDuplicateReadFilter();

        String seq_in = "ATCGATCG";
        String qual_in = "IIIIIIII";
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 10, seq_in.getBytes(),qual_in.getBytes());
        Assert.assertEquals(filter.test(read_in),true);

        read_in.setAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME,1);
        Assert.assertEquals(filter.test(read_in),false);

        read_in.setAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME,0);
        Assert.assertEquals(filter.test(read_in),true);

        read_in.clearAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME);
        read_in.setIsDuplicate(true);
        Assert.assertEquals(filter.test(read_in),false);

        read_in.setIsDuplicate(false);
        Assert.assertEquals(filter.test(read_in),true);
    }

}