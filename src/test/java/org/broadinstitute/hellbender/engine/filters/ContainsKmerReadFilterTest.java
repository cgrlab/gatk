package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerSmall;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

public class ContainsKmerReadFilterTest extends BaseTest {

    private SAMFileHeader header;
    private final int kSize = 11;
    private HopscotchSet<SVKmerSmall> kmerSet;

    @BeforeMethod
    public void before() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        String kmerRef = "ATCGAGCGCTAGCGATGGCGCGCGATCGCGCTAGCGCGCTAGC";
        SVKmerizer<SVKmerSmall> kmerizer = new SVKmerizer<>(kmerRef.getBytes(),kSize,new SVKmerSmall(kSize));
        ArrayList<SVKmerSmall> kmerList = new ArrayList<>();
        while (kmerizer.hasNext()) {
            kmerList.add(kmerizer.next());
        }
        kmerSet = new HopscotchSet<>(kmerList);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"GCGCGCGATCG", Boolean.TRUE},
                {"TAGCGATGGCGCGCGATCACGCTAG", Boolean.TRUE},
                {"TAGCGATGGCACGCGATCGAGCTAG", Boolean.FALSE},
                {"TAGCGA", Boolean.FALSE},
                {"", Boolean.FALSE}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testTest(String bases_in, Boolean test_out) throws Exception {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ContainsKmerReadFilter filter = new ContainsKmerReadFilter(ctx.broadcast(kmerSet),kSize);
        byte[] quals = bases_in.getBytes().clone();
        Arrays.fill(quals,(byte)'I');
        SAMUtils.fastqToPhred(quals);
        GATKRead read_in = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 10, bases_in.getBytes(), quals);
        boolean test_i = filter.test(read_in);
        Assert.assertEquals(test_out.booleanValue(),test_i);
    }

}