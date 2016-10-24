package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import java.io.File;

public class PathSeqKmerSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return PathSeqKmerSpark.class.getSimpleName();
    }

    @Test(groups = "spark")
    public void test() throws Exception {
        final File expectedFile = getTestFile("kmer.hss");

        final File ref = getTestFile("hg19mini.fasta");
        final File output = createTempFile("test", ".hss");
        if (!output.delete()) {
            Assert.fail();
        }

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument("reference", ref);
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        Assert.assertEquals(FileUtils.readFileToByteArray(output), FileUtils.readFileToByteArray(expectedFile));
    }

}
