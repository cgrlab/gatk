package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CRAMReferenceIntegrationTest extends CommandLineProgramTest{

    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools");

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }


    @Test(dataProvider="testingDataNoRef", expectedExceptions = UserException.class)
    public void testNoRef(String fileIn, String extOut) throws Exception {
        final File outFile = File.createTempFile(fileIn + ".", extOut);
        outFile.deleteOnExit();
        File ORIG_BAM = new File(TEST_DATA_DIR, fileIn);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @Test(dataProvider="testingDataNoRefMultipleInputs", expectedExceptions = UserException.class)
    public void testNoRefMulti(String fileIn1, String fileIn2, String extOut) throws Exception {
        final File outFile = File.createTempFile(fileIn1 + ".", extOut);
        outFile.deleteOnExit();
        File ORIG_BAM1 = new File(TEST_DATA_DIR, fileIn1);
        File ORIG_BAM2 = new File(TEST_DATA_DIR, fileIn2);
        final String[] args = new String[]{
                "--input" , ORIG_BAM1.getAbsolutePath(),
                "--input" , ORIG_BAM2.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @Test(dataProvider="testingDataWrongRef", expectedExceptions = UserException.class)
    public void testWrongRef(String fileIn, String extOut, String referenceFile) throws Exception {
        final File outFile = File.createTempFile(fileIn + ".", extOut);
        outFile.deleteOnExit();
        File ORIG_BAM = new File(TEST_DATA_DIR, fileIn);
        File reference = new File(TEST_DATA_DIR, referenceFile);
        final String[] args = new String[]{
                "--input" , ORIG_BAM.getAbsolutePath(),
                "--output", outFile.getAbsolutePath(),
                "-R", reference.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @DataProvider(name="testingDataNoRef")
    public Object[][] testingDataNoRef() {
        return new String[][]{
                {"cramtest.cram", ".sam"}
        };
    }

    @DataProvider(name="testingDataNoRefMultipleInputs")
    public Object[][] testingDataNoRefMulti() {
        return new String[][]{
                {"cramtest.sam", "cramtest.cram", ".bam"}
        };
    }

    @DataProvider(name="testingDataWrongRef")
    public Object[][] testingDataWrongRef() {
        return new String[][]{
                {"cramtest.cram", ".sam", "cramtest.fasta"},
        };
    }

}