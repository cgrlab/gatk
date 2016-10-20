package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * An argument collection for use with tools that accept one or more input files containing reads
 * (eg., BAM/SAM/CRAM files), and require at least one such input.
 */
public final class RequiredReadInputArgumentCollection extends ReadInputArgumentCollection {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "BAM/SAM/CRAM file containing reads", optional = false)
    public List<String> readFilesNames;

    @Override
    public List<File> getReadFiles() {
        ArrayList<File> ret = new ArrayList<>();
        for (String fn : readFilesNames) {
            ret.add(new File(fn));
        }
        return ret;
    }

    @Override
    public List<Path> getReadPaths() {
        ArrayList<Path> ret = new ArrayList<>();
        for (String fn : readFilesNames) {
            ret.add(BucketUtils.getPath(fn));
        }
        return ret;
    }

    @Override
    public List<String> getReadFilesNames() {
        return new ArrayList<>(readFilesNames);
    }
}
