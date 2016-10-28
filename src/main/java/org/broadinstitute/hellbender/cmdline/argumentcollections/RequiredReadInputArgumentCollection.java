package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
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
        try {
            ArrayList<Path> ret = new ArrayList<>();
            for (String fn : readFilesNames) {
                ret.add(IOUtils.getPath(fn));
            }
            return ret;
        } catch (IOException io) {
            throw new GATKException(io.getMessage());
        }
    }

    @Override
    public List<String> getReadFilesNames() {
        return new ArrayList<>(readFilesNames);
    }
}
