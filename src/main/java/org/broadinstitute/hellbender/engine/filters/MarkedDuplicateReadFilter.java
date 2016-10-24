package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;

public final class MarkedDuplicateReadFilter extends ReadFilter implements Serializable {

    private static final long serialVersionUID = 1L;

    //Filters out optical duplicate reads (marked with the OD tag)
    @Override
    public boolean test( final GATKRead read ) {
        if (read.hasAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)
                && read.getAttributeAsInteger(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME) != 0) {
            return false;
        }
        else if (read.isDuplicate()) {
            return false;
        }
        else {
            return true;
        }
    }
}
