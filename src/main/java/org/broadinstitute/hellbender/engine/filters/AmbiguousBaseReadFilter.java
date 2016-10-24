package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;

public final class AmbiguousBaseReadFilter extends ReadFilter implements Serializable {

    private static final long serialVersionUID = 1L;

    @Argument(fullName="ambigFilterFrac", shortName="ambigFilterFrac", optional=true)
    public float N_FRAC = 0.05f;

    public AmbiguousBaseReadFilter() { }

    public AmbiguousBaseReadFilter( final float n_frac ) { this.N_FRAC = n_frac; }

    //Filters out reads with more than a threshold number of N's
    @Override
    public boolean test( final GATKRead read ) {
        int num_N = 0;
        for (int i = 0; i < read.getLength(); i++) {
            if (read.getBase(i) == 'N') {
                num_N++;
            }
        }
        return (num_N < Math.ceil(read.getLength()*N_FRAC));
    }
}
