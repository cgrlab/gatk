package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public final class AmbiguousBaseReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;

    @Argument(fullName="ambigFilterFrac", shortName="ambigFilterFrac", optional=true)
    public float N_FRAC = 0.05f;

    public AmbiguousBaseReadFilter() { }

    public AmbiguousBaseReadFilter( final float n_frac ) { this.N_FRAC = n_frac; }

    //Filters out reads with more than a threshold number of N's
    @Override
    public boolean test( final GATKRead read ) {
        final int N_max = (int)(read.getLength()*N_FRAC);
        int num_N = 0;
        for (final byte base : read.getBases()) {
            if (!BaseUtils.isRegularBase(base)) {
                num_N++;
                if (num_N > N_max) {return false;}
            }
        }
        return num_N <= N_max;
    }
}
