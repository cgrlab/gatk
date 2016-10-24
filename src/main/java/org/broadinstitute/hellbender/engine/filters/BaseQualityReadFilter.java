package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import java.io.Serializable;

public final class BaseQualityReadFilter extends ReadFilter implements Serializable {

    private static final long serialVersionUID = 1L;

    @Argument(fullName="baseQualThresh", shortName="baseQualThresh", optional=true)
    public int BASE_QUALITY_THRESHOLD = 15;

    @Argument(fullName="baseQualFrac", shortName="baseQualFrac", optional=true)
    public float BASE_QUALITY_MAX_FRAC = 0.05f;

    public BaseQualityReadFilter( ) { }

    public BaseQualityReadFilter(final int quality_threshold, final float quality_frac) {
        this.BASE_QUALITY_THRESHOLD = quality_threshold;
        this.BASE_QUALITY_MAX_FRAC = quality_frac;
    }

    //Filters out reads with at least BASE_QUALITY_MAX_COUNT fraction of base qualities that are below BASE_QUALITY_THRESHOLD
    @Override
    public boolean test( final GATKRead read ) {
        short numBelowThreshold = 0;
        if (read.getBaseQualityCount() > 0) {
            for (int q : read.getBaseQualities()) {
                if (q < BASE_QUALITY_THRESHOLD)
                    numBelowThreshold++;
            }
        }
        return numBelowThreshold < Math.ceil(BASE_QUALITY_MAX_FRAC*read.getLength());
    }
}
