package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.List;

/**
 * Various utility functions helping calling structural variants.
 */
final class SVVariantCallerUtils {

    /**
     * A naive way to test if a chimeric alignment supports an inversion event.
     * Returning true does not necessarily mean there's actually an inversion.
     */
    @VisibleForTesting
    static boolean involvesStrandSwitch(final ChimericAlignment chimericAlignment) {
        return chimericAlignment.region1.forwardStrand != chimericAlignment.region2.forwardStrand;
    }

    /**
     * Test if a {@link BreakpointAllele} is an {@link org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion}
     */
    @VisibleForTesting
    static boolean isInversion(final BreakpointAllele allele) {
        try {
            final BreakpointAllele.BreakpointAlleleInversion invAllele = (BreakpointAllele.BreakpointAlleleInversion) allele;
            return invAllele.leftAlignedLeftBreakpoint.getContig().equals(invAllele.leftAlignedRightBreakpoint.getContig())
                    &&
                    (invAllele.getInversionType() == BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_3_TO_5 || invAllele.getInversionType() == BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_5_TO_3);
        } catch (final ClassCastException ccex) {
            return false;
        }
    }

    @VisibleForTesting
    static int getTotalHardClipping(final Cigar cigar) {
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        if (cigarElements.size() == 0) {
            return 0;
        }
        if (cigarElements.size() == 1) {
            return cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(0).getLength() : 0;
        }
        return (cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(0).getLength() : 0) +
                (cigarElements.get(cigarElements.size() - 1).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(cigarElements.size() - 1).getLength() : 0);
    }

    /**
     * Returns the number of clipped bases, including both soft and hard, represented in {@code forwardStrandCigar}
     * from the start or from the end
     * @param fromStart             from the start of the template or not
     * @param forwardStrandCigar    the {@link Cigar} in its forward strand representation
     */
    static int getNumClippedBases(final boolean fromStart, final Cigar forwardStrandCigar) {
        final List<CigarElement> elements = forwardStrandCigar.getCigarElements();
        if(elements.size()==1) return 0; // cannot be a giant clip
        final int offset = fromStart ? 1 : -1;

        int result = 0;
        int j = fromStart ? 0 : elements.size() - 1;
        CigarElement ce = forwardStrandCigar.getCigarElement(j);
        while (ce.getOperator().isClipping()) {
            result += ce.getLength();
            j += offset;
            if ( j < 0 || j >= forwardStrandCigar.getCigarElements().size() ) break;
            ce = forwardStrandCigar.getCigarElement(j);
        }
        return result;
    }

}
