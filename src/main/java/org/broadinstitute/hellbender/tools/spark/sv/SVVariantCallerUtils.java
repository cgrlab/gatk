package org.broadinstitute.hellbender.tools.spark.sv;

/**
 * Variants utility functions helping calling structural variants.
 */
final class SVVariantCallerUtils {

    /**
     * A naive way to test if a chimeric alignment supports an inversion event.
     * Returning true does not necessarily mean there's actually an inversion.
     */
    static boolean supportsInversionNaive(final ChimericAlignment chimericAlignment) {
        return chimericAlignment.region1.forwardStrand != chimericAlignment.region2.forwardStrand;
    }

    /**
     * Test if a {@link BreakpointAllele} is an {@link org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion}
     */
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

}
