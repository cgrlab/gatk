package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Summarize genotype statistics from all samples at the site level
 *
 * <p>This annotation collects several genotype-level statistics from all samples and summarizes them in the INFO field. The following statistics are collected:</p>
 * <ul>
 *     <li>Number of called chromosomes (should amount to ploidy * called samples)</li>
 *     <li>Number of no-called samples</li>
 *     <li>p-value from Hardy-Weinberg Equilibrium test</li>
 *     <li>Mean of all GQ values</li>
 *     <li>Standard deviation of all GQ values</li>
 * </ul>
 * <h3>Note</h3>
 * <p>These summaries can all be recomputed from the genotypes on the fly but it is a lot faster to add them here as INFO field annotations.</p>
 */

public final class GenotypeSummaries extends InfoFieldAnnotation {

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        Utils.nonNull(vc);
        if ( ! vc.hasGenotypes() ) {
            return null;
        }

        final Map<String,Object> returnMap = new HashMap<>();
        returnMap.put(GATKVCFConstants.NOCALL_CHROM_KEY, vc.getNoCallCount());

        final DescriptiveStatistics stats = new DescriptiveStatistics();
        for( final Genotype g : vc.getGenotypes() ) {
            if( g.hasGQ() ) {
                stats.addValue(g.getGQ());
            }
        }
        if( stats.getN() > 0L ) {
            returnMap.put(GATKVCFConstants.GQ_MEAN_KEY, String.format("%.2f", stats.getMean()));
            if( stats.getN() > 1L ) {
                returnMap.put(GATKVCFConstants.GQ_STDEV_KEY, String.format("%.2f", stats.getStandardDeviation()));
            }
        }

        return returnMap;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(
                GATKVCFConstants.NOCALL_CHROM_KEY,
                GATKVCFConstants.GQ_MEAN_KEY,
                GATKVCFConstants.GQ_STDEV_KEY);
    }
}
