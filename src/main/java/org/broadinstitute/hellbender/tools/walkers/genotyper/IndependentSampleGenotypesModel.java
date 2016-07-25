package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.ArrayList;
import java.util.List;

/**
 * This class delegates genotyping to allele count- and ploidy-dependent {@link GenotypeLikelihoodCalculator}s
 * under the assumption that sample genotypes are independent conditional on their population frequencies.
 */
public final class IndependentSampleGenotypesModel {
    private static final int DEFAULT_CACHE_PLOIDY_CAPACITY = 10;
    private static final int DEFAULT_CACHE_ALLELE_CAPACITY = 50;

    private final int cacheAlleleCountCapacity;
    private final int cachePloidyCapacity;

    private GenotypeLikelihoodCalculator[][] likelihoodCalculators;
    private final GenotypeLikelihoodCalculators calculators;

    public IndependentSampleGenotypesModel() { this(DEFAULT_CACHE_PLOIDY_CAPACITY, DEFAULT_CACHE_ALLELE_CAPACITY); }

    /**
     *  Initialize model with given maximum allele count and ploidy for caching
     */
    public IndependentSampleGenotypesModel(final int calculatorCachePloidyCapacity, final int calculatorCacheAlleleCapacity) {
        Utils.validateArg(calculatorCachePloidyCapacity >= 0, () -> "the ploidy provided cannot be negative: " + calculatorCachePloidyCapacity);
        Utils.validateArg(calculatorCacheAlleleCapacity >= 0, () -> "the maximum allele index provided cannot be negative: " + calculatorCacheAlleleCapacity);

        cachePloidyCapacity = calculatorCachePloidyCapacity;
        cacheAlleleCountCapacity = calculatorCacheAlleleCapacity;
        likelihoodCalculators = new GenotypeLikelihoodCalculator[calculatorCachePloidyCapacity][calculatorCacheAlleleCapacity];
        calculators = new GenotypeLikelihoodCalculators();
    }

    /**
     *
     * @param genotypingAlleles
     * @param data
     * @param <A>
     * @return
     */
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                            final GenotypingData<A> data) {

        Utils.nonNull(genotypingAlleles, "the allele cannot be null");
        Utils.nonNull(data, "the genotyping data cannot be null");

        // prepare data, get information necessary
        final AlleleListPermutation<A> permutation = data.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = AlleleLikelihoodMatrixMapper.newInstance(permutation);

        final int sampleCount = data.numberOfSamples();
        final PloidyModel ploidyModel = data.ploidyModel();
        final int alleleCount = genotypingAlleles.numberOfAlleles();

        // TODO: why not the following early return?
//        if(sampleCount==0) { // ever possible to be negative?
//            return null; // OK not null, but something similar
//        }

        // result container
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        // walk over all samples, genotype
        GenotypeLikelihoodCalculator likelihoodsCalculator = sampleCount > 0 ? getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount) : null;
        for (int i = 0; i < sampleCount; i++) {

            // get a new likelihoodsCalculator if this sample's ploidy differs from the previous sample's
            final int samplePloidy = ploidyModel.samplePloidy(i);
            if (samplePloidy != likelihoodsCalculator.ploidy()) {
                likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy, alleleCount);
            }

            final LikelihoodMatrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.apply(data.readLikelihoods().sampleMatrix(i));
            genotypeLikelihoods.add(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods);
    }

    /**
     *
     * @param samplePloidy
     * @param alleleCount
     * @return
     */
    private GenotypeLikelihoodCalculator getLikelihoodsCalculator(final int samplePloidy, final int alleleCount) {
        if (samplePloidy >= cachePloidyCapacity || alleleCount >= cacheAlleleCountCapacity) {
            return calculators.getInstance(samplePloidy, alleleCount);
        }
        final GenotypeLikelihoodCalculator result = likelihoodCalculators[samplePloidy][alleleCount];
        if (result != null) {
            return result;
        } else {
            final GenotypeLikelihoodCalculator newOne = calculators.getInstance(samplePloidy, alleleCount);
            likelihoodCalculators[samplePloidy][alleleCount] = newOne;
            return newOne;
        }
    }
}