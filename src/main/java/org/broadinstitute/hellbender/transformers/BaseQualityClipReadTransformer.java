package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.tools.ClipReads;
import org.broadinstitute.hellbender.utils.clipping.ClippingOp;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import static org.broadinstitute.hellbender.utils.read.ReadUtils.*;
import static org.broadinstitute.hellbender.utils.read.ReadUtils.setDeletionBaseQualities;

/**
 * Created by markw on 10/21/16.
 */

public class BaseQualityClipReadTransformer implements ReadTransformer {
    public static final long serialVersionUID = 1L;


    /**
     * The quality score based read clipper will be applied to the reads using this
     * quality score threshold.
     */
    @Argument(fullName = "qTrimmingThreshold", shortName = "QT", doc = "Quality score trimmer threshold", optional = true)
    private int qTrimmingThreshold = 15;

    /**
     * Clip bases from the read from
     * <p/>
     * argmax_x{ \sum{i = x + 1}^l (qTrimmingThreshold - qual)
     * <p/>
     * to the end of the read.  This is blatantly stolen from BWA.
     * <p/>
     * Walk through the read from the end (in machine cycle order) to the beginning, calculating the
     * running sum of qTrimmingThreshold - qual.  While we do this, we track the maximum value of this
     * sum where the delta > 0.  After the loop, clipPoint is either -1 (don't do anything) or the
     * clipping index in the read (from the end).
     *
     */
    @Override
    public GATKRead apply(GATKRead read) {
        GATKRead readClippedRightEnd = clipReadRightEnd(read);
        return clipReadLeftEnd(readClippedRightEnd);
    }

    private GATKRead clipReadRightEnd(GATKRead read) {
        int readLength = read.getLength();
        byte[] quals = read.getBaseQualities();

        int clipSum = 0, lastMax = -1, clipPoint = -1; // -1 means no clip
        for (int i = readLength - 1; i >= 0; i--) {
            clipSum += (qTrimmingThreshold - quals[i]);
            if (clipSum >= 0 && (clipSum >= lastMax)) {
                lastMax = clipSum;
                clipPoint = i;
            }
        }

        if (clipPoint != -1) {
            int newLength = clipPoint;
            final byte[] newBases = new byte[newLength];
            final byte[] newQuals = new byte[newLength];
            final int copyStart = 0;
            System.arraycopy(read.getBases(), copyStart, newBases, 0, newLength);
            System.arraycopy(read.getBaseQualities(), copyStart, newQuals, 0, newLength);
            final GATKRead clippedRead = read.copy();
            clippedRead.setBaseQualities(newQuals);
            clippedRead.setBases(newBases);
            return clippedRead;
        } else {
            return read;
        }
    }

    private GATKRead clipReadLeftEnd(GATKRead read) {
        int readLength = read.getLength();
        byte[] quals = read.getBaseQualities();

        int clipSum = 0, lastMax = -1, clipPoint = -1; // -1 means no clip
        for (int i = 0; i < readLength; i++) {
            clipSum += (qTrimmingThreshold - quals[i]);
            if (clipSum >= 0 && (clipSum >= lastMax)) {
                lastMax = clipSum;
                clipPoint = i;
            }
        }

        if (clipPoint != -1) {
            int newLength = readLength - 1 - clipPoint;
            final byte[] newBases = new byte[newLength];
            final byte[] newQuals = new byte[newLength];
            final int copyStart = clipPoint + 1;
            System.arraycopy(read.getBases(), copyStart, newBases, 0, newLength);
            System.arraycopy(read.getBaseQualities(), copyStart, newQuals, 0, newLength);
            final GATKRead clippedRead = read.copy();
            clippedRead.setBaseQualities(newQuals);
            clippedRead.setBases(newBases);
            return clippedRead;
        } else {
            return read;
        }
    }

}


