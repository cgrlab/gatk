package org.broadinstitute.hellbender.engine.filters;

import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerSmall;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.junit.Before;

import java.io.Serializable;

/**
 * Keep reads that contain at least one kmer from a broadcasted HopscotchSet of SVKmerSmall's
 */
public class ContainsKmerReadFilter extends ReadFilter implements Serializable {

    private static final long serialVersionUID = 1L;
    private Broadcast<HopscotchSet<SVKmerSmall>> kmerLibBroadcast;
    private int kSize;

    //Required otherwise compiler complains
    public ContainsKmerReadFilter() {    }

    public ContainsKmerReadFilter(Broadcast<HopscotchSet<SVKmerSmall>> broadcast, int kmer_size) {
        kmerLibBroadcast = broadcast;
        kSize = kmer_size;
    }

    //Filters out reads with at least BASE_QUALITY_MAX_COUNT fraction of base qualities that are below BASE_QUALITY_THRESHOLD
    @Override
    public boolean test( final GATKRead read ) {
        SVKmerizer<SVKmerSmall> kmers = new SVKmerizer<>(read.getBases(),kSize,new SVKmerSmall(kSize));
        while (kmers.hasNext()) {
            if (kmerLibBroadcast.value().contains(kmers.next())) {return true;}
        }
        return false;
    }
}
