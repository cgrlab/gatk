package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Tuple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.tools.ClipReads;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import htsjdk.samtools.SAMUtils;

import java.util.*;

//Masks low-complexity sequences with 'N' and base qualities with '#'

public class DUSTReadTransformer implements ReadTransformer {
    public static final long serialVersionUID = 1L;
    private final Logger logger = LogManager.getLogger(DUSTReadTransformer.class);

    @Argument(fullName="Mask PHRED score", shortName="dustPhred", optional=true)
    public short MASK_PHRED = 15;

    @Argument(fullName="Window size", shortName="dustW", optional=true)
    public int W = 64;

    @Argument(fullName="Score threshold", shortName="dustT", optional=true)
    public float T = 20;

    private static final int NUM_TRIPLETS = 64;

    private class PStruct {
        public int start = 0;
        public int finish = 0;
        public float score = 0;
    }

    @Override
    public GATKRead apply(final GATKRead read) {

        byte[] q = read.getBases();

        //Find convert to 2-bit and assign random bases to N's
        //TODO: reject reads with too many N's
        Random rng = new Random();
        rng.setSeed(0);
        for (int i = 0; i < q.length; i++) {
            if (q[i] == 'N') {
                q[i] = (byte)rng.nextInt(4);
            } else {
                if (q[i] == 'A' || q[i] == 'a') {
                    q[i] = 0x0;
                } else if (q[i] == 'T' || q[i] == 't') {
                    q[i] = 0x1;
                } else if (q[i] == 'C' || q[i] == 'c') {
                    q[i] = 0x2;
                } else if (q[i] == 'G' || q[i] == 'g') {
                    q[i] = 0x3;
                } else {
                    logger.warn("Invalid base " + (char) q[i] + " in read " + read.getName() + ", assigning to random base.");
                    q[i] = (byte)rng.nextInt(4);
                }
            }
        }

        //Start main algorithm
        List<Tuple<Integer,Integer>> res = new ArrayList<>();
        LinkedList<PStruct> P = new LinkedList<>();
        int rv = 0;
        int rw = 0;
        short[] cv = new short[NUM_TRIPLETS];
        short[] cw = new short[NUM_TRIPLETS];
        LinkedList<Integer> w = new LinkedList<>();
        int L = 0;
        int wstart;

        for (int wfinish = 2; wfinish < read.getLength(); wfinish++) {
            wstart = Math.max(wfinish - W + 1, 0);

            //Begin : SAVE_MASKED_REGIONS
            if (!P.isEmpty()) {
                PStruct p1 = P.getLast();
                if (p1.start < wstart) {
                    int l = res.size();
                    if (l > 0) {
                        Tuple<Integer, Integer> rt = res.get(l - 1);
                        if (p1.start < rt.b + 1) {
                            Tuple<Integer, Integer> new_rt = new Tuple<>(rt.a, Math.max(rt.b, p1.finish));
                            res.set(l - 1, new_rt);
                        } else {
                            res.add(new Tuple<Integer, Integer>(p1.start, p1.finish));
                        }
                    } else {
                        res.add(new Tuple<Integer, Integer>(p1.start, p1.finish));
                    }
                    while (p1 != null && p1.start < wstart) {
                        P.removeLast();
                        p1 = P.isEmpty() ? null : P.getLast();
                    }
                }
            }
            //End : SAVE_MASKED_REGIONS

            int t = triplet(q[wfinish-2],q[wfinish-1],q[wfinish]);

            //Begin: SHIFT_WINDOW
            int s;
            if (w.size() >= W - 2) {
                s = (int)w.pop();
                cw[s]--;
                rw -= cw[s];
                if (L > w.size()) {
                    L--;
                    cv[s]--;
                    rv -= cv[s];
                }
            }
            w.add(t);
            L++;
            rw += cw[t];
            cw[t]++;
            rv += cv[t];
            cv[t]++;
            if (cv[t]*10 > 2*T) {
                do {
                    s = (int)w.get(w.size()-L);
                    cv[s]--;
                    rv -= cv[s];
                    L--;
                } while (s != t);
            }
            //End: SHIFT_WINDOW

            if (rw*10 > L*T) {
                //Begin: FIND_PERFECT

                short[] c = cv.clone();
                int r = rv;
                ListIterator<PStruct> iter = P.listIterator();
                PStruct perf = iter.hasNext() ? iter.next() : null;
                float max_score = 0;
                for (int i = w.size() - L - 1; i >= 0; i--) {
                    t = w.get(i);
                    r += c[t];
                    c[t]++;
                    float new_score = ((float)r)/(w.size() - i - 1);
                    if (new_score*10 > T) {
                        while (perf != null && perf.start >= i + wstart) {
                            max_score = Math.max(max_score,perf.score);
                            perf = iter.hasNext() ? iter.next() : null;
                        }
                        if (new_score >= max_score) {
                            max_score = new_score;
                            PStruct new_perf = new PStruct();
                            new_perf.start = i + wstart;
                            new_perf.finish = w.size() + 1 + wstart;
                            new_perf.score = new_score;
                            if (iter.hasPrevious()) {
                                iter.previous();
                            } else {
                                iter = P.listIterator();
                            }
                            iter.add(new_perf);
                            if (iter.hasNext()) {iter.next();}
                        }
                    }
                }

                //End: FIND_PERFECT
            }
        }

        wstart = Math.max(0,read.getLength() - W + 1);

        PStruct p1;
        while (!P.isEmpty()) {

            //Begin: SAVE_MASKED_REGIONS (TODO: copied code from above)
            p1 = P.getLast();
            if (p1.start < wstart) {
                int l = res.size();
                if (l > 0) {
                    Tuple<Integer,Integer> rt = res.get(l - 1);
                    if (p1.start < rt.b + 1) {
                        Tuple<Integer,Integer> new_rt = new Tuple<>(rt.a,Math.max(rt.b,p1.finish));
                        res.set(l - 1,new_rt);
                    } else {
                        res.add(new Tuple<Integer,Integer>(p1.start,p1.finish));
                    }
                } else {
                    res.add(new Tuple<Integer,Integer>(p1.start,p1.finish));
                }
                while (p1 != null && p1.start < wstart) {
                    P.removeLast();
                    p1 = P.isEmpty() ? null : P.getLast();
                }
            }
            //End: SAVE_MASKED_REGIONS

            wstart++;
        }

        //Mask the intervals and base qualities
        byte[] q_new = read.getBases();
        for (Tuple<Integer,Integer> interval : res) {
            for (int i = interval.a; i <= interval.b; i++) {
                q_new[i] = 'N';
            }
        }
        read.setBases(q_new);

        if (read.getBaseQualityCount() == read.getLength()) {
            byte[] qual_new = read.getBaseQualities();
            for (Tuple<Integer,Integer> interval : res) {
                for (int i = interval.a; i <= interval.b; i++) {
                    qual_new[i] = (byte)MASK_PHRED;
                }
            }
            read.setBaseQualities(qual_new);
        }

        return read;
    }

    //Returns integer 0-63 of a given triplet
    private int triplet(byte b1, byte b2, byte b3) {
        return (b1 | (b2 << 2) | (b3 << 4));
    }

}
