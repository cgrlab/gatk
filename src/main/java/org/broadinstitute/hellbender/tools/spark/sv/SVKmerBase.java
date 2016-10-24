package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

//Superclass for SVKmer and SVKmerSmall

public abstract class SVKmerBase {

    public enum Base {
        A(0L),
        C(1L),
        G(2L),
        T(3L);

        public final long value;

        Base(final long value) { this.value = value; }
    }

    // Lookup table for reverse-complementing each possible byte value.
    // Each pair of bits represents a base, so you have to reverse bits pairwise and then invert all bits.
    // This is most quickly and easily done with a lookup table.
    protected static final long[] BYTEWISE_REVERSE_COMPLEMENT;
    static {
        BYTEWISE_REVERSE_COMPLEMENT = new long[256];
        for ( int idx = 0; idx != 256; ++idx ) {
            BYTEWISE_REVERSE_COMPLEMENT[idx] = reverseComplementByteValueAsLong(idx);
        }
    }

    protected abstract void serialize(final Kryo kryo, final Output output );
    public abstract SVKmerBase successor(final Base base, final int kSize );
    public abstract SVKmerBase predecessor(final Base base, final int kSize );
    public abstract SVKmerBase reverseComplement( final int kSize );
    public abstract SVKmerBase canonical( final int kSize );
    public abstract Base firstBase( final int kSize );
    public abstract Base lastBase();
    public abstract String toString( final int kSize );

    protected static long reverseComplementByteValueAsLong( final int bIn ) {
        // this turns the 8 bits [b1 b2 b3 b4 b5 b6 b7 b8] into [~b7 ~b8 ~b5 ~b6 ~b3 ~b4 ~b1 ~b2]
        return ~(((bIn & 3) << 6) | (((bIn >> 2) & 3) << 4) | (((bIn >> 4) & 3) << 2) | ((bIn >> 6) & 3)) & 0xffL;
    }

    protected static int fnvLong( final int start, final long toHash ) {
        return fnvInt(fnvInt(start, (int)(toHash >> 32)), (int)toHash);
    }

    protected static int fnvInt( int start, final int toHash ) {
        final int mult = 16777619;
        start ^= (toHash >> 24) & 0xff;
        start *= mult;
        start ^= (toHash >> 16) & 0xff;
        start *= mult;
        start ^= (toHash >> 8) & 0xff;
        start *= mult;
        start ^= toHash & 0xff;
        start *= mult;
        return start;
    }

}
