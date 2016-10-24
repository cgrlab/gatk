package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * An immutable SVKmerSmall. This class is the same as SVKmer but uses a 1 long instead of 2 longs to store the kmer.
 * K must be between 1 and 31
 * Canonicalization is unimplemented for even K.
 */
@DefaultSerializer(SVKmerSmall.Serializer.class)
public class SVKmerSmall extends SVKmerBase implements Comparable<SVKmerSmall> {
    // these are treated as K-bit unsigned integers
    private final long valLow; // Kmer bits

    // Lookup table for reverse-complementing each possible byte value.
    // Each pair of bits represents a base, so you have to reverse bits pairwise and then invert all bits.
    // This is most quickly and easily done with a lookup table.
    private static final long[] BYTEWISE_REVERSE_COMPLEMENT;
    static {
        BYTEWISE_REVERSE_COMPLEMENT = new long[256];
        for ( int idx = 0; idx != 256; ++idx ) {
            BYTEWISE_REVERSE_COMPLEMENT[idx] = reverseComplementByteValueAsLong(idx);
        }
    }

    /**
     *  Makes an empty SVKmerSmall.  If you call toString on it, it'll look like poly-A.
     */
    public SVKmerSmall( final int kSize ) {
        Utils.validateArg(kSize >= 1 && kSize < 32, "K must be between 1 and 31.");
        valLow = 0;
    }

    public SVKmerSmall( final SVKmerSmall that ) { this.valLow = that.valLow; }

    private SVKmerSmall( final long valLow ) { this.valLow = valLow; }

    protected SVKmerSmall( final Kryo kryo, final Input input ) {
        valLow = input.readLong();
    }

    protected void serialize( final Kryo kryo, final Output output ) {
        output.writeLong(valLow);
    }

    /**
     * Returns a new SVKmerSmall that's like this one, but with its leading base discarded and a new one added to the end.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.successor(SVKmerSmall.Base.C,5).toString(5) is "CTGAC".
     * @param base must be 0, 1, 2, or 3, corresponding to A, C, G, or T.
     */
    public final SVKmerSmall successor( final Base base, final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        // note we multiply kSize by two in SVKmerSmall because we no longer divide the bits into two longs
        final long mask = (1L << kSize*2) - 1L;
        // move all the bits up two places, OR in the pair of successor bits at the bottom, and mask to kSize bits
        final long newV2 = ((valLow << 2) | (base.value & 3L)) & mask;
        return new SVKmerSmall(newV2);
    }

    /**
     * Returns a new SVKmerSmall that's like this one, but with its trailing base discarded and a new one added to the start.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.predecessor(SVKmerSmall.Base.T,5).toString(5) is "TACTG".
     * @param base must be 0, 1, 2, or 3, corresponding to A, C, G, or T.
     */
    public final SVKmerSmall predecessor( final Base base, final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        // note we multiply kSize by two in SVKmerSmall because we no longer divide the bits into two longs
        final long mask = (1L << kSize*2) - 1L;
        // move all the bits down two places, OR in the bottom two bits from valHigh at the top, and mask to kSize bits
        final long newV2 = ((valLow >> 2) | (base.value << (kSize*2-2))) & mask;
        return new SVKmerSmall(newV2);
    }

    /**
     * Returns a new SVKmerSmall that's the reverse-complement of this one.
     * E.g., if kmer.toString(5) is "ACTGA", then kmer.rc(5).toString(5) is "TCAGT".
     */
    public final SVKmerSmall reverseComplement( final int kSize ) {
        // bit hack to make a long value with the kSize least significant bits set to 1
        // note we multiply kSize by two in SVKmerSmall because we no longer divide the bits into two longs
        final long mask = (1L << kSize*2) - 1L;
        // number of unused bits at the top
        final int compK = 64 - kSize*2;
        // move the significant bits up to the top, reverse complement, and mask to kSize bits.
        return new SVKmerSmall(reverseComplement(valLow << compK) & mask);
    }

    /**
     * Returns a SVKmerSmall that is a canonical representation of this one.
     * An odd-K SVKmerSmall is in canonical form if its middle base is A or C.
     * The reverse-complement of a non-canonical SVKmerSmall is a canonical SVKmerSmall, and vice versa.  (Think about it.)
     * Canonical form is not defined for even-K Kmers (too expensive to compute routinely).
     */
    public SVKmerSmall canonical( final int kSize ) {
        Utils.validateArg( (kSize & 1) != 0, "K must be odd to canonicalize.");
        // for odd-size kmers, the high bit of the middle base is in least significant position in valHigh.
        // test its value by ANDing with 1.  if it's zero the middle base is A or C and we're good to go.
        if ( ((valLow >> kSize) & 1L) == 0 ) return this;
        // middle base is G or T.  reverse complement.
        return reverseComplement(kSize);
    }

    public final Base firstBase( final int kSize ) { return Base.values()[(int)(valLow >> (kSize*2-2))]; }
    public final Base lastBase() { return Base.values()[(int)(valLow & 3)]; }

    @Override
    public boolean equals( final Object obj ) {
        return obj instanceof SVKmerSmall && equals((SVKmerSmall)obj);
    }

    public final boolean equals( final SVKmerSmall that ) {
        return this.valLow == that.valLow;
    }

    @Override
    public final int hashCode() {
        // 32-bit FNV-1a algorithm
        return fnvLong((int)2166136261L, valLow);
    }

    /**
     * SVKmerSmall comparison is consistent with equals.
     * It's also the same as the lexicographic ordering you'd get using toString on the Kmers.
     */
    @Override
    public final int compareTo( final SVKmerSmall that ) {
        return Long.compare(this.valLow, that.valLow);
    }

    /**
     * Not an override.  An SVKmerSmall doesn't know what K is, so it has to be supplied.
     */
    public final String toString( final int kSize ) {
        final StringBuilder sb = new StringBuilder(kSize);

        // we'll produce the string in reverse order and reverse it at the end
        long val = valLow;
        for ( int nChars = kSize/2; nChars > 0; --nChars ) {
            // grab the two least significant bits to index into the BASE_CHARS array
            sb.append(BaseUtils.BASE_CHARS[(int)val & 3]);
            // roll the whole mess down two bits
            val >>= 2;
        }
        // this for loop will have one more iteration than the previous one if kSize is odd
        for ( int nChars = (kSize+1)/2; nChars > 0; --nChars ) {
            // grab two lowest bits to index into array
            sb.append(BaseUtils.BASE_CHARS[(int)val & 3]);
            // move 'em down
            val >>= 2;
        }
        // we built the string in least-significant to most-significant bit order.  reverse it now.
        return sb.reverse().toString();
    }


    // Reverse-complement a long by taking the reverse-complement of each of its bytes in reverse order.
    private static long reverseComplement( long val ) {
        // process val one byte at a time
        long result = BYTEWISE_REVERSE_COMPLEMENT[(int)val & 0xFF]; // handle the low-order byte
        int nBytes = 8;
        while ( --nBytes != 0 ) { // pre-decrementing:  we'll go through the loop 7 times
            // rotate down by a byte
            val >>= 8;
            // rotate up by a byte and OR in the reverse complement of the next byte
            result = (result << 8) | BYTEWISE_REVERSE_COMPLEMENT[(int)val & 0xFF];
        }
        return result;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SVKmerSmall> {
        @Override
        public void write(final Kryo kryo, final Output output, final SVKmerSmall svKmer ) {
            svKmer.serialize(kryo, output);
        }

        @Override
        public SVKmerSmall read(final Kryo kryo, final Input input, final Class<SVKmerSmall> klass ) {
            return new SVKmerSmall(kryo, input);
        }
    }

}
