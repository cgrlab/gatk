package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Iterator over successive Kmers from a sequence of characters.
 * Silently skips over parts of the sequence that has characters other than A, C, G, or T.
 */
public class SVKmerizer<KType extends SVKmerBase> implements Iterator<KType> {
    protected final CharSequence seq;
    protected final int kSize;
    protected int idx = 0;
    protected KType nextKmer;

    public SVKmerizer( final byte[] seq, final int kSize, KType kmer ) {
        this(new ASCIICharSequence(seq), kSize, kmer );
    }

    public SVKmerizer( final CharSequence seq, final int kSize, KType kmer ) {
        this.seq = seq;
        this.kSize = kSize;
        this.nextKmer = nextKmer(kmer, 0);
    }

    protected SVKmerizer( final int kSize, final CharSequence seq, KType kmer ) {
        this.seq = seq;
        this.kSize = kSize;
    }

    @Override
    public boolean hasNext() {
        return nextKmer != null;
    }

    @Override
    public KType next() {
        if ( nextKmer == null ) throw new NoSuchElementException("Kmerization sequence exhausted.");
        final KType result = nextKmer;
        nextKmer = nextKmer(nextKmer, kSize-1);
        return result;
    }

    public static <KType extends SVKmerBase> KType toKmer( final CharSequence seq, KType kmer ) {
        final SVKmerizer<KType> sk = new SVKmerizer<KType>(seq, seq.length(), kmer);
        Utils.validateArg(sk.hasNext(), () -> "Can't make a SVKmer from '"+seq+"'");
        return sk.next();
    }

    public static <KType extends SVKmerBase> KType toKmer( final byte[] seq, KType kmer ) {
        return toKmer(new ASCIICharSequence(seq),kmer);
    }

    public static <KType extends SVKmerBase> Stream<KType> stream( final CharSequence seq, final int kSize, KType kmer ) {
        return StreamSupport.stream(((Iterable<KType>)() -> new SVKmerizer<KType>(seq, kSize, kmer)).spliterator(), false);
    }

    public static <KType extends SVKmerBase> Stream<KType> stream( final byte[] seq, final int kSize, KType kmer ) {
        return stream(new ASCIICharSequence(seq),kSize,kmer);
    }

    @SuppressWarnings("unchecked")
    protected KType nextKmer( KType tmpKmer, int validBaseCount ) {
        final int len = seq.length();
        while ( idx < len ) {
            switch ( seq.charAt(idx) ) {
                case 'a': case 'A': tmpKmer = (KType)tmpKmer.successor(SVKmerBase.Base.A, kSize); break; //TODO: unchecked casts
                case 'c': case 'C': tmpKmer = (KType)tmpKmer.successor(SVKmerBase.Base.C, kSize); break;
                case 'g': case 'G': tmpKmer = (KType)tmpKmer.successor(SVKmerBase.Base.G, kSize); break;
                case 't': case 'T': tmpKmer = (KType)tmpKmer.successor(SVKmerBase.Base.T, kSize); break;
                default: validBaseCount = -1;
            }
            idx += 1;

            if ( ++validBaseCount == kSize ) return tmpKmer;
        }
        return null;
    }

    // a shim to turn a byte array into a character sequence by treating the bytes as ASCII characters
    public static final class ASCIICharSequence implements CharSequence {
        ASCIICharSequence( final byte[] bytes ) { this.bytes = bytes; }

        @Override public int length() { return bytes.length; }

        @Override public char charAt( final int index ) { return (char)(bytes[index] & 0xff); }

        @Override public CharSequence subSequence( final int start, final int end ) {
            return new ASCIICharSequence(Arrays.copyOfRange(bytes, start, end));
        }

        @Override public String toString() { return new StringBuilder(this).toString(); }

        final byte[] bytes;
    }
}
