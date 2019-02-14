/**
 * Module for reconstructing string problems.
 */

/**
 * Break down string into all of its possible kmers
 *
 * @param strand - sequence to decompose
 * @param k - lenght of kmer
 *
 * @return string array of all kmers composing the strand
 *
 */
// TODO find a better place to put this - don't want it in DNA, at least for now.
export function stringComposition(strand: string, k: number): string[] {
    const kmers: string[] = [];

    // tslint:disable-next-line: no-increment-decrement
    for (let i = 0; i <= strand.length - k; i++) {
        kmers.push(strand.slice(i, i + k));
    }

    return kmers;
 }
