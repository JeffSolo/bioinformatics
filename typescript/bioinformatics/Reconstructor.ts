/**
 * Module for reconstructing string problems.
 */

// TODO: possibly find a better place to put all of these, in a class somewhere?
/**
 * Break down string into all of its possible kmers
 *
 * @param strand - sequence to decompose
 * @param k - length of kmer
 *
 * @return string array of all kmers composing the strand
 *
 */
export function stringComposition(strand: string, k: number): string[] {
    const kmers: string[] = [];

    // tslint:disable-next-line: no-increment-decrement
    for (let i = 0; i <= strand.length - k; i++) {
        kmers.push(strand.slice(i, i + k));
    }

    return kmers;
 }

/**
 * Combine kmers into single string, when given in order
 *
 * @param strands - strands to compile
 *
 * @return compiled strand
 *
 */
export function genomePath(strands: string[]): string {
  if (strands.length === 0) { throw new Error('"strands" can not be empty'); }
  // @ts-ignore: won't be undefined due to previous check
  let out: string = strands.shift();

  strands.forEach((strand) => {
    const overlap = strand.slice(0, -1);
    if (out.endsWith(overlap)) {
      out += strand.slice(-1);
    }
  });

  return out;
}
