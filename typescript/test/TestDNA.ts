import assert from 'assert';
import DNA from '../bioinformatics/DNA';

describe('DNA class', () => {
    it('should get DNA complement', () => {
        const dna: DNA = new DNA('ACGT');
        assert.equal(dna.complement(), 'TGCA');
    });

    it('should get DNA reverse complement', () => {
       const dna: DNA = new DNA('ACGTTCGA');
       assert.equal(dna.reverse_complement(), 'TCGAACGT');
    });

});
