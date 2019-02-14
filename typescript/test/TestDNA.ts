import assert from 'assert';
import DNA from '../bioinformatics/DNA';

describe('DNA class', () => {
    it('should decompose strand', () => {
        const dna = new DNA('ACGT');
        assert.strictEqual(dna.complement(), 'TGCA');
    });

    it('should get DNA reverse complement', () => {
       const dna = new DNA('ACGTTCGA');
       assert.strictEqual(dna.reverse_complement(), 'TCGAACGT');
    });

});
