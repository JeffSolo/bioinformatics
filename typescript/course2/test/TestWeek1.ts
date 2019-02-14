import assert from 'assert';
import {stringComposition} from '../../bioinformatics/Reconstructor';

describe('Week 1 test', () => {
    it('decompose strand', () => {
        const input = 'CAATCCAAC';
        const size = 5;
        const expected = ['CAATC', 'AATCC', 'ATCCA', 'TCCAA', 'CCAAC'];
        assert.deepEqual(stringComposition(input, size), expected);
    });
});
