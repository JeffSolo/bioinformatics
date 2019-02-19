import assert from 'assert';
import fs from 'fs';
import {IDataFile, parseDataFile} from '../../bioinformatics/CourseHelper';
import {stringComposition} from '../../bioinformatics/Reconstructor';

const baseTestInputPath = './typescript/course2/test/inputs';
const baseTestOutputPath = './typescript/course2/test/outputs';
let inputPath = '';
let outputPath = '';

describe('Week 1 Tests', () => {
  describe('String Composition Tests', () => {
    before(() => {
      inputPath =  `${baseTestInputPath}/kmerComposition`;
      outputPath = `${baseTestOutputPath}/kmerComposition`;
    });

    ['sample.txt', 'test1.txt', 'test2.txt', 'test3.txt', 'test4.txt'].forEach((fileName) => {
      it(`should decompose strand in file ${fileName}`, () => {
        const {body, header}: IDataFile = parseDataFile(`${inputPath}/${fileName}`, {hasHeader: true});
        const actualOutput = stringComposition(body, Number(header));
        const expectedOutput = fs.readFileSync(`${outputPath}/${fileName}`).toString().split('\n');

        assert.deepEqual(actualOutput, expectedOutput);
      });
    });
  });
});
