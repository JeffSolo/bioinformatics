// import DNA from '../bioinformatics/DNA';
import {IDataFile, parseDataFile, saveToFile} from '../bioinformatics/CourseHelper';
import {genomePath, stringComposition} from '../bioinformatics/Reconstructor';

const dataPath = './datasets/week1';
const outPath = './output/week1';

const {body: sequence, header: k}: IDataFile = parseDataFile(`${dataPath}/dataset_197_3.txt`, {hasHeader: true});
saveToFile(`${outPath}/stringComp.txt`, stringComposition(sequence, Number(k)), true);

const {body: strands}: IDataFile = parseDataFile(`${dataPath}/dataset_198_3.txt`, {});
saveToFile(`${outPath}/genomePath.txt`, genomePath(strands.split(' ')), true);
