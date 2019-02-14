import * as fs from 'fs';

export function parseDataFile(
    filePath: string,
    hasHeader: boolean = false,
    hasFooter: boolean = false,
    joinCharacter: string = ''
): object {
  // tslint:disable-next-line:non-literal-fs-path
  const fileContent: string[] = fs.readFileSync(filePath).toString().split('\n');
  // TODO ensure that length of array is between > 0 and <= 3
  const footer: string | undefined = hasFooter ? fileContent.pop() : undefined;
  const body: string | undefined = fileContent.pop();
  const header: string | undefined = hasHeader ? fileContent.pop() : undefined;
  // TODO could check that array is empty

  return {body, header, footer};
}

// Since we're just printing output, let it be an array of anything, or a string, or a number
// tslint:disable-next-line:no-any
export function printFormattedOutput(answer: any[] | string | number, joiner: string = ' '): void {
  if (typeof(answer) === 'object') {
    console.log(answer.join(joiner));
  } else {
    console.log(answer);
  }
}

export function parseParameters(parametersString: string, splitCharacter: string = ' '): string[] {
  return parametersString.split(splitCharacter);
}
