import * as fs from 'fs';

export function parseDataFile(
    filePath: string,
    hasHeader: boolean = false,
    hasFooter: boolean = false,
    joinCharacter: string = ''
): object {
  // tslint:disable-next-line:non-literal-fs-path
  const fileContent: string[] = fs.readFileSync(filePath).toString().split('\n');
  const header: string | undefined = hasHeader ? fileContent.shift() : undefined;
  const footer: string | undefined = hasFooter ? fileContent.pop() : undefined;
  const body: string = fileContent.join(joinCharacter);

  return {body, header, footer};
}

// We're just printing output, so don't care about type here
// tslint:disable-next-line:no-any
export function printFormattedOutput(answer: any[] | any, joiner: string = ' '): void {
  if (typeof(answer) === 'object') {
    console.log(answer.join(joiner));
  } else {
    console.log(answer);
  }
}

export function parseParameters(parametersString: string, splitCharacter: string = ' '): string[] {
  return parametersString.split(splitCharacter);
}
