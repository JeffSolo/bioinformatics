import * as fs from 'fs';

/**
 * Parses data file into body, header, and footer.
 *
 * @remarks Header and footer can only be one line.
 *
 * @parameter filePath - File we want to read data from
 * @parameter asHeader - Whether or not file has a header
 * @parameter hasFooter - Whether or not file has a footer, assumes 1 line
 * @parameter join_character - What character to separate each body line with when turning into a single string
 *
 * @returns object consisting of {body, header, footer}
 */
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

/**
 * Format and print output to console.
 *
 * @parameter output - Info we want to print
 * @parameter joiner - Spacing character for printing array items
 */
// tslint:disable-next-line:no-any
export function printFormattedOutput(output: any[] | any, joiner: string = ' '): void {
  if (typeof(output) === 'object') {
    console.log(output.join(joiner));
  } else {
    console.log(output);
  }
}

/**
 * Parse a string of parameters into array.
 *
 * @parameter parameterString - Info we want to print
 * @parameter splitCharacter - Spacing character for printing array items
 *
 * @returns String array of parameters
 */
export function parseParameters(parametersString: string, splitCharacter: string = ' '): string[] {
  return parametersString.split(splitCharacter);
}
