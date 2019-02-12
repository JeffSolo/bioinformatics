export default class DNA {
    private static complement: object = {
        A: 'T',
        T: 'A',
        C: 'G',
        G: 'C',
    };

    public sequence: string;

    constructor(sequence: string) {
        this.sequence = sequence;
    }

    public complement(): string {
        let comp: string = '';
        const array: string[] = this.sequence.split('');

        [...array].forEach((c: string) => comp += DNA.complement[c]);

        return comp;
    }

    public reverse_complement() : string {
        return this.complement().split('').reverse().join('');
    }
}