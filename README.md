# bin-to-dna
Encoding Algorithm described by Bornholt et al. in "A DNA-Based Archival Storage System".

This is a Rust implementation of the algorithm for transforming binary data into DNA from the paper "A DNA-Based Archival Storage System". It creates a Huffman tree with a base of 3 given an input string, and encodes a string into `Trit`s. The `Trit`s are then transformed into DNA based on the rotating nucleotide-trit table.

While the paper describes an algorithm that always outputs a `Trit` with a length of 5-6 for a given input letter, this implementation may output a much shorter `Trit` if the input string is short. This is because this implementation does not include unicodes that were not seen in the input when building the Huffman tree.

## Usage
```
$ git clone https://github.com/jeenalee/bin-to-dna.git
$ cd bin-to-dna
$ cargo run "[string to encode in DNA]"
```

## Example
```
jeena:bin-to-dna Jeena$ cargo run "hello world"
    Finished debug [unoptimized + debuginfo] target(s) in 0.0 secs
     Running `target/debug/balanced-bin-to-dna hello\ world`
These are the trits:
	[One, One, One, One, Zero, Zero, Two, Two, One, Two, One, Zero, Two, One, One, Two, One, Two, One, One, Zero, Two, One, Zero, One]
These are the nucleotides:
	[Guanine, Adenine, Guanine, Adenine, Cytosine, Guanine, Cytosine, Adenine, Guanine, Cytosine, Thymine, Adenine, Thymine, Cytosine, Thymine, Guanine, Adenine, Thymine, Cytosine, Thymine, Adenine, Thymine, Cytosine, Guanine, Adenine]
These are the trits decoded from nucleotides:
	[One, One, One, One, Zero, Zero, Two, Two, One, Two, One, Zero, Two, One, One, Two, One, Two, One, One, Zero, Two, One, Zero, One]
These are the bytes decoded from trits:
	[104, 101, 108, 108, 111, 32, 119, 111, 114, 108, 100]
This is the decoded nucleotides converted to a string:
	"hello world"
```
