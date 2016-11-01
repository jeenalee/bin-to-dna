extern crate min_max_heap;

use min_max_heap::*;
use std::cmp::Ordering;
use std::env;
use std::collections::HashMap;
use std::str;

/// `BalancedTree` is a type that uses the weight to create a balanced
/// Huffman encoding tree.
#[derive(Debug, Eq)]
struct BalancedTree {
    weight: usize,
    tree: Tree,
}

/// `Tree` can be either a Node or a Leaf.
#[derive(Debug, PartialEq, Eq)]
enum Tree {
    /// `Node` holds three branches as this is a Huffman encoding tree
    /// with a base of three.
    Node {
        two: Option<Box<BalancedTree>>,
        one: Option<Box<BalancedTree>>,
        zero: Option<Box<BalancedTree>>,
    },
    /// `Leaf` is the end of a tree branch and holds the byte for the
    /// path from the root to the leaf of the Huffman encoding tree.
    Leaf(u8),
}

/// `Trit` is like a bit, but for a base of three.
#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
enum Trit {
    Two,
    One,
    Zero,
}

/// `Nucleotide` is the four variants that make up DNA.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Nucleotide {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
}


fn main() {
    if let Some(to_encode) = env::args().nth(1) {
        let tree = BalancedTree::new(&to_encode);
        let encoding_map = tree.encoding_map();
        let trits = encode_to_trits(&encoding_map, &to_encode);
        println!("These are the trits:\n\t{:?}", trits);
        let nucleotides = encode_to_nucleotides(&trits);
        println!("These are the nucleotides:\n\t{:?}", nucleotides);
        let decoded_nucleotides = decode_nucleotides(&nucleotides).expect(
            "Nucleotides were not decoded.");
        println!("These are the trits decoded from nucleotides:\n\t{:?}", decoded_nucleotides);
        let decoded_trits = tree.decode(&decoded_nucleotides).expect(
            "Trits were not decoded.");
        println!("These are the bytes decoded from trits:\n\t{:?}", decoded_trits);
        let decoded_str = str::from_utf8(&decoded_trits).expect(
            "u8 was not converted into a str.");
        println!("This is the decoded nucleotides converted to a string:\n\t{:?}", decoded_str);
    }
}

/// Takes a slice of bytes and returns a vector of length of
/// std::u8::MAX that maps a byte to count.
fn count_bytes(bytes: &[u8]) -> Vec<usize> {
    let mut counts = vec![0; std::u8::MAX as usize];
    for &b in bytes {
        counts[b as usize] += 1;
    }
    counts
}

/// Removes tuples in which the byte count is zero, and returns a
/// sorted vector in which the tuples are sorted by count from
/// smallest to largest.
fn label_and_filter_out_zero(v: Vec<usize>) -> Vec<(u8, usize)> {
    let mut processed_vec: Vec<_> = v.into_iter().enumerate()
        .map(|(index, value)| (index as u8, value))
        .filter(|&(_, x)| x != 0).collect();
    sort_by_byte_count(&mut processed_vec);
    processed_vec
}

/// Sorts a vector by byte count and returns the same vector.
fn sort_by_byte_count(i: &mut Vec<(u8, usize)>) {
    i.sort_by(|a, b| {
        if a.1 < b.1 {
            // Sort from smallest to largest.
            std::cmp::Ordering::Less
        } else if a.1 > b.1 {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });
}

impl BalancedTree {
    /// Creates a `BalancedTree` from an input.
    fn new(input: &str) -> BalancedTree {
        let b = input.as_bytes();
        let c = count_bytes(b);
        let v = label_and_filter_out_zero(c);
        let mut heap = create_min_max_heap(v);
        BalancedTree::new_helper(&mut heap)
    }

    /// Recursively builds a balanced Huffman encoding tree with the
    /// base of three.
    fn new_helper(trees: &mut MinMaxHeap<BalancedTree>) -> BalancedTree {
        assert!(trees.len() >= 1);
        if trees.len() == 1 {
            trees.pop_min().unwrap()
        } else {
            let min1 = trees.pop_min();
            let min2 = trees.pop_min();
            let min3 = trees.pop_min();

            let get_weight = |t: &Option<BalancedTree>| t.as_ref().map_or_else(|| 0, |t| t.weight);
            let weight = get_weight(&min1) + get_weight(&min2) + get_weight(&min3);

            let tree = BalancedTree {
                weight: weight,
                tree: Tree::Node {
                    two: match min1 {
                        Some(m) => Some(Box::new(m)),
                        None => None,
                    },
                    one: match min2 {
                        Some(m) => Some(Box::new(m)),
                        None => None,
                    },
                    zero: match min3 {
                        Some(m) => Some(Box::new(m)),
                        None => None,
                    },
                }
            };
            trees.push(tree);
            BalancedTree::new_helper(trees)
        }
    }

    /// Makes an encoding map based on a Huffman tree.
    fn encoding_map(&self) -> HashMap<u8, Vec<Trit>> {
        let mut code = Vec::new();
        let mut tree_map = HashMap::new();
        self.encoding_map_helper(&mut code, &mut tree_map);
        tree_map
    }

    /// Recursively builds an encoding tree map by walking each node
    /// in a Huffman tree.
    fn encoding_map_helper(&self, code: &mut Vec<Trit>, map: &mut HashMap<u8, Vec<Trit>>) {
        match self {
            &BalancedTree { ref tree, ..} => match tree {
                &Tree::Leaf(byte) => { map.insert(byte, code.clone()); },
                &Tree::Node { ref two, ref one, ref zero }  => {
                    // Visit Node.two.
                    if let &Some(ref t) = two {
                        code.push(Trit::Two);
                        t.encoding_map_helper(code, map);
                        code.pop();
                    }
                    if let &Some(ref o) = one {
                        // Visit Node.one.
                        code.push(Trit::One);
                        o.encoding_map_helper(code, map);
                        code.pop();
                    }
                    if let &Some(ref z) = zero {
                        // Visit Node.zero.
                        code.push(Trit::Zero);
                        z.encoding_map_helper(code, map);
                        code.pop();
                    }
                },
            }
        }
    }

    /// Decodes the slice of `Trit`s by walking the tree, and returns
    /// the decoded byte and the rest of the `Trit`s.
    fn decode_byte<'a, 'b>(&'a self, trits: &'b [Trit]) -> Result<(u8, &'b [Trit]), &'static str> {
        match self.tree {
            Tree::Leaf(byte) => Ok((byte, trits)),
            Tree::Node { ref two, ref one, ref zero } => {
                match (trits.first().cloned(), two, one, zero) {
                    (Some(Trit::Two), &Some(ref two), _, _) => two.decode_byte(&trits[1..]),
                    (Some(Trit::One), _, &Some(ref one), _) => one.decode_byte(&trits[1..]),
                    (Some(Trit::Zero), _, _, &Some(ref zero)) => zero.decode_byte(&trits[1..]),
                    (None, _, _, _) => Err(("Unexpectedly reached the end of the tree.")),
                    _ => Err(("Bad encoding of trits. The trit doesn't correspond to bytes in the tree.")),

                }
            }
        }
    }

    /// Iterates through the given trits and returns the vector of
    /// decoded bytes.
    fn decode<'a, 'b>(&'a self, mut trits: &'b [Trit]) -> Result<Vec<u8>, &'static str> {
        let mut decoded = Vec::new();
        while trits.len() > 0 {
            let (byte, remaining_trits) = try!(self.decode_byte(trits));
            decoded.push(byte);
            trits = remaining_trits;
        }
        Ok(decoded)
    }
}

/// Based on the input and the encoding tree map, returns the encoded
/// input.
fn encode_to_trits(encoding_map: &HashMap<u8, Vec<Trit>>, input: &String) -> Vec<Trit> {
    let bytes = input.as_bytes();
    let mut output = Vec::new();
    for byte in bytes {
        // Map needs to retain ownership. So we iterate through it and clone it.
        output.extend(encoding_map[byte].iter().cloned());
    }
    output
}

/// Creates the initial `MinMaxHeap` of the `BalancedTree`s, in which
/// each byte is a `Leaf`.
fn create_min_max_heap(v: Vec<(u8, usize)>) -> MinMaxHeap<BalancedTree> {
    let mut trees: MinMaxHeap<BalancedTree> = MinMaxHeap::new();
    for (byte, count) in v {
        let n = BalancedTree {
            weight: count,
            tree: Tree::Leaf(byte),
        };
        trees.push(n);
    }
    trees
}

/// Finds the current `Nucleotide` based on a given `Trit` and the previous
/// `Nucleotide`.
fn find_current_nucleotide(trit: Trit, previous_nucleotide: Nucleotide)
                        -> Nucleotide {
    use Trit::*;
    use Nucleotide::*;

    match (previous_nucleotide, trit) {
        (Adenine, Zero) => Cytosine,
        (Adenine, One) => Guanine,
        (Adenine, Two) => Thymine,

        (Cytosine, Zero) => Guanine,
        (Cytosine, One) => Thymine,
        (Cytosine, Two) => Adenine,

        (Guanine, Zero) => Thymine,
        (Guanine, One) => Adenine,
        (Guanine, Two) => Cytosine,

        (Thymine, Zero) => Adenine,
        (Thymine, One) => Cytosine,
        (Thymine, Two) => Guanine,
    }
}

/// Encodes `Trit`s to `Nucleotide`s.
fn encode_to_nucleotides(trits: &[Trit]) -> Vec<Nucleotide> {
    let mut nucleotides = Vec::new();

    // The paper does not mention what the first `Nucleotide` will
    // be. This algorithm assumes the first `Nucleotide`'s previous
    // `Nucleotide` is Adenine.
    let mut prev_nucleotide = Nucleotide::Adenine;

    for &trit in trits {
        let current_nucleotide = find_current_nucleotide(trit, prev_nucleotide);
        nucleotides.push(current_nucleotide);
        prev_nucleotide = current_nucleotide;
    }
    nucleotides
}

/// Finds the current `Trit` based on a given `Nucleotide` and the
/// previous `Nucleotide`.
fn find_current_trit(current_nucleotide: Nucleotide, previous_nucleotide: Nucleotide)
                     -> Result<Trit, ()> {
    use Trit::*;
    use Nucleotide::*;

    match (previous_nucleotide, current_nucleotide) {
        (Adenine, Cytosine) => Ok(Zero),
        (Adenine, Guanine) => Ok(One),
        (Adenine, Thymine) => Ok(Two),

        (Cytosine, Guanine) => Ok(Zero),
        (Cytosine, Thymine) => Ok(One),
        (Cytosine, Adenine) => Ok(Two),

        (Guanine, Thymine) => Ok(Zero),
        (Guanine, Adenine) => Ok(One),
        (Guanine, Cytosine) => Ok(Two),

        (Thymine, Adenine) => Ok(Zero),
        (Thymine, Cytosine) => Ok(One),
        (Thymine, Guanine) => Ok(Two),

        _ => Err(()),
    }
}

/// Decodes `Nucleotide`s to `Trit`s.
fn decode_nucleotides(nucleotides: &[Nucleotide]) -> Result<Vec<Trit>, ()> {
    let mut trits = Vec::new();

    // The paper does not mention what the first `Nucleotide` will
    // be. This algorithm assumes the first `Nucleotide`'s previous
    // `Nucleotide` is Adenine.
    let mut prev_nucleotide = Nucleotide::Adenine;

    for &nucleotide in nucleotides {
        let current_trit = try!(find_current_trit(nucleotide, prev_nucleotide));
        trits.push(current_trit);
        prev_nucleotide = nucleotide;
    }
    Ok(trits)
}

impl Ord for BalancedTree {
    fn cmp(&self, other: &BalancedTree) -> Ordering {
        self.weight.cmp(&other.weight)
    }
}

impl PartialOrd for BalancedTree {
    fn partial_cmp(&self, other: &BalancedTree) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for BalancedTree {
    fn eq(&self, other: &BalancedTree) -> bool {
        self.weight == other.weight
    }
}
