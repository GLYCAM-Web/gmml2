### Overlap resolution
The overlap resolution in the glycoprotein builder is more complex than in carbohydrates, as it can involve multiple glycans as well as the option to move amino acid sidechains.

#### Default shapes and sampling
Glycoprotein overlaps are resolved in one of two ways, by
1. Using the most likely rotamers and default bond angles whenever possible. This is what we call the default shape
2. Varying the shape of glycans according to statistical likelihood of rotamers and bond angles

The same algorithm is used in either case. The difference is achieved in adjusting how far randomly generated preferred angles may deviate from the default angles, as well as how far from the preferred angle the algorithm is allowed to search for an optimal angle.

Constraining the range of search to only part of the entire rotamer range is important for promoting more realistic distributions of angles over multiple samples.

#### Overlap rejection threshold
Some systems are hard to fully resolve with 0 overlap. The overlap rejection threshold lets the algorithm consider a system fully resolved so long as the highest force acting on a single atom is below this threshold. What level is acceptable is likely application-specific.

#### Random descent
The algorithm performs a random descent, which continues until no further improvement can be found. As a first step, we generate preferred shapes for all glycans, and adjust linkages of all glycans in order to resolve any overlaps best possible Then, we figure out which glycans contain overlaps. If glycans overlap with eachother, these will be moved in concert.

Next, for each concert or glycan with overlaps, a new shape will be generated at random, and the concert/glycan linkages adjusted. If the new shape post-adjustment results in lower overlap than before, this shape will be used. Otherwise, the concert/glycan will be reset to its previous shape.
After all concerts/glycans have been iterated over, the algorithm again checks whether any have overlaps. If they do, it continues until either

1. No overlaps remain, or
2. The algorithm has run for X iterations without finding a global improvement, where X is defined as persist cycles in the input

#### Adjusting sidechains
In certain cases, sidechains close to a glycan attachment point can make resolving all overlap impossible without being moved out of the way. Enabling sidechain movement can help with this. It can also allow for a greater variety of glycan shapes in general.

This option slightly adjusts the algorithm. After setting concerts/glycan to a random shape, linkages will be adjusted to resolve overlaps with other glycans, and protein sans movable parts of their sidechains. Then, any movable sidechains now overlapping with glycans will be adjusted to the shape that results in the least overlap. When multiple shapes are possible, one is picked at random.

Lastly, the glycan linkages are adjusted as normal, resolving overlaps with other glycans and all parts of the protein.

The algorithm continues as described before, until all overlaps have been resolved or we go a number of iterations without improvement.

Users generally prefer modifying as little of the protein as possible. For this reason, sidechains are restored to their initial shape after the final iteration. Currently, all are restored, after which overlaps are checked, and overlapping sidechains moved again. What should probably be done is only restore sidechains which can be restored without causing overlap. A concern was that modified shapes of certain sidechains would block other sidechains from being restored to their initial shape in some cases, unless done in the specific order. Since the order is not saved, restoring all of them was chosen as a first version of the implementation.

#### Freezing glycosite conformations
On the note that users prefer modifying as little of the protein as possible, each glycosite linkage contains two bond angles defined by the protein (χ1, χ2), and two bond angles defined by the glycan (φ, ψ). Enabling this option will freeze both χ1 and χ2 to their initial angles, while allowing φ and ψ to move freely.

A drawback of this option is that we are not figuring out *which* conformer(s) the protein angles belong to. Thus, the φ, ψ range of motion is possibly using the wrong conformer. Another, more general complication is that the ‘N’ atom the glycan is attached to is sometimes incorrectly switched with the ‘O’ atom in the pdb file itself. We should detect both this and the intial conformer of the protein angles in a future version.

#### Deleting unresolvable glycans
Another option is to delete glycans whose overlaps cannot be fully resolved. This is done *after* the full algorithm has run its course. One glycan with overlap will be deleted, at which point the algorithm will restart, run its full course again, and make another check for deleting a glycan. Repeat until no glycans with overlap remain.
