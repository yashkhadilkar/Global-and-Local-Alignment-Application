
# Global and Local Alignment Application

A comprehensive Python application for performing global and local sequence alignments using dynamic programming algorithms. The application features a user-friendly graphical interface built with Tkinter.




## Problem Statement

Sequence alignment is a fundamental problem in bioinformatics where we need to compare biological sequences (DNA, RNA, or proteins) to identify similarities, differences, and evolutionary relationships.

However, most existing tools: 

- Lack educational value: Commercial tools are black boxes that don't show the underlying algorithm
- Don't show the process: Users see results but not how the dynamic programming matrix is built
## Solution

This application solves these problems by providing:

- Visual learning tool: See exactly how the dynamic programming matrix is filled and how optimal paths are traced
- Educational transparency: Perfect for teaching bioinformatics algorithms or understanding how sequence alignment works
- User-friendly interface: No command-line expertise required
## Alignment Types

- Global Alignment: Implements the Needleman-Wunsch algorithm for end-to-end sequence alignment
- Local Alignment: Implements the Smith-Waterman algorithm for finding optimal local alignments
## Installation Prerequisites


- Python 3.8 or higher
- NumPy
- Tkinter (usually included with Python)
## Usage

1. Launch the application: Run python tkinter_example.py
2. Enter sequences: Input your nucleotide sequences in the provided text fields
3. Set scoring parameters:
- Match Score: Points awarded for matching characters
- Mismatch Score: Points deducted for non-matching characters
- Gap Score: Penalty for introducing gaps
4. Choose alignment type: Click either "Global Align" or "Local Align"
5. View results: The application displays:
- Color-coded alignment matrix with optimal path highlighted in blue
- Aligned sequences showing matches, mismatches, and gaps

## Scoring Parameters
### Typical Values
- Match Score: +1 to +5 (positive values)
- Mismatch Score: -1 to -3 (negative values)
- Gap Penalty: -1 to -5 (negative values)

### Guidelines
- Higher match scores favor longer alignments
- More negative mismatch scores penalize substitutions
- More negative gap penalties discourage insertions/deletions
## Contributing

Contributions are always welcome! Please feel free to submit pull requests.

## License

[MIT](https://choosealicense.com/licenses/mit/)


## Acknowledgements

 - Dynamic programming principles in bioinformatics
 - Arnav Amin
