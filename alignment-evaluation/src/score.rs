use serde::{Serialize, Deserialize};
use spam_align::{Alignment, PositionAlignment};
use itertools::Itertools;
use spam_align::align::micro_alignment::Site;
use num_integer::binomial;
use std::ops::Add;

#[derive(Default, PartialEq, Copy, Clone, Debug, Serialize, Deserialize)]
pub struct AlignmentScores {
    pub sum_of_pairs: f64,
    pub column_score: f64,
}


pub fn compute_scores(
    ref_alignment: &Alignment,
    test_alignment: &Alignment,
    ignore_symbol_case: bool,
) -> AlignmentScores {
    let seq_cnt = test_alignment.aligned_data.len();
    let mut correctly_aligned = 0;
    let mut correct_column_count = 0;
    let test_column_cnt = test_alignment
        .aligned_data
        .first()
        .expect("Empty alignment")
        .data
        .len();

    let mut gaps_per_seq = vec![0; seq_cnt];
    let mut gaps_per_col = vec![0; test_column_cnt];
    let mut correct_in_column = vec![false; test_column_cnt];

    for column in 0..test_column_cnt {
        gaps_per_col.iter_mut().for_each(|el| *el = 0);
        correct_in_column.iter_mut().for_each(|el| *el = false);
        for (seq_id, seq) in test_alignment.aligned_data.iter().enumerate() {
            if seq.data[column] == Alignment::GAP_CHARACTER {
                gaps_per_seq[seq_id] += 1;
                gaps_per_col[seq_id] += 1;
            }
        }
        for ((a_idx, seq_a), (b_idx, seq_b)) in test_alignment
            .aligned_data
            .iter()
            .enumerate()
            .tuple_combinations()
        {
            let a = seq_a.data[column];
            let b = seq_b.data[column];

            if a == Alignment::GAP_CHARACTER || b == Alignment::GAP_CHARACTER {
                continue;
            }

            let site_a = Site {
                seq: a_idx,
                pos: column - gaps_per_seq[a_idx],
            };
            let site_b = Site {
                seq: b_idx,
                pos: column - gaps_per_seq[b_idx],
            };

            let symbols_aligned =
                ignore_symbol_case || a.is_ascii_uppercase() && b.is_ascii_uppercase();

            if symbols_aligned
                && ref_alignment.pos_aligned(site_a, site_b) == PositionAlignment::Correct
            {
                correct_in_column[a_idx] = true;
                correct_in_column[b_idx] = true;
                correctly_aligned += 1;
            }
        }
        let correct_in_column_sum: usize = correct_in_column.iter().copied().map(usize::from).sum();
        // TODO Note that this only works for ref alignments with core blocks where every seq
        // is aligned in this position
        if correct_in_column_sum == seq_cnt {
            correct_column_count += 1;
        }
    }

    let sum_of_pairs = correctly_aligned as f64 / true_site_pair_count(ref_alignment) as f64;
    let column_score = correct_column_count as f64 / ref_alignment.core_block_columns() as f64;
    AlignmentScores {
        sum_of_pairs,
        column_score,
    }
}


fn true_site_pair_count(alignment: &Alignment) -> usize {
    let core_block_data = alignment.core_block_data();
    let aligned_seq_len = core_block_data[0].data.len();
    let aligned_pos_per_col = (0..aligned_seq_len).map(|pos| {
        core_block_data.iter().fold(
            0,
            |acc, seq| {
                if seq.data[pos] == b'-' {
                    acc
                } else {
                    acc + 1
                }
            },
        )
    });
    let max_site_pair_count =
        aligned_pos_per_col.fold(0, |acc, aligned_count| acc + binomial(aligned_count, 2));

    max_site_pair_count
}


impl Add for AlignmentScores {
    type Output = AlignmentScores;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.sum_of_pairs += rhs.sum_of_pairs;
        self.column_score += rhs.column_score;
        self
    }
}