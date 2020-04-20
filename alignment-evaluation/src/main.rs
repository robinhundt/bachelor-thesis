use std::ffi::OsStr;
use std::fmt::Display;
use std::fs::{DirEntry, File};
use std::iter::FromIterator;
use std::ops::{Add, Not};
use std::path::PathBuf;
use std::process::Command;
use std::sync::Mutex;
use std::time::{Duration, Instant};
use std::{fmt, fs};

use anyhow::Result;
use bali_score::{balibase, compute_scores, fasta, AlignmentScores};
use indicatif::ProgressIterator;
use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use spam_align::align::{align, AlignProgress};
use spam_align::spaced_word::{read_patterns_from_file, Pattern};
use spam_align::{read_fasta, write_as_fasta};

fn main() -> Result<()> {
    // compute_results_for_balibase(
    //     AlignmentProgram::MafftAuto,
    //     &["--quiet", "--auto"],
    //     vec![],
    //     "../datasets/bb3_release".into(),
    //     "../evaluation-data/".into(),
    // )
    // .unwrap();
    // compute_results_for_balibase(
    //     AlignmentProgram::MafftFast,
    //     &["--quiet", "--retree", "1", "--maxiterate", "0"],
    //     vec![],
    //     "../datasets/bb3_release".into(),
    //     "../evaluation-data/".into(),
    // )
    // .unwrap();
    //
    // compute_results_for_balibase(
    //     AlignmentProgram::Dialign,
    //     &["-fa"],
    //     vec![("DIALIGN2_DIR", "../dialign_package/dialign2_dir")],
    //     "../datasets/bb3_release".into(),
    //     "../evaluation-data/".into(),
    // )
    // .unwrap();

    // for pattern_set_path in fs::read_dir("../pattern_sets/data")?.progress_cunt(31) {
    //     let pattern_set_path = pattern_set_path?;
    //     let pattern_set = read_patterns_from_file(pattern_set_path.path())?;
    //     let pattern_set = PatternSet {
    //         patterns: pattern_set,
    //         name: pattern_set_path
    //             .path()
    //             .file_stem()
    //             .unwrap()
    //             .to_string_lossy()
    //             .to_owned()
    //             .to_string(),
    //     };
    //     compute_results_for_balibase(
    //         AlignmentProgram::SpamAlign(pattern_set),
    //         &vec![],
    //         vec![],
    //         "../datasets/bb3_release".into(),
    //         "../evaluation-data/".into(),
    //     )
    //     .unwrap();
    // }

    Ok(())
}

enum AlignmentProgram {
    MafftAuto,
    MafftFast,
    Dialign,
    SpamAlign(PatternSet),
}

struct PatternSet {
    patterns: Vec<Pattern>,
    name: String,
}

impl Display for AlignmentProgram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignmentProgram::MafftAccurate => write!(f, "Mafft-Accurate"),
            AlignmentProgram::MafftFast => write!(f, "Mafft-Fast"),
            AlignmentProgram::Dialign => write!(f, "Dialign"),
            AlignmentProgram::SpamAlign(pattern_set) => {
                write!(f, "spam-align-max-dim-{}", pattern_set.name)
            }
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct EvalResult {
    #[serde(flatten)]
    scores: AlignmentScores,
    time_ms: u128,
}

fn compute_results_for_balibase(
    program: AlignmentProgram,
    args: &[&str],
    envs: Vec<(&str, &str)>,
    balibase_path: PathBuf,
    out_path: PathBuf,
) -> Result<()> {
    let balibase_folders = ["RV11", "RV12", "RV20", "RV30", "RV40", "RV50"];
    // let balibase_folders = ["RV20"];
    let balibase_folders = balibase_folders.iter().map(|folder| {
        let mut path = balibase_path.clone();
        path.push(folder);
        path
    });
    for folder in balibase_folders {
        let mut out_folder_path = out_path.clone();
        out_folder_path.push(program.to_string());
        out_folder_path.push(folder.file_name().unwrap());

        fs::create_dir_all(&out_folder_path)?;
        let mut bb_files: (Vec<_>, Vec<_>) = fs::read_dir(folder)?
            .map(Result::unwrap)
            .map(|dir_entry: DirEntry| dir_entry.path())
            .filter(|path| {
                let fasta_or_xml = path.extension() == Some(OsStr::new("tfa"))
                    || path.extension() == Some(OsStr::new("xml"));
                let is_not_bbs = path
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with("BBS")
                    .not();
                fasta_or_xml && is_not_bbs
            })
            .partition(|path| path.extension() == Some(OsStr::new("tfa")));

        bb_files.0.sort();
        bb_files.1.sort();
        let bb_files = bb_files;

        let eval_results = Mutex::new(vec![]);

        bb_files.par_iter().for_each(|(fasta_file, xml_file)| {
            let mut out_path = out_folder_path.clone();
            out_path.push(fasta_file.file_name().unwrap());
            let duration = match &program {
                AlignmentProgram::MafftAuto | AlignmentProgram::MafftFast => {
                    run_mafft(fasta_file, &out_path, args).unwrap()
                }
                AlignmentProgram::Dialign => {
                    run_dialign(fasta_file, &out_path, args, envs.clone()).unwrap()
                }
                AlignmentProgram::SpamAlign(pattern_set) => {
                    run_spam_align(fasta_file, &out_path, pattern_set).unwrap()
                }
            };
            let test_alignment = fasta::parse(&out_path).unwrap();
            let ref_alignment = balibase::parse(xml_file).unwrap();
            let scores = compute_scores(&ref_alignment, &test_alignment, false);
            let eval_result = EvalResult {
                scores,
                time_ms: duration.as_millis(),
            };

            eval_results.lock().unwrap().push(eval_result.clone());
            let eval_results_file_name =
                format!("{}.json", out_path.file_name().unwrap().to_str().unwrap());
            out_path.set_file_name(eval_results_file_name);
            let eval_results_file = File::create(out_path).unwrap();
            serde_json::to_writer_pretty(eval_results_file, &eval_result).unwrap();
        });
        let mut aggr_eval_results = eval_results
            .into_inner()
            .unwrap()
            .into_iter()
            .fold(EvalResult::default(), |aggr, el| aggr + el);
        let num_files = bb_files.0.len();
        aggr_eval_results.time_ms /= num_files as u128;
        aggr_eval_results.scores.column_score /= num_files as f64;
        aggr_eval_results.scores.sum_of_pairs /= num_files as f64;
        out_folder_path.push("aggr.json");
        let aggr_file = File::create(out_folder_path)?;
        serde_json::to_writer_pretty(aggr_file, &aggr_eval_results)?;
    }
    Ok(())
}

fn run_mafft(fasta_in: &PathBuf, fasta_out: &PathBuf, args: &[&str]) -> Result<Duration> {
    let out_file = File::create(fasta_out)?;
    let now = Instant::now();
    let exit_status = Command::new("mafft")
        .args(args)
        .arg(fasta_in)
        .stdout(out_file)
        .spawn()?
        .wait()?;
    let alignment_time = now.elapsed();
    assert!(exit_status.success());
    Ok(alignment_time)
}

fn run_dialign(
    fasta_in: &PathBuf,
    fasta_out: &PathBuf,
    args: &[&str],
    envs: Vec<(&str, &str)>,
) -> Result<Duration> {
    let now = Instant::now();
    let exit_status = Command::new("../dialign_package/src/dialign2-2")
        .args(args)
        .arg("-fn")
        .arg(fasta_out)
        .arg(fasta_in)
        .envs(envs)
        .spawn()?
        .wait()?;
    let alignment_time = now.elapsed();
    assert!(exit_status.success());
    Ok(alignment_time)
}

fn run_spam_align(
    fasta_in: &PathBuf,
    fasta_out: &PathBuf,
    pattern_set: &PatternSet,
) -> Result<Duration> {
    let now = Instant::now();
    let mut input = read_fasta(fasta_in)?;
    align(&mut input, &pattern_set.patterns, AlignProgress::Hide);
    write_as_fasta(fasta_out, &input)?;
    Ok(now.elapsed())
}

impl Add for EvalResult {
    type Output = EvalResult;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.scores = self.scores + rhs.scores;
        self.time_ms += rhs.time_ms;
        self
    }
}
