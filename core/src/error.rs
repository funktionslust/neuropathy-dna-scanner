//! Error enum + stable exit code mapping.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("input file is not a valid CRAM or BAM: {0}")]
    InvalidAlignmentFormat(String),

    #[error("required index ({0}) is missing or unreadable")]
    MissingIndex(String),

    #[error("reference build mismatch: CRAM expects GRCh38, got {got}")]
    ReferenceBuildMismatch { got: String },

    #[error(
        "reference MD5 mismatch on chr17: CRAM header expects {expected}, reference has {got}"
    )]
    ReferenceMd5Mismatch { expected: String, got: String },

    #[error("input file is truncated or missing the EOF marker")]
    TruncatedInput,

    #[error("zero mapped reads in {region}")]
    EmptyRegion { region: String },

    #[error("invalid configuration: {0}")]
    InvalidConfig(String),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Genuinely unexpected internal-state errors. Not a catch-all: add a
    /// specific variant when one would be more informative.
    #[error("internal error: {0}")]
    Internal(String),
}

impl Error {
    /// Exit code for this error variant. Part of the public contract.
    pub fn exit_code(&self) -> i32 {
        match self {
            Error::InvalidAlignmentFormat(_) => 10,
            Error::MissingIndex(_) => 11,
            Error::ReferenceBuildMismatch { .. } => 12,
            Error::ReferenceMd5Mismatch { .. } => 13,
            Error::TruncatedInput => 14,
            Error::EmptyRegion { .. } => 15,
            Error::InvalidConfig(_) => 17,
            Error::Io(_) => 20,
            Error::Internal(_) => 99,
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    /// Every `Error` variant paired with its expected exit code. Adding a new
    /// variant means adding one row here; `exit_codes_are_stable` and
    /// `display_is_single_line_and_typical_length` both iterate this list.
    fn all_variants_with_codes() -> Vec<(Error, i32)> {
        vec![
            (Error::InvalidAlignmentFormat("input.foo".into()), 10),
            (Error::MissingIndex(".bai".into()), 11),
            (
                Error::ReferenceBuildMismatch {
                    got: "GRCh37".into(),
                },
                12,
            ),
            (
                Error::ReferenceMd5Mismatch {
                    expected: "1d3085cf3e275ecb01ed2b30db8b6e2c".into(),
                    got: "9b4c5b8e1a3f8d7e6c5b4a3f2e1d0c9b".into(),
                },
                13,
            ),
            (Error::TruncatedInput, 14),
            (
                Error::EmptyRegion {
                    region: "chr17:15229777-15265079".into(),
                },
                15,
            ),
            (
                Error::InvalidConfig("coverage_floor_refuse must be positive".into()),
                17,
            ),
            (
                Error::Io(std::io::Error::new(std::io::ErrorKind::NotFound, "fixture")),
                20,
            ),
            (Error::Internal("unexpected".into()), 99),
        ]
    }

    /// Typical-case bound for `Display` length. The longest realistic case is
    /// `ReferenceMd5Mismatch` at ~133 chars; 200 leaves headroom for paths.
    const MAX_DISPLAY_LEN: usize = 200;

    #[test]
    fn exit_codes_are_stable() {
        for (err, expected) in all_variants_with_codes() {
            assert_eq!(err.exit_code(), expected, "{err:?}");
        }
    }

    #[test]
    fn display_is_single_line_and_typical_length() {
        for (err, _) in all_variants_with_codes() {
            let s = err.to_string();
            assert!(
                s.len() <= MAX_DISPLAY_LEN,
                "Display too long ({} chars, max {}): {}",
                s.len(),
                MAX_DISPLAY_LEN,
                s
            );
            assert!(!s.contains('\n'), "Display contains newline: {s}");
            assert_eq!(
                s,
                s.trim(),
                "Display has leading/trailing whitespace: {s:?}"
            );
        }
    }

    #[test]
    fn io_error_propagates_via_question_mark() {
        // Verifies `#[from] std::io::Error` is wired correctly: `?` must
        // convert into `Error::Io` automatically.
        fn fallible() -> Result<String> {
            let s = std::fs::read_to_string("/definitely/does/not/exist/nds-test")?;
            Ok(s)
        }
        let err = fallible().unwrap_err();
        assert!(matches!(err, Error::Io(_)));
        assert_eq!(err.exit_code(), 20);
    }
}
