//! Supported-subtypes matrix: single source of truth.
//!
//! Every user-facing surface that advertises what neuropathy-dna-scanner
//! screens for — the text report, the web app, the README — reads from the
//! [`SUPPORTED`] constant via one of the three `render_*` functions.
//! Adding, removing, or reordering a row is a public-API change.

use std::fmt::Write;

/// Disclaimer text embedded verbatim in every generated report, in the
/// README, and in the web application's result panel.
pub const DISCLAIMER: &str = r#"Research-grade screening tool only.

This tool is NOT a medical device, NOT a diagnostic test, and NOT a substitute for clinical genetic testing. It currently screens only for the subtypes listed as Supported in the matrix above. It does NOT screen for any other form of CMT, and it may miss atypical or mosaic cases of the subtypes it does screen for.

A positive result requires clinical confirmation by MLPA, array-CGH, or sequencing at a certified diagnostic laboratory before any medical decision.

A negative result does NOT exclude CMT. It excludes only the subtypes marked Supported at the version of the tool you ran. Anyone with symptoms suggestive of inherited neuropathy (foot drop, high arches, distal weakness, hammer toes, family history) should see a neurologist or clinical geneticist regardless of this tool's output.

The author is not a clinician and provides no medical advice. Use at your own discretion. By using this tool you acknowledge that you understand these limitations."#;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatrixEntry {
    pub subtype: &'static str,
    pub gene_mechanism: &'static str,
    pub share_of_cmt: &'static str,
    pub status: SubtypeStatus,
    pub roadmap: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubtypeStatus {
    Supported,
    SupportedSideEffect,
    NotSupportedRoadmap,
    OutOfScope,
}

impl SubtypeStatus {
    /// Plain-text label for the "current" status column. Source of truth for
    /// the markdown and HTML variants. `NotSupportedRoadmap` and `OutOfScope`
    /// share a label here — the distinction surfaces in the roadmap column.
    fn label_plain(self) -> &'static str {
        match self {
            SubtypeStatus::Supported => "Supported",
            SubtypeStatus::SupportedSideEffect => "Supported (side-effect of depth-ratio method)",
            SubtypeStatus::NotSupportedRoadmap | SubtypeStatus::OutOfScope => "Not supported",
        }
    }

    /// Markdown-pipe variant of the status label. Supported statuses are
    /// wrapped in `**bold**`; "Not supported" is plain.
    pub fn current_label_md(self) -> String {
        match self {
            SubtypeStatus::Supported | SubtypeStatus::SupportedSideEffect => {
                format!("**{}**", self.label_plain())
            }
            SubtypeStatus::NotSupportedRoadmap | SubtypeStatus::OutOfScope => {
                self.label_plain().to_string()
            }
        }
    }
}

pub const SUPPORTED: &[MatrixEntry] = &[
    MatrixEntry {
        subtype: "CMT1A",
        gene_mechanism: "PMP22 1.4 Mb duplication (CNV)",
        share_of_cmt: "~50%",
        status: SubtypeStatus::Supported,
        roadmap: "--",
    },
    MatrixEntry {
        subtype: "HNPP",
        gene_mechanism: "PMP22 deletion at same locus (CNV)",
        share_of_cmt: "reciprocal of CMT1A",
        status: SubtypeStatus::SupportedSideEffect,
        roadmap: "--",
    },
    MatrixEntry {
        subtype: "PTLS",
        gene_mechanism: "RAI1 3.5 Mb duplication (17p11.2 CNV)",
        share_of_cmt: "not CMT; ~1:20,000",
        status: SubtypeStatus::Supported,
        roadmap: "--",
    },
    MatrixEntry {
        subtype: "SMS",
        gene_mechanism: "RAI1 3.5 Mb deletion (17p11.2 CNV)",
        share_of_cmt: "not CMT; ~1:15,000",
        status: SubtypeStatus::SupportedSideEffect,
        roadmap: "--",
    },
    MatrixEntry {
        subtype: "YUHAL",
        gene_mechanism: "PMP22+RAI1 ~6 Mb contiguous duplication",
        share_of_cmt: "very rare",
        status: SubtypeStatus::Supported,
        roadmap: "--",
    },
    MatrixEntry {
        subtype: "CMTX1",
        gene_mechanism: "GJB1 point mutations (X-linked)",
        share_of_cmt: "~10%",
        status: SubtypeStatus::NotSupportedRoadmap,
        roadmap: "v0.2",
    },
    MatrixEntry {
        subtype: "CMT1B",
        gene_mechanism: "MPZ point mutations",
        share_of_cmt: "~4%",
        status: SubtypeStatus::NotSupportedRoadmap,
        roadmap: "v0.3",
    },
    MatrixEntry {
        subtype: "CMT2A",
        gene_mechanism: "MFN2 point mutations",
        share_of_cmt: "~4%",
        status: SubtypeStatus::NotSupportedRoadmap,
        roadmap: "v0.3",
    },
    MatrixEntry {
        subtype: "CMT1E",
        gene_mechanism: "PMP22 point mutations (not the duplication)",
        share_of_cmt: "~1%",
        status: SubtypeStatus::NotSupportedRoadmap,
        roadmap: "v0.3",
    },
    MatrixEntry {
        subtype: "CMT4 family",
        gene_mechanism: "GDAP1, SH3TC2, NDRG1, ... (recessive)",
        share_of_cmt: "rare",
        status: SubtypeStatus::NotSupportedRoadmap,
        roadmap: "v0.3",
    },
    MatrixEntry {
        subtype: "CMT2 long tail",
        gene_mechanism: "NEFL, HSPB1, GARS1, RAB7A, ...",
        share_of_cmt: "rare",
        status: SubtypeStatus::NotSupportedRoadmap,
        roadmap: "v0.3",
    },
    MatrixEntry {
        subtype: "CMT-DI (dominant intermediate)",
        gene_mechanism: "requires phase information",
        share_of_cmt: "rare",
        status: SubtypeStatus::OutOfScope,
        roadmap: "Out of scope",
    },
    MatrixEntry {
        subtype: "CMT with repeat expansions",
        gene_mechanism: "requires ExpansionHunter-class tooling",
        share_of_cmt: "very rare",
        status: SubtypeStatus::OutOfScope,
        roadmap: "Out of scope",
    },
    MatrixEntry {
        subtype: "Mosaic CMT",
        gene_mechanism: "requires statistical model (parascopy class)",
        share_of_cmt: "rare",
        status: SubtypeStatus::OutOfScope,
        roadmap: "Out of scope",
    },
];

/// Render the matrix as a markdown-pipe table. Works in any terminal and any
/// markdown viewer.
pub fn render_text() -> String {
    let mut out = String::new();
    out.push_str("| CMT subtype | Gene / mechanism | Approx. share of CMT cases | v0.1 (current) | Roadmap |\n");
    out.push_str("|---|---|---|---|---|\n");
    for entry in SUPPORTED {
        writeln!(
            out,
            "| {} | {} | {} | {} | {} |",
            entry.subtype,
            entry.gene_mechanism,
            entry.share_of_cmt,
            entry.status.current_label_md(),
            entry.roadmap,
        )
        .unwrap();
    }
    out
}

/// Markdown table for the README. Named separately from [`render_text`] so
/// the two callers can diverge without breaking source compatibility.
pub fn render_markdown() -> String {
    render_text()
}

/// HTML `<table>` variant with a `class="cmt-matrix"` hook for CSS.
pub fn render_html() -> String {
    let mut out = String::new();
    out.push_str("<table class=\"cmt-matrix\">\n");
    out.push_str("  <thead>\n");
    out.push_str("    <tr>\n");
    out.push_str("      <th>CMT subtype</th>\n");
    out.push_str("      <th>Gene / mechanism</th>\n");
    out.push_str("      <th>Approx. share of CMT cases</th>\n");
    out.push_str("      <th>v0.1 (current)</th>\n");
    out.push_str("      <th>Roadmap</th>\n");
    out.push_str("    </tr>\n");
    out.push_str("  </thead>\n");
    out.push_str("  <tbody>\n");
    for entry in SUPPORTED {
        out.push_str("    <tr>\n");
        writeln!(out, "      <td>{}</td>", html_escape(entry.subtype)).unwrap();
        writeln!(out, "      <td>{}</td>", html_escape(entry.gene_mechanism)).unwrap();
        writeln!(out, "      <td>{}</td>", html_escape(entry.share_of_cmt)).unwrap();
        writeln!(out, "      <td>{}</td>", html_status_cell(entry.status)).unwrap();
        writeln!(out, "      <td>{}</td>", html_escape(entry.roadmap)).unwrap();
        out.push_str("    </tr>\n");
    }
    out.push_str("  </tbody>\n");
    out.push_str("</table>\n");
    out
}

/// HTML-escape the four characters that matter for cell content.
fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// HTML status cell. Supported statuses are wrapped in `<strong>`; plain
/// text for "Not supported".
fn html_status_cell(status: SubtypeStatus) -> String {
    match status {
        SubtypeStatus::Supported | SubtypeStatus::SupportedSideEffect => {
            format!("<strong>{}</strong>", status.label_plain())
        }
        SubtypeStatus::NotSupportedRoadmap | SubtypeStatus::OutOfScope => {
            status.label_plain().to_string()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn supported_has_expected_shape() {
        assert_eq!(SUPPORTED.len(), 14);
        let supported_count = SUPPORTED
            .iter()
            .filter(|e| {
                matches!(
                    e.status,
                    SubtypeStatus::Supported | SubtypeStatus::SupportedSideEffect
                )
            })
            .count();
        assert_eq!(supported_count, 5, "CMT1A, HNPP, PTLS, SMS, YUHAL");
        let subtypes: Vec<&str> = SUPPORTED.iter().map(|e| e.subtype).collect();
        assert_eq!(&subtypes[..5], &["CMT1A", "HNPP", "PTLS", "SMS", "YUHAL"]);
    }

    #[test]
    fn status_labels_are_stable() {
        assert_eq!(SubtypeStatus::Supported.label_plain(), "Supported");
        assert_eq!(
            SubtypeStatus::SupportedSideEffect.label_plain(),
            "Supported (side-effect of depth-ratio method)"
        );
        assert_eq!(
            SubtypeStatus::NotSupportedRoadmap.label_plain(),
            "Not supported"
        );
        assert_eq!(SubtypeStatus::OutOfScope.label_plain(), "Not supported");

        assert_eq!(SubtypeStatus::Supported.current_label_md(), "**Supported**");
        assert_eq!(
            SubtypeStatus::SupportedSideEffect.current_label_md(),
            "**Supported (side-effect of depth-ratio method)**"
        );
        assert_eq!(
            SubtypeStatus::NotSupportedRoadmap.current_label_md(),
            "Not supported"
        );
        assert_eq!(
            SubtypeStatus::OutOfScope.current_label_md(),
            "Not supported"
        );
    }

    const EXPECTED_TEXT: &str = "\
| CMT subtype | Gene / mechanism | Approx. share of CMT cases | v0.1 (current) | Roadmap |
|---|---|---|---|---|
| CMT1A | PMP22 1.4 Mb duplication (CNV) | ~50% | **Supported** | -- |
| HNPP | PMP22 deletion at same locus (CNV) | reciprocal of CMT1A | **Supported (side-effect of depth-ratio method)** | -- |
| PTLS | RAI1 3.5 Mb duplication (17p11.2 CNV) | not CMT; ~1:20,000 | **Supported** | -- |
| SMS | RAI1 3.5 Mb deletion (17p11.2 CNV) | not CMT; ~1:15,000 | **Supported (side-effect of depth-ratio method)** | -- |
| YUHAL | PMP22+RAI1 ~6 Mb contiguous duplication | very rare | **Supported** | -- |
| CMTX1 | GJB1 point mutations (X-linked) | ~10% | Not supported | v0.2 |
| CMT1B | MPZ point mutations | ~4% | Not supported | v0.3 |
| CMT2A | MFN2 point mutations | ~4% | Not supported | v0.3 |
| CMT1E | PMP22 point mutations (not the duplication) | ~1% | Not supported | v0.3 |
| CMT4 family | GDAP1, SH3TC2, NDRG1, ... (recessive) | rare | Not supported | v0.3 |
| CMT2 long tail | NEFL, HSPB1, GARS1, RAB7A, ... | rare | Not supported | v0.3 |
| CMT-DI (dominant intermediate) | requires phase information | rare | Not supported | Out of scope |
| CMT with repeat expansions | requires ExpansionHunter-class tooling | very rare | Not supported | Out of scope |
| Mosaic CMT | requires statistical model (parascopy class) | rare | Not supported | Out of scope |
";

    #[test]
    fn render_text_snapshot() {
        assert_eq!(render_text(), EXPECTED_TEXT);
    }

    #[test]
    fn render_markdown_snapshot() {
        // v0.1: render_markdown() is byte-identical to render_text().
        assert_eq!(render_markdown(), EXPECTED_TEXT);
    }

    const EXPECTED_HTML: &str = "\
<table class=\"cmt-matrix\">
  <thead>
    <tr>
      <th>CMT subtype</th>
      <th>Gene / mechanism</th>
      <th>Approx. share of CMT cases</th>
      <th>v0.1 (current)</th>
      <th>Roadmap</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>CMT1A</td>
      <td>PMP22 1.4 Mb duplication (CNV)</td>
      <td>~50%</td>
      <td><strong>Supported</strong></td>
      <td>--</td>
    </tr>
    <tr>
      <td>HNPP</td>
      <td>PMP22 deletion at same locus (CNV)</td>
      <td>reciprocal of CMT1A</td>
      <td><strong>Supported (side-effect of depth-ratio method)</strong></td>
      <td>--</td>
    </tr>
    <tr>
      <td>PTLS</td>
      <td>RAI1 3.5 Mb duplication (17p11.2 CNV)</td>
      <td>not CMT; ~1:20,000</td>
      <td><strong>Supported</strong></td>
      <td>--</td>
    </tr>
    <tr>
      <td>SMS</td>
      <td>RAI1 3.5 Mb deletion (17p11.2 CNV)</td>
      <td>not CMT; ~1:15,000</td>
      <td><strong>Supported (side-effect of depth-ratio method)</strong></td>
      <td>--</td>
    </tr>
    <tr>
      <td>YUHAL</td>
      <td>PMP22+RAI1 ~6 Mb contiguous duplication</td>
      <td>very rare</td>
      <td><strong>Supported</strong></td>
      <td>--</td>
    </tr>
    <tr>
      <td>CMTX1</td>
      <td>GJB1 point mutations (X-linked)</td>
      <td>~10%</td>
      <td>Not supported</td>
      <td>v0.2</td>
    </tr>
    <tr>
      <td>CMT1B</td>
      <td>MPZ point mutations</td>
      <td>~4%</td>
      <td>Not supported</td>
      <td>v0.3</td>
    </tr>
    <tr>
      <td>CMT2A</td>
      <td>MFN2 point mutations</td>
      <td>~4%</td>
      <td>Not supported</td>
      <td>v0.3</td>
    </tr>
    <tr>
      <td>CMT1E</td>
      <td>PMP22 point mutations (not the duplication)</td>
      <td>~1%</td>
      <td>Not supported</td>
      <td>v0.3</td>
    </tr>
    <tr>
      <td>CMT4 family</td>
      <td>GDAP1, SH3TC2, NDRG1, ... (recessive)</td>
      <td>rare</td>
      <td>Not supported</td>
      <td>v0.3</td>
    </tr>
    <tr>
      <td>CMT2 long tail</td>
      <td>NEFL, HSPB1, GARS1, RAB7A, ...</td>
      <td>rare</td>
      <td>Not supported</td>
      <td>v0.3</td>
    </tr>
    <tr>
      <td>CMT-DI (dominant intermediate)</td>
      <td>requires phase information</td>
      <td>rare</td>
      <td>Not supported</td>
      <td>Out of scope</td>
    </tr>
    <tr>
      <td>CMT with repeat expansions</td>
      <td>requires ExpansionHunter-class tooling</td>
      <td>very rare</td>
      <td>Not supported</td>
      <td>Out of scope</td>
    </tr>
    <tr>
      <td>Mosaic CMT</td>
      <td>requires statistical model (parascopy class)</td>
      <td>rare</td>
      <td>Not supported</td>
      <td>Out of scope</td>
    </tr>
  </tbody>
</table>
";

    #[test]
    fn render_html_snapshot() {
        assert_eq!(render_html(), EXPECTED_HTML);
    }

    #[test]
    fn html_escape_handles_special_chars() {
        assert_eq!(html_escape("plain"), "plain");
        assert_eq!(html_escape("a & b"), "a &amp; b");
        assert_eq!(html_escape("<tag>"), "&lt;tag&gt;");
        assert_eq!(html_escape("a&b<c>d"), "a&amp;b&lt;c&gt;d");
    }
}
