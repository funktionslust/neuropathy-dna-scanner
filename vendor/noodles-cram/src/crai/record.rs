//! CRAM index record and fields.

use noodles_core::Position;

/// A CRAM index record.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Record {
    reference_sequence_id: Option<usize>,
    alignment_start: Option<Position>,
    alignment_span: usize,
    offset: u64,
    landmark: u64,
    slice_length: u64,
}

impl Record {
    /// Creates a CRAM index record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    /// ```
    pub fn new(
        reference_sequence_id: Option<usize>,
        alignment_start: Option<Position>,
        alignment_span: usize,
        offset: u64,
        landmark: u64,
        slice_length: u64,
    ) -> Self {
        Self {
            reference_sequence_id,
            alignment_start,
            alignment_span,
            offset,
            landmark,
            slice_length,
        }
    }

    /// Returns the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.reference_sequence_id(), Some(0));
    /// ```
    pub fn reference_sequence_id(&self) -> Option<usize> {
        self.reference_sequence_id
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.alignment_start(), Position::new(10946));
    /// ```
    pub fn alignment_start(&self) -> Option<Position> {
        self.alignment_start
    }

    /// Returns the alignment span.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.alignment_span(), 6765);
    /// ```
    pub fn alignment_span(&self) -> usize {
        self.alignment_span
    }

    /// Returns the offset of the container from the start of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.offset(), 17711);
    /// ```
    pub fn offset(&self) -> u64 {
        self.offset
    }

    /// Returns the offset of the slice from the start of the container.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.landmark(), 233);
    /// ```
    pub fn landmark(&self) -> u64 {
        self.landmark
    }

    /// Returns the size of the slice in bytes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.slice_length(), 317811);
    /// ```
    pub fn slice_length(&self) -> u64 {
        self.slice_length
    }
}
