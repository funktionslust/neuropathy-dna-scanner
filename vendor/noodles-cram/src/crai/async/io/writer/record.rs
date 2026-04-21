use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::crai::Record;

pub(super) async fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const LINE_FEED: u8 = b'\n';
    const UNMAPPED: i32 = -1;

    let reference_sequence_id = if let Some(id) = record.reference_sequence_id() {
        id.to_string()
    } else {
        UNMAPPED.to_string()
    };

    let alignment_start = record
        .alignment_start()
        .map(usize::from)
        .unwrap_or_default();

    let s = format!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        reference_sequence_id,
        alignment_start,
        record.alignment_span(),
        record.offset(),
        record.landmark(),
        record.slice_length()
    );

    writer.write_all(s.as_bytes()).await?;
    writer.write_all(&[LINE_FEED]).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_core::Position;

        let record = Record::new(Some(0), Position::new(10946), 6765, 17711, 233, 317811);

        let mut buf = Vec::new();
        write_record(&mut buf, &record).await?;

        let expected = b"0\t10946\t6765\t17711\t233\t317811\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
