use std::ffi::OsString;
use std::path::{Path, PathBuf};

/// Length of a string that can represent any [`u32`].
pub(crate) const U32_STR_LEN: usize = {
    let mut n = u32::MAX;
    let mut digits = 0;
    while n > 0 {
        digits += 1;
        n /= 10;
    }
    digits
};

/// Returns `file_name` padded with 0s until reaching length `pad_len`.
///
/// Returns `None` if `file_name` is invalid (see [`Path::file_name()`]).
///
/// Extensions characters do not count towards `pad_len`.
/// ```
pub(crate) fn pad_file_name(file_name: &str, pad_len: usize) -> Option<OsString> {
    // This is only safe if we take str and convert to path ourselves
    let path = Path::new(file_name);
    let stem = path.file_stem()?.to_string_lossy();
    let mut padded: OsString = format!("{stem:0>pad_len$}").into();
    if let Some(ext) = path.extension() {
        padded.push(".");
        padded.push(ext);
    }
    Some(padded)
}

pub(crate) fn file_path(outdir: impl AsRef<Path>, subfolder: &str, ext: &str, time_step: u32) -> Option<PathBuf> {
    let padded = pad_file_name(
        &format!("{time_step}.{ext}"),
        U32_STR_LEN
    )?;
    Some(outdir.as_ref().join(subfolder).join(padded))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pad_file_name() {
        assert_eq!(pad_file_name("1.txt", 2), Some(OsString::from("01.txt")));
    }
}