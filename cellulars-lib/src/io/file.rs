use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};

/// Length of a string that can represent any [`u32`].
pub const U32_STR_LEN: usize = {
    let mut n = u32::MAX;
    let mut digits = 0;
    while n > 0 {
        digits += 1;
        n /= 10;
    }
    digits
};

/// Returns `path` padded with 0s until reaching length `pad_len`.
///
/// Returns `None` if `path` points to ".." or "\".
///
/// Extensions characters do not count towards `pad_len`.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
/// use cellulars_lib::io::file::pad_file_path;
///
/// let p = pad_file_path("1.txt", 2).unwrap();
/// assert_eq!(p, PathBuf::from("01.txt"));
/// ```
pub fn pad_file_path(path: impl AsRef<Path>, pad_len: usize) -> Option<PathBuf> {
    let path = path.as_ref();
    let ext = path.extension().unwrap_or(OsStr::new(""));
    let file_name = path.file_name()?;
    let mut padded = OsString::new();
    dbg!(file_name.len() - ext.len() - 1);
    if pad_len > file_name.len() - ext.len() - 1 {
        for _ in 0..(pad_len - (file_name.len() - ext.len() - 1)) {
            padded.push("0");
        }
    }
    padded.push(file_name);
    Some(path.with_file_name(padded))
}