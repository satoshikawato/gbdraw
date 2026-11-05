//! The SVG writer: the text rules of `gbdraw.render.group_content` and
//! `gbdraw.svg.path_data`.

use std::fmt::Write;

/// ElementTree's attribute escaping (`_escape_attrib`).
fn escape_attribute(out: &mut String, value: &str) {
    if !value.bytes().any(|b| matches!(b, b'&' | b'<' | b'>' | b'"' | b'\r' | b'\n' | b'\t')) {
        out.push_str(value);
        return;
    }
    for c in value.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\r' => out.push_str("&#13;"),
            '\n' => out.push_str("&#10;"),
            '\t' => out.push_str("&#09;"),
            c => out.push(c),
        }
    }
}

/// `<tag a="1" b="2" />`: attributes sorted by name, empty values left out.
pub fn empty_element(out: &mut String, tag: &str, attributes: &mut [(&str, &str)]) {
    attributes.sort_unstable_by(|a, b| a.0.cmp(b.0));
    out.push('<');
    out.push_str(tag);
    for (name, value) in attributes.iter() {
        if value.is_empty() {
            continue;
        }
        out.push(' ');
        out.push_str(name);
        out.push_str("=\"");
        escape_attribute(out, value);
        out.push('"');
    }
    out.push_str(" />");
}

/// A path number: 3 decimals in the shortest form, no trailing zeros, no `-0`.
pub fn path_number(out: &mut String, value: f64) {
    let start = out.len();
    write!(out, "{value:.3}").expect("a String accepts every write");
    let kept = out[start..].trim_end_matches('0').len();
    out.truncate(start + kept);
    if out.ends_with('.') {
        out.pop();
    }
    if &out[start..] == "-0" {
        out.truncate(start);
        out.push('0');
    }
}

/// One `x,y` pair of path data.
pub fn path_point(out: &mut String, x: f64, y: f64) {
    path_number(out, x);
    out.push(',');
    path_number(out, y);
}

/// Python's `format(value, f".{precision}g")` of a finite number.
pub fn format_g(out: &mut String, value: f64, precision: usize) {
    if value == 0.0 {
        out.push_str(if value.is_sign_negative() { "-0" } else { "0" });
        return;
    }
    let scientific = format!("{:.*e}", precision - 1, value);
    let (mantissa, exponent) = scientific.split_once('e').expect("`e` format writes an exponent");
    let exponent: i32 = exponent.parse().expect("`e` format writes an integer exponent");
    let trimmed = |text: &str, out: &mut String| {
        out.push_str(if text.contains('.') { text.trim_end_matches('0').trim_end_matches('.') } else { text });
    };
    if exponent < -4 || exponent >= precision as i32 {
        trimmed(mantissa, out);
        let sign = if exponent < 0 { '-' } else { '+' };
        write!(out, "e{sign}{:02}", exponent.abs()).expect("a String accepts every write");
    } else {
        let decimals = (precision as i32 - 1 - exponent) as usize;
        trimmed(&format!("{value:.decimals$}"), out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn number(value: f64) -> String {
        let mut out = String::new();
        path_number(&mut out, value);
        out
    }

    fn g(value: f64, precision: usize) -> String {
        let mut out = String::new();
        format_g(&mut out, value, precision);
        out
    }

    #[test]
    fn path_numbers_match_python() {
        for (value, text) in [(0.0, "0"), (-0.0, "0"), (-0.0001, "0"), (1.5, "1.5"), (2.0, "2"), (0.0625, "0.062"), (0.1875, "0.188"), (-12.5625, "-12.562"), (1000.0, "1000")] {
            assert_eq!(number(value), text, "{value}");
        }
    }

    #[test]
    fn g_format_matches_python() {
        for (value, precision, text) in [(0.0, 12, "0"), (-0.0, 12, "-0"), (1e-50, 12, "1e-50"), (3.4e-120, 12, "3.4e-120"), (64.871, 12, "64.871"), (1350.0, 12, "1350"), (1e12, 12, "1e+12"), (999999999999.0, 12, "999999999999"), (0.0001, 12, "0.0001"), (0.00001, 12, "1e-05"), (0.64871, 6, "0.64871"), (1.0, 6, "1"), (123456.7, 6, "123457"), (1234567.0, 6, "1.23457e+06")] {
            assert_eq!(g(value, precision), text, "{value}");
        }
    }

    #[test]
    fn elements_sort_escape_and_skip_empty_values() {
        let mut out = String::new();
        let mut attributes = [("id", "a&b"), ("d", "M 0,0"), ("fill", ""), ("data-x", "\"<\t>\"")];
        empty_element(&mut out, "path", &mut attributes);
        assert_eq!(out, "<path d=\"M 0,0\" data-x=\"&quot;&lt;&#09;&gt;&quot;\" id=\"a&amp;b\" />");
    }
}
