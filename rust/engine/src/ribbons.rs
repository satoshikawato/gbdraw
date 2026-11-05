//! The ribbons of one pairwise comparison: the children of
//! `gbdraw.render.groups.linear.pairwise_match.PairWiseMatchGroup` for matches
//! without collinearity, orthogroup or feature-attachment data.

use crate::svg::{empty_element, format_g, path_point};
use crate::Unsupported;
use serde_json::Value;
use std::fmt::Write;

/// The ten row columns, in the order of `_REQUIRED_ROW_ATTRIBUTES`.
const COLUMNS: [&str; 10] = [
    "qstart", "qend", "sstart", "send", "identity", "alignment_length", "evalue", "bitscore", "mismatches", "gap_opens",
];
const QSTART: usize = 0;
const QEND: usize = 1;
const SSTART: usize = 2;
const SEND: usize = 3;
const IDENTITY: usize = 4;

enum Column {
    Integers(Vec<i64>),
    Floats(Vec<f64>),
}

impl Column {
    fn number(&self, row: usize) -> f64 {
        match self {
            Column::Integers(values) => values[row] as f64,
            Column::Floats(values) => values[row],
        }
    }

    /// `_required_attribute_text` of one cell.
    fn write_text(&self, row: usize, out: &mut String) {
        out.clear();
        match self {
            Column::Integers(values) => write!(out, "{}", values[row]).expect("a String accepts every write"),
            Column::Floats(values) if values[row].is_finite() => format_g(out, values[row], 12),
            Column::Floats(_) => out.push(' '),
        }
    }
}

fn unsupported(reason: &str) -> Unsupported {
    Unsupported(reason.to_string())
}

fn read_columns(params: &Value, bytes: &[u8]) -> Result<(usize, Vec<Column>), Unsupported> {
    let rows = params["rows"].as_u64().ok_or_else(|| unsupported("no row count"))? as usize;
    if bytes.len() != rows * 8 * COLUMNS.len() {
        return Err(unsupported("column bytes do not match the row count"));
    }
    let mut columns = Vec::with_capacity(COLUMNS.len());
    for (name, data) in COLUMNS.iter().zip(bytes.chunks_exact((rows * 8).max(1))) {
        let words = data.chunks_exact(8).map(|word| word.try_into().expect("chunks of 8 bytes"));
        columns.push(match params["columns"][name].as_str() {
            Some("i") => Column::Integers(words.map(i64::from_le_bytes).collect()),
            Some("f") => Column::Floats(words.map(f64::from_le_bytes).collect()),
            _ => return Err(Unsupported(format!("column {name} is not numeric"))),
        });
    }
    Ok((rows, columns))
}

/// `_hex_to_rgb`: `#RGB` or `#RRGGBB`.
fn hex_color(text: &str) -> Result<[i64; 3], Unsupported> {
    let digits = text.trim().strip_prefix('#').ok_or_else(|| unsupported("not a hex color"))?;
    if !digits.bytes().all(|b| b.is_ascii_hexdigit()) || !matches!(digits.len(), 3 | 6) {
        return Err(unsupported("not a hex color"));
    }
    let channel = |index: usize| {
        let pair = if digits.len() == 3 { digits[index..=index].repeat(2) } else { digits[index * 2..index * 2 + 2].to_string() };
        i64::from_str_radix(&pair, 16).expect("two hex digits")
    };
    Ok([channel(0), channel(1), channel(2)])
}

fn number(params: &Value, pointer: &str) -> Result<f64, Unsupported> {
    params.pointer(pointer).and_then(Value::as_f64).ok_or_else(|| Unsupported(format!("no number at {pointer}")))
}

fn text<'a>(params: &'a Value, pointer: &str) -> Result<&'a str, Unsupported> {
    params.pointer(pointer).and_then(Value::as_str).ok_or_else(|| Unsupported(format!("no text at {pointer}")))
}

/// `_match_draw_order_key` without collinearity data: a stable sort, as `sorted`.
fn draw_order(rows: usize, columns: &[Column]) -> Vec<usize> {
    let keys: Vec<[f64; 5]> = (0..rows)
        .map(|row| {
            let at = |column: usize| columns[column].number(row);
            let query_span = (at(QEND) - at(QSTART)).abs();
            let subject_span = (at(SEND) - at(SSTART)).abs();
            [-(query_span * subject_span), -(query_span + subject_span), at(IDENTITY), at(QSTART), at(SSTART)]
        })
        .collect();
    let mut order: Vec<usize> = (0..rows).collect();
    order.sort_by(|a, b| keys[*a].partial_cmp(&keys[*b]).expect("finite keys"));
    order
}

/// The children text of a comparison group. `params` is the JSON the Python
/// group writes; `columns` holds `COLUMNS` in order, each `rows` little-endian
/// 8-byte values.
pub fn render(params: &str, columns: &[u8]) -> Result<String, Unsupported> {
    let params: Value = serde_json::from_str(params).map_err(|error| Unsupported(error.to_string()))?;
    let (rows, columns) = read_columns(&params, columns)?;
    for row in 0..rows {
        if [QSTART, QEND, SSTART, SEND, IDENTITY].iter().any(|column| !columns[*column].number(row).is_finite()) {
            return Err(unsupported("a position or identity is not finite"));
        }
    }

    let (query_x, query_scale, query_y) = (number(&params, "/query/x")?, number(&params, "/query/px_per_bp")?, number(&params, "/query/y")?);
    let (subject_x, subject_scale, subject_y) = (number(&params, "/subject/x")?, number(&params, "/subject/px_per_bp")?, number(&params, "/subject/y")?);
    let (minimum, maximum) = (number(&params, "/identity_domain/0")?, number(&params, "/identity_domain/1")?);
    let (low, high) = (hex_color(text(&params, "/min_color")?)?, hex_color(text(&params, "/max_color")?)?);
    let (opacity, stroke, stroke_width) = (text(&params, "/fill_opacity")?, text(&params, "/stroke")?, text(&params, "/stroke_width")?);
    let track = text(&params, "/track_id")?;
    let curve = params["style"].as_str() == Some("curve");
    let style = if curve { "curve" } else { "ribbon" };
    let tension = if curve { number(&params, "/curve_tension")? } else { 0.0 };
    let query_id = text(&params, "/endpoints/data-query-record-id")?;
    let query_index = text(&params, "/endpoints/data-query-record-index")?;
    let subject_id = text(&params, "/endpoints/data-subject-record-id")?;
    let subject_index = text(&params, "/endpoints/data-subject-record-index")?;

    let mut out = String::with_capacity(rows * 640);
    let (mut d, mut match_id, mut fill, mut factor_text) = (String::new(), String::new(), String::new(), String::new());
    let mut cells: [String; 10] = Default::default();
    for (index, &row) in draw_order(rows, &columns).iter().enumerate() {
        let identity = columns[IDENTITY].number(row);
        let factor = if minimum >= maximum {
            if identity >= maximum { 1.0 } else { 0.0 }
        } else {
            f64::max(0.0, f64::min(1.0, (identity - minimum) / (maximum - minimum)))
        };
        let channel = |i: usize| (low[i] as f64 + (high[i] - low[i]) as f64 * factor) as i64;
        fill.clear();
        write!(fill, "#{:02x}{:02x}{:02x}", channel(0), channel(1), channel(2)).expect("a String accepts every write");

        let on_query = |position: f64| query_x + position * query_scale;
        let on_subject = |position: f64| subject_x + position * subject_scale;
        let (q0, q1) = (on_query(columns[QSTART].number(row)), on_query(columns[QEND].number(row)));
        let (s0, s1) = (on_subject(columns[SSTART].number(row)), on_subject(columns[SEND].number(row)));
        d.clear();
        d.push_str("M ");
        path_point(&mut d, q0, query_y);
        if curve {
            let dy = subject_y - query_y;
            let (c1, c2) = (query_y + dy * tension, subject_y - dy * tension);
            d.push('C');
            path_point(&mut d, q0, c1);
            d.push(' ');
            path_point(&mut d, s0, c2);
            d.push(' ');
            path_point(&mut d, s0, subject_y);
            d.push_str(" L");
            path_point(&mut d, s1, subject_y);
            d.push('C');
            path_point(&mut d, s1, c2);
            d.push(' ');
            path_point(&mut d, q1, c1);
            d.push(' ');
            path_point(&mut d, q1, query_y);
        } else {
            d.push('L');
            path_point(&mut d, q1, query_y);
            d.push_str(" L");
            path_point(&mut d, s1, subject_y);
            d.push('L');
            path_point(&mut d, s0, subject_y);
        }
        d.push_str(" z");

        match_id.clear();
        write!(match_id, "{track}_match{}", index + 1).expect("a String accepts every write");
        factor_text.clear();
        format_g(&mut factor_text, factor, 6);
        for (cell, column) in cells.iter_mut().zip(&columns) {
            column.write_text(row, cell);
        }
        let mut attributes: [(&str, &str); 23] = [
            ("d", &d),
            ("fill", &fill),
            ("fill-opacity", opacity),
            ("stroke", stroke),
            ("stroke-width", stroke_width),
            ("data-pairwise-match-style", style),
            ("data-identity-factor", &factor_text),
            ("data-gbdraw-match-id", &match_id),
            ("data-match-kind", "pairwise"),
            ("data-query-record-index", query_index),
            ("data-subject-record-index", subject_index),
            ("data-query-record-id", query_id),
            ("data-subject-record-id", subject_id),
            ("data-qstart", &cells[0]),
            ("data-qend", &cells[1]),
            ("data-sstart", &cells[2]),
            ("data-send", &cells[3]),
            ("data-identity", &cells[4]),
            ("data-alignment-length", &cells[5]),
            ("data-evalue", &cells[6]),
            ("data-bitscore", &cells[7]),
            ("data-mismatches", &cells[8]),
            ("data-gap-opens", &cells[9]),
        ];
        empty_element(&mut out, "path", &mut attributes);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const PARAMS: &str = r##"{"rows":2,"columns":{"qstart":"i","qend":"i","sstart":"i","send":"i","identity":"f","alignment_length":"i","evalue":"f","bitscore":"f","mismatches":"i","gap_opens":"i"},
        "track_id":"comparison1","query":{"x":10.0,"px_per_bp":0.5,"y":0.0},"subject":{"x":20.0,"px_per_bp":0.25,"y":60.0},"identity_domain":[0.0,100.0],
        "min_color":"#FFE7E7","max_color":"#FF7272","fill_opacity":"1.0","stroke":"none","stroke_width":"0.0","style":"ribbon","curve_tension":0.5,
        "endpoints":{"data-query-record-index":"0","data-subject-record-index":"1","data-query-record-id":"A&B","data-subject-record-id":"C"}}"##;

    fn columns() -> Vec<u8> {
        let mut bytes = Vec::new();
        for row in [[1i64, 3], [11, 4], [2, 8], [12, 9]] {
            bytes.extend(row.iter().flat_map(|v| v.to_le_bytes()));
        }
        bytes.extend([90.5f64, 64.871].iter().flat_map(|v| v.to_le_bytes()));
        bytes.extend([10i64, 2].iter().flat_map(|v| v.to_le_bytes()));
        bytes.extend([1e-50f64, 0.0].iter().flat_map(|v| v.to_le_bytes()));
        bytes.extend([1350.0f64, 40.5].iter().flat_map(|v| v.to_le_bytes()));
        for row in [[1i64, 0], [0, 1]] {
            bytes.extend(row.iter().flat_map(|v| v.to_le_bytes()));
        }
        bytes
    }

    #[test]
    fn draws_the_larger_match_first_with_python_attribute_order() {
        let text = render(PARAMS, &columns()).unwrap();
        assert!(text.starts_with("<path d=\"M 10.5,0L15.5,0 L23,60L20.5,60 z\" data-alignment-length=\"10\" data-bitscore=\"1350\" data-evalue=\"1e-50\""), "{text}");
        assert!(text.contains("data-gbdraw-match-id=\"comparison1_match1\" data-identity=\"90.5\" data-identity-factor=\"0.905\""), "{text}");
        assert!(text.contains("data-query-record-id=\"A&amp;B\""), "{text}");
        assert_eq!(text.matches("<path").count(), 2);
    }

    #[test]
    fn refuses_what_python_must_report() {
        assert!(render(&PARAMS.replace("#FF7272", "red"), &columns()).is_err());
        assert!(render(PARAMS, &columns()[8..]).is_err());
        let mut bytes = columns();
        bytes[4 * 16..4 * 16 + 8].copy_from_slice(&f64::NAN.to_le_bytes());
        assert!(render(PARAMS, &bytes).is_err());
    }
}
