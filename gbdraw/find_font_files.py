#!/usr/bin/env python
# coding: utf-8

import logging
import subprocess
from fontTools.ttLib import TTFont

# Logging setup
logger = logging.getLogger(__name__)

def find_font_paths(font_families):
    font_paths = []
    for font_family in font_families:
        try:
            fclist_out = subprocess.run(['fc-list', f':family={font_family}'], capture_output=True, text=True, check=True).stdout
            font_paths.extend(line.strip() for line in fclist_out.splitlines() if line.strip())
        except subprocess.CalledProcessError as e:
            message = f"Error querying font family '{font_family}': {e}"
            #logger.info(message)
    return font_paths

def process_font_line(line):
    parts = line.split(':')
    if len(parts) < 3:
        return None, None, []
    return parts[0].strip(), parts[1].strip(), parts[2].strip("style=").split(',')

def get_font_dict(font_families, styles):
    font_dict = {}
    font_file_paths = find_font_paths(font_families)

    for font_line in font_file_paths:
        file_path, font_family, available_styles = process_font_line(font_line)
        if not font_family:
            continue
        font_dict.setdefault(font_family, {})
        for style in styles:
            if style in available_styles:
                font_dict[font_family][style] = file_path

    # Print messages for missing fonts or styles
    for font in font_families:
        if font not in font_dict:
            message = f"Font '{font}' could not be found on the system."
            #logger.info(message)
        else:
            for style in styles:
                if style not in font_dict[font]:
                    message = f"Style '{style}' could not be found for the font '{font}'."
                    #logger.info(message)

    return font_dict

def _get_kerning_value(font, left_glyph_index, right_glyph_index):
    """
    Get kerning adjustment value for a glyph pair from kern or GPOS tables.
    
    Args:
        font: TTFont object
        left_glyph_index: Glyph index of the left glyph
        right_glyph_index: Glyph index of the right glyph
    
    Returns:
        int: Kerning adjustment value in font design units, or 0 if not found
    """
    if left_glyph_index is None or right_glyph_index is None:
        return 0
    
    try:
        # Try kern table first (simpler, faster, used in older fonts)
        if "kern" in font:
            kern_table = font["kern"]
            # kern table can have multiple subtables
            if hasattr(kern_table, "kernTables"):
                for subtable in kern_table.kernTables:
                    # Check if this subtable has horizontal kerning (bit 0 of coverage)
                    if subtable.coverage & 1:
                        # Look up the glyph pair
                        pair = (left_glyph_index, right_glyph_index)
                        if hasattr(subtable, "kernTable") and pair in subtable.kernTable:
                            return subtable.kernTable[pair]
    except Exception as e:
        logger.debug(f"Error reading kern table: {e}")
    
    try:
        # Try GPOS table (modern OpenType fonts)
        if "GPOS" in font:
            gpos_table = font["GPOS"]
            if hasattr(gpos_table, "table") and gpos_table.table:
                lookup_list = gpos_table.table.LookupList
                if lookup_list:
                    # Iterate through lookups to find pair adjustment
                    for lookup in lookup_list.Lookup:
                        if lookup.LookupType == 2:  # Pair adjustment lookup
                            for subtable in lookup.SubTable:
                                # Format 1: Pair adjustment positioning
                                if hasattr(subtable, "PairSets") and hasattr(subtable, "Coverage"):
                                    # Check if left glyph is in coverage
                                    coverage = subtable.Coverage
                                    if hasattr(coverage, "glyphs") and left_glyph_index in coverage.glyphs:
                                        # Find the index of left glyph in coverage
                                        try:
                                            left_idx = coverage.glyphs.index(left_glyph_index)
                                            if left_idx < len(subtable.PairSets):
                                                pair_set = subtable.PairSets[left_idx]
                                                if pair_set:
                                                    for pair_value_record in pair_set:
                                                        if hasattr(pair_value_record, "SecondGlyph"):
                                                            if pair_value_record.SecondGlyph == right_glyph_index:
                                                                if hasattr(pair_value_record, "Value1"):
                                                                    value_record = pair_value_record.Value1
                                                                    if hasattr(value_record, "XAdvance"):
                                                                        return value_record.XAdvance
                                        except (ValueError, IndexError, AttributeError):
                                            pass
                                    # Format 2: Class-based pair adjustment
                                    elif hasattr(subtable, "ClassDef1") and hasattr(subtable, "ClassDef2"):
                                        # Get class for left and right glyphs
                                        class_def1 = subtable.ClassDef1
                                        class_def2 = subtable.ClassDef2
                                        class1 = class_def1.get(left_glyph_index, 0) if hasattr(class_def1, "get") else 0
                                        class2 = class_def2.get(right_glyph_index, 0) if hasattr(class_def2, "get") else 0
                                        if hasattr(subtable, "Class1Record"):
                                            if class1 < len(subtable.Class1Record):
                                                class1_record = subtable.Class1Record[class1]
                                                if hasattr(class1_record, "Class2Record"):
                                                    if class2 < len(class1_record.Class2Record):
                                                        class2_record = class1_record.Class2Record[class2]
                                                        if hasattr(class2_record, "Value1"):
                                                            value_record = class2_record.Value1
                                                            if hasattr(value_record, "XAdvance"):
                                                                return value_record.XAdvance
    except Exception as e:
        logger.debug(f"Error reading GPOS table: {e}")
    
    return 0

def get_text_bbox_size_pixels(font_path, text, font_size, dpi):
    """
    This function directly parses the font file using the `fontTools.ttLib` library
    to read the metrics (dimensional information) of individual glyphs (characters).
    It calculates the width and height of the entire text string by considering
    each character's advance width, side bearings, and vertical extents (yMax, yMin).
    The width calculation accounts for kerning adjustments from the font's kern
    table (legacy fonts) or GPOS table (modern OpenType fonts).

    Args:
        font_path (str): Path to the font file (e.g., .ttf) to be used.
        text (str): The text string for which to calculate the size.
        font_size (int or float): The font size (in points).
        dpi (int): Dots Per Inch (DPI). Used to convert font units into pixels.

    Returns:
        tuple[float, float]:
            text_width_pixels (float): The calculated width of the text in pixels,
                including kerning adjustments.
            text_height_pixels (float): The calculated height of the text in pixels.
    """
    font = TTFont(font_path)
    hmtx = font['hmtx']
    cmap = font['cmap']
    glyph = font['glyf']
    ymax = font['head'].yMax
    ymin = font['head'].yMin
    t = cmap.getcmap(3, 1).cmap
    total_width = 0
    ymaxes = []
    ymins = []
    car_count = 0
    text = str(text)
    rsb_previous = 0
    previous_glyph_index = None
    for char in text:
        ymax = 0
        ymin = 0
        glyph_index = t[ord(char)]
        try: 
            ymax = glyph[glyph_index].yMax
        except AttributeError:
            ymax = 0
        try:             
            ymin = glyph[glyph_index].yMin
        except AttributeError:
            ymin = 0
        try:             
            xmin = glyph[glyph_index].xMin
        except AttributeError:
            xmin = 0
        try:             
            xmax = glyph[glyph_index].xMax
        except AttributeError:
            xmax = 0
        ymaxes.append(ymax)
        ymins.append(ymin)
        
        # Apply kerning adjustment between consecutive glyphs
        if previous_glyph_index is not None:
            kerning_value = _get_kerning_value(font, previous_glyph_index, glyph_index)
            total_width += kerning_value
        
        advance_width, lsb = hmtx[glyph_index]
        if xmax == 0:
            advance_width -= rsb_previous
        if lsb > 0:
            total_width -= lsb
        else:
            total_width += lsb
        rsb = (xmax -advance_width)
        total_width += advance_width
        if rsb > 0:
            total_width -= rsb
        else:
            total_width += rsb
        if car_count == 0:
            total_width += lsb
        else:
            pass
        rsb_previous = rsb
        previous_glyph_index = glyph_index
        car_count +=1

    units_per_em = font['head'].unitsPerEm
    pixel = dpi * font_size / 72
    text_width_pixels = (total_width/units_per_em) * pixel
    
    text_height_pixels = ((max(ymaxes) + abs(min(ymins)))/units_per_em) * pixel
    return text_width_pixels, text_height_pixels
