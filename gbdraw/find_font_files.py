#!/usr/bin/env python
# coding: utf-8

import logging
import subprocess
from fontTools.ttLib import TTFont

logger = logging.getLogger()


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
def get_text_bbox_size_pixels(font_path, text, font_size, dpi):
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
        advance_width, lsb = hmtx[glyph_index]
        rsb = (advance_width - (lsb + xmax - xmin))
        total_width += advance_width
        if car_count == 0:
            total_width -= lsb
        else:
            pass
        car_count +=1
    total_width -= rsb  
    units_per_em = font['head'].unitsPerEm
    pixel = dpi * font_size / 72
    text_width_pixels = (total_width/units_per_em) * pixel
    
    text_height_pixels = ((max(ymaxes) + abs(min(ymins)))/units_per_em) * pixel
    return text_width_pixels, text_height_pixels
