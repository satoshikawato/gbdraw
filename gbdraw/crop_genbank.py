#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
from io import StringIO  
import re  
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (
    SeqFeature, FeatureLocation, CompoundLocation,
    BeforePosition, AfterPosition
)

def _get_args():
    parser = argparse.ArgumentParser(description='Crop genbank file. ')
    parser.add_argument(
        '-i',
        '--input',
        help='Genbank/DDBJ flatfile (required)',
        type=str,
        required=True)
    parser.add_argument(
        "--output",
        "-o",
        "--out",
        "--output",
        metavar="FILE",
        help="output Genbank file",
        required=True)
    parser.add_argument(
        "-s",
        "--start",
        type=int,
        help="start position (1-based, inclusive)",
        required=True)
    parser.add_argument(
        "-e",
        "--end",
        type=int,
        help="end position (1-based, inclusive)",
        required=True)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


def gbk_to_seqrecord(in_gbk):
    """Load the first instance of GenBank as SerRecord
    Args:
    """
    records = SeqIO.parse(in_gbk, 'genbank')
    record = next(records)
    return record


def check_start_end_coords(record, start, end):
    """
    Converts 1-based inclusive start/end to 0-based exclusive slice indices.
    Also performs boundary checks.
    """
    record_len = len(record.seq)
    
    # Convert 1-based inclusive start to 0-based slice start
    start_0 = start - 1
    if start_0 < 0:
        start_0 = 0
    
    # 1-based inclusive end is the 0-based exclusive slice end
    end_0 = end
    if end_0 > record_len:
        end_0 = record_len

    if start_0 >= end_0 and (start != 1 or end != record_len): # Allow full sequence
        raise ValueError(f"Start position ({start}) must be less than end position ({end}).")
        
    return start_0, end_0

def _crop_and_shift_location(loc_before, loc_current, loc_next, crop_start_0, crop_end_0):
    """
    Takes a single FeatureLocation (loc) and crops/shifts it
    relative to the crop window (crop_start_0, crop_end_0).
    Returns a new FeatureLocation with < and > boundaries if partial.
    
    This version relies on Biopython's Position arithmetic, 
    which correctly preserves existing < and > markers when shifting.
    """
    new_seq_len = crop_end_0 - crop_start_0

    strand = loc_current.strand
    if strand == 1:
        if loc_before and loc_before.end < crop_start_0:
            final_start_pos = BeforePosition(max(0, loc_current.start - crop_start_0))
        elif loc_current.start < crop_start_0:
            final_start_pos = BeforePosition(max(0, loc_current.start - crop_start_0))
        else:
            final_start_pos = max(0, loc_current.start - crop_start_0)

        if loc_next and loc_next.start > crop_end_0:
            final_end_pos = AfterPosition(min(new_seq_len, loc_current.end - crop_start_0))
        elif loc_current.end > crop_end_0:
            final_end_pos = AfterPosition(min(new_seq_len, loc_current.end - crop_start_0))
        else:
            final_end_pos = min(new_seq_len, loc_current.end - crop_start_0)
        return FeatureLocation(
            final_start_pos,
            final_end_pos,
            strand=loc_current.strand
        )
    else:
        if loc_before and loc_before.start > crop_end_0:
            final_end_pos = AfterPosition(min(new_seq_len, loc_current.end - crop_start_0))
        else:
            final_end_pos = min(new_seq_len, loc_current.end - crop_start_0)

        if loc_next and loc_next.end < crop_start_0:
            final_start_pos = BeforePosition(max(0, loc_current.start - crop_start_0))
        else:
            final_start_pos = max(0, loc_current.start - crop_start_0)
    
        return FeatureLocation(
            final_start_pos,
            final_end_pos,
            strand=loc_current.strand
        )

def crop_and_shift_features(original_features, start_0, end_0):
    """
    Iterates over old features, crops/shifts them relative to the
    crop window (start_0, end_0), and returns a new list of features.
    """
    new_features_list = []
    for old_feature in original_features:

        if not (old_feature.location.start < end_0 and old_feature.location.end > start_0):
            continue

        old_loc = old_feature.location
        new_location = None 
        
        if isinstance(old_loc, CompoundLocation):
            new_parts = []
            num_parts = len(old_loc.parts)
            for n in range(num_parts):
                part = old_loc.parts[n]
                
                if not (part.start < end_0 and part.end > start_0):
                    continue
                
                part_before = old_loc.parts[n - 1] if n > 0 else None
                part_next = old_loc.parts[n + 1] if n < (num_parts - 1) else None
                
                new_part = _crop_and_shift_location(part_before, part, part_next, start_0, end_0)
                new_parts.append(new_part)
            
            if not new_parts:
                continue 
            
            if len(new_parts) == 1:
                new_location = new_parts[0] 
            else:
                new_location = CompoundLocation(new_parts, operator=old_loc.operator)
        
        else: 
            part_before = None
            part_next = None
            new_location = _crop_and_shift_location(part_before, old_loc, part_next, start_0, end_0)
        
        new_feature = SeqFeature(
            location=new_location,
            type=old_feature.type,
            qualifiers=old_feature.qualifiers, 
            id=old_feature.id
        )
        new_features_list.append(new_feature)
        
    return new_features_list


def write_cropped_genbank(record_to_write, out_filepath, original_accession, start_1based, end_1based):
    """
    Writes the SeqRecord to a file, modifying the ACCESSION line
    to include the cropped REGION.
    """
    temp_handle = StringIO()  
    SeqIO.write(record_to_write, temp_handle, "genbank")  
    content = temp_handle.getvalue()  
    
    # Add REGION to ACCESSION line
    content = re.sub(  
        rf'ACCESSION\s+{original_accession}',  
        f'ACCESSION   {original_accession} REGION: {start_1based}..{end_1based}',  
        content  
    )  

    with open(out_filepath, 'w') as f:  
        f.write(content)

def main():
    # 1. Setup
    args = _get_args()
    in_gbk = args.input
    out_gbk = args.output
    start_1 = args.start
    end_1 = args.end
    
    # 2. Load and Validate
    record = gbk_to_seqrecord(in_gbk)
    start_0, end_0 = check_start_end_coords(record, start_1, end_1)

    # 3. Crop Sequence and Create New Record
    new_seq = record.seq[start_0:end_0]
    new_record = SeqRecord(
        new_seq,
        id=record.id,
        name=record.name,
        description=record.description,
        dbxrefs=record.dbxrefs.copy(),
        annotations=record.annotations
    )
    
    # 4. Process Features (Extracted Function)
    new_record.features = crop_and_shift_features(record.features, start_0, end_0)

    # 5. Write Output (Extracted Function)
    original_accession = record.annotations['accessions'][0]
    write_cropped_genbank(new_record, out_gbk, original_accession, start_1, end_1)

if __name__ == "__main__":
    main()
    
