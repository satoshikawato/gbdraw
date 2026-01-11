
umamba activate gbdraw-0.5.0-test

IN=NC_012920
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_separate_strands -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --track_type middle --suppress_gc --suppress_skew --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_separate_strands -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --track_type tuckin --suppress_gc --suppress_skew --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_separate_strands -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --track_type spreadout --suppress_gc --suppress_skew --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type middle
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type tuckin
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type spreadout
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_separate_strands_inner_labels -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --track_type middle --suppress_gc --suppress_skew --separate_strands --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_separate_strands_inner_labels -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --track_type tuckin --suppress_gc --suppress_skew --separate_strands --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_separate_strands_inner_labels -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --track_type spreadout --suppress_gc --suppress_skew --separate_strands --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_inner_labels -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type middle --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_inner_labels -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type tuckin --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_inner_labels -f svg,png --block_stroke_width 2 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type spreadout --allow_inner_labels

IN=AP027280
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type tuckin
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type middle
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type spreadout
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type tuckin --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type middle --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type spreadout --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_separate_strands -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type tuckin --suppress_gc --suppress_skew --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_separate_strands -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle --suppress_gc --suppress_skew --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_separate_strands -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type spreadout --suppress_gc --suppress_skew --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_separate_strands_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type tuckin --suppress_gc --suppress_skew --separate_strands --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_separate_strands_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle --suppress_gc --suppress_skew --separate_strands --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_separate_strands_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type spreadout --suppress_gc --suppress_skew --separate_strands --allow_inner_labels

IN=LC738868
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_separate_strands -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type tuckin --suppress_gc --suppress_skew --separate_strands -t custom_color_table.tsv -d modified_default_colors.tsv
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_separate_strands -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle --suppress_gc --suppress_skew --separate_strands -t custom_color_table.tsv -d modified_default_colors.tsv
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_separate_strands -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type spreadout --suppress_gc --suppress_skew --separate_strands -t custom_color_table.tsv -d modified_default_colors.tsv
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type tuckin -t custom_color_table.tsv -d modified_default_colors.tsv
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type middle -t custom_color_table.tsv -d modified_default_colors.tsv
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type spreadout -t custom_color_table.tsv -d modified_default_colors.tsv
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_separate_strands_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type tuckin --suppress_gc --suppress_skew --separate_strands -t custom_color_table.tsv -d modified_default_colors.tsv --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_separate_strands_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle --suppress_gc --suppress_skew --separate_strands -t custom_color_table.tsv -d modified_default_colors.tsv --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_separate_strands_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type spreadout --suppress_gc --suppress_skew --separate_strands -t custom_color_table.tsv -d modified_default_colors.tsv --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type tuckin -t custom_color_table.tsv -d modified_default_colors.tsv --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type middle -t custom_color_table.tsv -d modified_default_colors.tsv --allow_inner_labels
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_inner_labels -f svg --block_stroke_width 1 --block_stroke_color gray --show_labels --suppress_gc --suppress_skew --track_type spreadout -t custom_color_table.tsv -d modified_default_colors.tsv --allow_inner_labels


IN=NC_001623
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle  --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type tuckin  --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type spreadout  --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle  
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type tuckin  
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type spreadout 

IN=NC_001623
gbdraw linear -i ${IN}.gb -o ${IN}_linear_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels  --separate_strands
gbdraw linear -i ${IN}.gb -o ${IN}_linear -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels


IN=AP027078
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle  --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type tuckin  --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type spreadout  --separate_strands
gbdraw circular --gbk ${IN}.gb -o ${IN}_middle -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle  
gbdraw circular --gbk ${IN}.gb -o ${IN}_tuckin -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type tuckin  
gbdraw circular --gbk ${IN}.gb -o ${IN}_spreadout -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type spreadout 


gbdraw circular --gbk NC_010162.gb  -f svg --show_labels --label_whitelist NC_010162.whitelist.tsv  --separate_strands -w 10000 -s 1000 -t NC_010162.feature-specific_table.tsv --outer_label_x_radius_offset 0.90 --outer_label_y_radius_offset 1.0 --palette edelweiss --species "<i>Sorangium cellulosum</i>" --strain "So ce56"