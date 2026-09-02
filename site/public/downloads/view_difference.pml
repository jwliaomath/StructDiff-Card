reinitialize
load reference.pdb, reference
load mobile_aligned_difference.pdb, mobile
hide everything
show cartoon, reference or mobile
color gray80, reference
spectrum b, blue_white_red, mobile and b >= 0, minimum=0, maximum=5
color violet, mobile and b < 0
set cartoon_transparency, 0.45, reference
bg_color white
orient
