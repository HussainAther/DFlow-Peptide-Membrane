import math

def effective_length_with_tilt(length, tilt_angle_deg):
    angle_rad = math.radians(tilt_angle_deg)
    return length * math.cos(angle_rad)

def check_tilt_fit(peptide_length, membrane_thickness, angle_deg, margin=1):
    eff_len = effective_length_with_tilt(peptide_length, angle_deg)
    return eff_len >= (membrane_thickness + margin)

