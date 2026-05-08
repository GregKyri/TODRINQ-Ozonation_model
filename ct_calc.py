def calculate_ct(watertemp):

    def interpolate(y1, y2, x1, x2, x):
        if x1 == x2:
            return y1
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)

    # Define CT tables
    virus_ct_table = {
        0.5: {1: 0.23, 5: 0.15, 10: 0.13, 15: 0.08, 20: 0.06, 25: 0.04},
        1.0: {1: 0.46, 5: 0.30, 10: 0.26, 15: 0.16, 20: 0.13, 25: 0.08},
        1.5: {1: 0.68, 5: 0.45, 10: 0.38, 15: 0.23, 20: 0.19, 25: 0.12},
        2.0: {1: 0.90, 5: 0.60, 10: 0.50, 15: 0.30, 20: 0.25, 25: 0.15},
        2.5: {1: 1.40, 5: 0.90, 10: 0.80, 15: 0.50, 20: 0.40, 25: 0.25},
        3.0: {1: 1.80, 5: 1.20, 10: 1.00, 15: 0.60, 20: 0.50, 25: 0.30},
    }

    giard_ct_table = {
        0.5: {1: 0.48, 5: 0.32, 10: 0.23, 15: 0.16, 20: 0.12, 25: 0.08},
        1.0: {1: 0.97, 5: 0.63, 10: 0.48, 15: 0.32, 20: 0.24, 25: 0.16},
        1.5: {1: 1.50, 5: 0.95, 10: 0.72, 15: 0.48, 20: 0.36, 25: 0.24},
        2.0: {1: 1.90, 5: 1.30, 10: 0.95, 15: 0.63, 20: 0.48, 25: 0.32},
        2.5: {1: 2.40, 5: 1.60, 10: 1.20, 15: 0.79, 20: 0.60, 25: 0.40},
        3.0: {1: 2.90, 5: 1.90, 10: 1.40, 15: 0.95, 20: 0.72, 25: 0.48},
    }

    ecoli_ct_table = {
        0.5: {1: 0.24, 5: 0.16, 10: 0.12, 15: 0.08, 20: 0.06, 25: 0.04},
        1.0: {1: 0.49, 5: 0.32, 10: 0.24, 15: 0.16, 20: 0.12, 25: 0.08},
        1.5: {1: 0.75, 5: 0.48, 10: 0.36, 15: 0.24, 20: 0.18, 25: 0.12},
        2.0: {1: 0.95, 5: 0.65, 10: 0.48, 15: 0.32, 20: 0.24, 25: 0.16},
        2.5: {1: 1.20, 5: 0.80, 10: 0.60, 15: 0.40, 20: 0.30, 25: 0.20},
        3.0: {1: 1.45, 5: 0.95, 10: 0.70, 15: 0.48, 20: 0.36, 25: 0.24},
    }

    temperatures = [1, 5, 10, 15, 20, 25]

    watertemp = max(1, min(25, watertemp))

    temp_lower = max(t for t in temperatures if t <= watertemp)
    temp_upper = min(t for t in temperatures if t >= watertemp)

    def interpolate_ct(ct_table, log_removal):

        available_logs = sorted(ct_table.keys())

        # Bound log removal to available range
        log_removal = max(min(log_removal, max(available_logs)), min(available_logs))

        log_lower = max(l for l in available_logs if l <= log_removal)
        log_upper = min(l for l in available_logs if l >= log_removal)

        # Interpolate temperature first
        ct_lower_temp_lowlog = interpolate(
            ct_table[log_lower][temp_lower],
            ct_table[log_lower][temp_upper],
            temp_lower,
            temp_upper,
            watertemp
        )

        ct_lower_temp_highlog = interpolate(
            ct_table[log_upper][temp_lower],
            ct_table[log_upper][temp_upper],
            temp_lower,
            temp_upper,
            watertemp
        )

        # Interpolate log removal
        return interpolate(
            ct_lower_temp_lowlog,
            ct_lower_temp_highlog,
            log_lower,
            log_upper,
            log_removal
        )

    # Desired log removals
    giardialog = 1.5
    viruslog = 1.8
    ecolilog = 2.1

    ct_giardia = interpolate_ct(giard_ct_table, giardialog)
    ct_virus = interpolate_ct(virus_ct_table, viruslog)
    ct_ecoli = interpolate_ct(ecoli_ct_table, ecolilog)

    return max(ct_giardia, ct_virus, ct_ecoli)